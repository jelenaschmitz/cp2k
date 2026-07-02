# Exchange–Correlation in RIDFT

### How RIDFT computes E_xc, and why it reuses the GPW/GAPW machinery unchanged

This document records the design decision and the one-line code change that make the
exchange–correlation (XC) treatment of **RIDFT** explicit. It is the RIDFT-specific
companion to [`xc_energy_gpw_gapw.md`](xc_energy_gpw_gapw.md), which describes the XC
energy in GPW and GAPW in detail. Read that first; this file only covers what is (and
is not) different for RIDFT.

The short version: **RIDFT does not have its own XC code, and it needs none.** It
evaluates XC with the *same* routine, on the *same* grid, as GPW, and it already builds
every grid quantity that routine needs. The only thing RIDFT changes in the Kohn–Sham
build is the Hartree/nuclear electrostatics (replaced by an RI integral plus the analytic
`build_core_erf_ridft` term). XC is left completely untouched. **The conclusion of this
work was therefore that no code change is required** — only this documentation.

---

## 1. The question this answers

RIDFT's purpose is to avoid the plane-wave Poisson solve for the Hartree potential and
the grid collocation of the sharp nuclear Gaussians, replacing them with

- an **RI integral** for the electron–electron Hartree term, and
- the analytic **`build_core_erf_ridft`** term for the electron–nucleus and
  nucleus–nucleus erf electrostatics (see
  [`nuclear_contributions_core_ae.md`](nuclear_contributions_core_ae.md)).

A natural worry is: *if RIDFT zeroes the grid Hartree, does it also lose XC?* and
*does RIDFT need its own XC grid?* The answer to both is no, and this document explains
why, so nobody re-implements something that already works.

---

## 2. What was already true (verified in code)

RIDFT already shares the GPW XC path end to end. No guard excludes RIDFT from it:

1. **The XC energy and potential are computed by the shared routine.**
   In `qs_ks_build_kohn_sham_matrix` ([`src/qs_ks_methods.F`](../src/qs_ks_methods.F#L792)),
   `qs_vxc_create` is called unconditionally (it is *not* wrapped in `IF (.NOT. ridft)`):

   ```fortran
   CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho_struct, xc_section=xc_section, &
                      vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, exc=energy%exc, &
                      edisp=edisp, dispersion_env=qs_env%dispersion_env, &
                      just_energy=just_energy_xc)
   ```

   This fills `energy%exc` exactly as in GPW (the grid integral of the XC energy
   density — see [`xc_energy_gpw_gapw.md`](xc_energy_gpw_gapw.md) §3).

2. **The density grid that XC needs already exists for RIDFT — verified.**
   `ridft` appears nowhere in the density-collocation path
   (`qs_rho_methods.F`, `qs_collocate_density.F`, `qs_scf.F`). In
   `qs_rho_update_rho` ([`src/qs_rho_methods.F:441`](../src/qs_rho_methods.F#L441)),
   RIDFT is none of harris / semi-empirical / DFTB / xTB / LRIGPW / RIGPW, so it takes
   the final `ELSE` branch → `qs_rho_update_rho_low`. There
   ([`qs_rho_methods.F:542-550`](../src/qs_rho_methods.F#L542)) `calculate_rho_elec` is
   called with `soft_valid = gapw = .FALSE.`, collocating the full valence density onto
   `rho_r` **and** `rho_g` on the auxbas grid, followed by
   `rho_r_valid = .TRUE., rho_g_valid = .TRUE.`. So RIDFT has exactly the grid quantities
   `qs_vxc_create` consumes: `rho_r` (real space), `rho_g` (G-space, for GGA gradients),
   and the `auxbas_pw_pool` / `xc_pw_pool` of the standard `pw_env`. `rho_struct` is the
   full valence density ([`qs_ks_methods.F:774`](../src/qs_ks_methods.F#L774)). Nothing is
   missing and nothing needs to be created.

3. **The XC potential is integrated into the KS matrix.**
   The only thing RIDFT zeroes is the *Hartree* real-space potential, in
   `sum_up_and_integrate` ([`src/qs_ks_utils.F:1325-1328`](../src/qs_ks_utils.F#L1325)):

   ```fortran
   IF (dft_control%qs_control%ridft) THEN
      CALL pw_zero(v_rspace)          ! v_rspace = v_hartree_rspace  (NOT the XC potential)
   END IF
   ```

   `v_rspace` here is `v_hartree_rspace`, retrieved a few lines above. The XC potential
   lives in a *separate* array, `v_rspace_new`, and is still integrated into the KS
   matrix at [`qs_ks_utils.F:1486`](../src/qs_ks_utils.F#L1486). So zeroing the Hartree
   grid potential prevents double counting of the Hartree term (RIDFT adds its RI
   Hartree separately) **without** affecting XC.

4. **The GAPW one-centre (atomic) XC correction is not added.**
   `calculate_vxc_atom` (which fills `energy%exc1`) is gated by
   `IF (gapw .OR. gapw_xc)` ([`qs_ks_methods.F:809`](../src/qs_ks_methods.F#L809)).
   RIDFT is a distinct `method_id` (`do_method_ridft`), so `gapw` and `gapw_xc` are
   both `.FALSE.` — the atomic term is correctly skipped. RIDFT's XC is the soft
   cell-grid term only, which is the whole XC for a GPW-style (non-augmented) density.

The reason this had "never worked" in testing is mundane: every RIDFT test input uses
`&XC_FUNCTIONAL NONE` (e.g.
`tests/.../17_GTH_WITH_RIDFT_Test/h2_rt.inp`), so `qs_vxc_create` returns early with
`exc = 0` and the XC path is never exercised. Switch the functional on and the existing
path computes XC correctly.

---

## 3. The decision

Because the XC computation is already correct and shared, the requirements were:

- **Use exactly the same XC as GPW/GAPW — do not change the XC calculation itself.**
- **Reuse the existing auxbas/xc PW grid (the grid around the input cell).** Do not
  build a new grid type.
- **XC is the soft cell-grid term only** (no GAPW atomic partitioning for RIDFT).
- **Out-of-cell point charges are handled elsewhere** — only through the RI `v_H` and
  the core (`build_core_erf_ridft`), never on the grid. That is a separate, already
  tracked task (the GTH / true-point-charge work), not part of the XC change.

So there is **no new XC code and no new grid**. The investigation's outcome is that the
existing code already satisfies every requirement.

---

## 4. The change: none (and why a guard was tried and then removed)

**Final state: no source change.** `qs_ks_methods.F` and the XC path are left exactly as
upstream.

During the work a single behaviour-neutral guard was briefly added at the XC entry to
"make the intent explicit":

```fortran
IF (dft_control%qs_control%ridft) THEN
   CPASSERT(.NOT. (gapw .OR. gapw_xc))   ! always true: method_id is exclusive
END IF
```

It was **removed** because it is pointless and slightly misleading:

- The assertion can never fail — `method_id` is mutually exclusive, so `ridft` already
  implies `.NOT. (gapw .OR. gapw_xc)`. The `IF` does no work.
- Worse, an `IF (ridft)` sitting in the XC region *implies* RIDFT needs special XC
  handling, when the whole point (verified in §2) is that it does not. Leaving it in
  would invite a future reader to add RIDFT-specific XC logic that is not wanted.

The faithful implementation of "use exactly the same XC as GPW/GAPW, do not change the
XC calculation" is therefore the empty diff. RIDFT energies are, by construction,
produced by the unchanged shared `qs_vxc_create` on the reused auxbas/xc grid.

---

## 5. What was explicitly NOT done (and why)

- **No GAPW atomic XC for RIDFT.** That would need the full GAPW atomic environment
  (`rho_atom_set`, `init_gapw_basis_set`, `init_rho0`, set up in
  [`qs_environment.F:1281`](../src/qs_environment.F#L1281) only for
  gapw/gapw_xc) plus `GAPW_1C` basis sets in every input. Out of scope; RIDFT uses the
  soft cell-grid term only.
- **No change to `calc_rho_tot_gspace` / `calculate_rho_core`.** Keeping the nuclear
  Gaussian core off the grid (so out-of-cell point charges never get collocated) is a
  separate concern handled by the RI `v_H` + `build_core_erf_ridft` path. It is tracked
  with the GTH / point-charge work, not here, to keep this change purely about XC and
  strictly behaviour-neutral.
- **No new XC grid object.** Reusing the existing auxbas/xc PW grid was the chosen
  approach.

---

## 6. How to verify

1. **Regression (XC off).** Because the source is unchanged, the existing RIDFT test
   (`17_GTH_WITH_RIDFT_Test`, `&XC_FUNCTIONAL NONE`) is bit-identical by construction.
2. **XC on.** Set a real functional (e.g. `&XC_FUNCTIONAL PBE`) in a RIDFT input and a
   matching GPW input on the same cell/cutoff/basis. RIDFT's `energy%exc` must equal the
   GPW `energy%exc` for the same density (the routine and grid are identical). Any
   remaining total-energy difference is the Hartree/nuclear RI fit error already
   documented for RIDFT, not an XC difference. This is the real test of the claim in §2 —
   that RIDFT's grid quantities (`rho_r`, `rho_g`) feed `qs_vxc_create` correctly.

---

### One-paragraph takeaway

RIDFT computes exchange–correlation by calling the same `qs_vxc_create` on the same
cell grid as GPW, producing `energy%exc` and the XC KS contribution with no RIDFT-specific
XC code. It already builds the grid quantities that routine needs — `rho_r` and `rho_g`
are collocated and marked valid in `qs_rho_update_rho_low`, on the auxbas grid. RIDFT only
replaces the Hartree/nuclear electrostatics; the `pw_zero` in `sum_up_and_integrate`
zeroes the Hartree grid potential, not the XC potential. The correct implementation is
therefore the **empty diff**: a guard was tried and removed because it asserts something
always true and wrongly implies RIDFT needs special XC handling. No new grid is created —
the existing auxbas/xc grid around the input cell is reused — and out-of-cell point
charges are left to the RI `v_H` + `build_core_erf_ridft` path.
