# The Exchange–Correlation Energy in CP2K (GPW and GAPW)

### A guided tour of how `E_xc` is computed, from the density on the grid down to the Fortran

This document explains, step by step, how CP2K evaluates the **exchange–correlation
(XC) energy contribution** to the total Kohn–Sham energy in the two standard
Gaussian-based schemes:

- **GPW** — Gaussian and Plane Waves (`dft_control%qs_control%gapw = .FALSE.`)
- **GAPW** — Gaussian Augmented Plane Waves (`gapw = .TRUE.`), and the closely related
  `gapw_xc` variant.

It is written to be read top-to-bottom and points at the exact source lines that do
the work. The companion document
[`nuclear_contributions_core_ae.md`](nuclear_contributions_core_ae.md) covers the
nuclear/Hartree side; this one is purely about `E_xc`.

The energy slot we are filling is `energy%exc` (and, for GAPW, the additional
`energy%exc1`). Both are defined in
[`src/qs_energy_types.F`](../src/qs_energy_types.F#L41) and folded into
`energy%total` in [`src/qs_ks_methods.F`](../src/qs_ks_methods.F#L1648-L1655).

---

## 1. The physical problem

In Kohn–Sham DFT the XC energy is a functional of the electron density (and, for GGA,
its gradient; for meta-GGA, the kinetic-energy density $\tau$):

$$
E_{xc}[\rho] \;=\; \int \varepsilon_{xc}\big(\rho(\mathbf r),\,\nabla\rho(\mathbf r),\,\tau(\mathbf r)\big)\, \mathrm d\mathbf r .
$$

Here $\varepsilon_{xc}$ is the **XC energy density** (energy per unit volume). The
whole job of the code is to:

1. build $\rho$ (and $\nabla\rho$, $\tau$ as needed) on a numerical grid,
2. evaluate $\varepsilon_{xc}$ pointwise (this is the "functional"),
3. **integrate** $\varepsilon_{xc}$ over space to get the scalar $E_{xc}$,
4. and (for the SCF) also evaluate the **potential** $v_{xc} = \delta E_{xc}/\delta\rho$,
   which is integrated against the basis functions to build the KS matrix.

The grid used differs between GPW and GAPW. That is the whole story:

- **GPW** does the integral on a single, uniform **plane-wave (real-space) grid**.
- **GAPW** does it on the plane-wave grid for a *smooth* density **plus** a set of
  per-atom **spherical (Lebedev × radial) grids** for the sharp part near each nucleus.

---

## 2. The two ingredients: the functional evaluator and the integrator

Before tracing GPW/GAPW separately, it helps to see the two reusable pieces both
schemes share.

### 2.1 The `xc_rho_set` / `xc_derivative_set` machinery

Everything goes through
[`xc_rho_set_and_dset_create`](../src/xc/xc.F#L839) in
[`src/xc/xc.F`](../src/xc/xc.F). It takes the density (`rho_r` in real space, `rho_g`
in G-space for gradients, `tau` for meta-GGA) and:

- builds an `xc_rho_set` holding $\rho$, $|\nabla\rho|$, $\tau$, etc. on the grid;
- builds an `xc_derivative_set` (`deriv_set`) — a container of "derivatives" of the
  functional with respect to its arguments. The **0th-order derivative is the energy
  density** $\varepsilon_{xc}$ itself; the 1st-order derivatives are the pieces of the
  potential.

The actual functional (PBE, LYP, libxc, …) is evaluated by
[`xc_functionals_eval`](../src/xc/xc_derivatives.F) (called from
[`vxc_of_r_new`](../src/xc/xc_atom.F#L113) for the atomic path, and inside
`xc_rho_set_and_dset_create` for the grid path). After evaluation, the energy density
lives in the order-0 derivative, retrieved with the **empty descriptor**
`[INTEGER::]`.

### 2.2 Two entry points: energy-only vs. energy+potential

In [`src/xc/xc.F`](../src/xc/xc.F) there are two public routines, selected by the
`just_energy` flag upstream:

| Routine | Purpose | Returns |
|---|---|---|
| [`xc_exc_calc`](../src/xc/xc.F#L1459) | energy only | scalar `exc` |
| [`xc_vxc_pw_create1`](../src/xc/xc.F#L149) → [`xc_vxc_pw_create`](../src/xc/xc.F#L1144) | energy **and** potential | `exc` + `vxc_rho`, `vxc_tau` |

Both must give the **identical** `exc`; the comment at
[`xc.F:1402`](../src/xc/xc.F#L1402) ("this has to be kept consistent with
`xc_exc_calc`") flags exactly this requirement.

---

## 3. GPW: the plane-wave grid path

### 3.1 Where it is launched

The driver is `qs_ks_build_kohn_sham_matrix` in
[`src/qs_ks_methods.F`](../src/qs_ks_methods.F). For a normal (non-GAPW_XC) run the XC
call is at [`qs_ks_methods.F:792`](../src/qs_ks_methods.F#L792):

```fortran
CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho_struct, xc_section=xc_section, &
                   vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, exc=energy%exc, &
                   just_energy=just_energy_xc, ...)
```

So **`energy%exc` is filled directly by `qs_vxc_create`**. `rho_struct` here is the
total valence density represented on the plane-wave grid.

### 3.2 Inside `qs_vxc_create`

[`src/qs_vxc.F`](../src/qs_vxc.F#L95), `qs_vxc_create`:

1. **Fetches the density grids** ([`qs_vxc.F:158`](../src/qs_vxc.F#L158)) — `rho_r`
   (real space), `rho_g` (G-space, used so GGAs can take the gradient analytically by
   multiplying by `i·G`), and `tau_r` for meta-GGA.
2. **Handles grid choice** ([`qs_vxc.F:226-257`](../src/qs_vxc.F#L226)): the XC grid
   (`xc_pw_pool`) may be finer than the density grid (`auxbas_pw_pool`). If so the
   density is transferred to the finer XC grid first (`uf_grid` branch). This is the
   `XC_GRID` keyword in the input.
3. **Adds NLCC core density** if present ([`qs_vxc.F:260`](../src/qs_vxc.F#L260)) — the
   nonlinear core correction density participates in `E_xc`.
4. **Calls the evaluator** ([`qs_vxc.F:272-284`](../src/qs_vxc.F#L272)):
   - if `just_energy`: `exc = xc_exc_calc(...)`,
   - else: `xc_vxc_pw_create1(..., exc=exc, ...)`.
5. **SIC / scaling** ([`qs_vxc.F:326`](../src/qs_vxc.F#L326)): `exc = exc*my_scaling`
   applies any self-interaction-correction scaling. For a plain functional
   `my_scaling = 1`.

The returned `exc` becomes `energy%exc`.

### 3.3 The actual integral — `xc_exc_calc`

[`src/xc/xc.F:1459`](../src/xc/xc.F#L1459). This is the cleanest place to see the
energy being formed:

```fortran
CALL xc_rho_set_and_dset_create(..., deriv_order=0, ..., calc_potential=.FALSE.)
deriv => xc_dset_get_derivative(deriv_set, [INTEGER::])      ! the 0-th deriv = e_0
CALL xc_derivative_get(deriv, deriv_data=e_0)               ! e_0(i,j,k) = eps_xc

CALL smooth_cutoff(pot=e_0, rho=rho_set%rho, ...)           ! low-density tapering

exc = accurate_sum(e_0)*rho_r(1)%pw_grid%dvol               ! E_xc = sum(eps) * dV
IF (... PW_MODE_DISTRIBUTED) CALL ...%group%sum(exc)        ! MPI reduction
```

Step by step ([`xc.F:1481-1503`](../src/xc/xc.F#L1481)):

1. **Evaluate the functional to order 0** → `e_0` is the array of $\varepsilon_{xc}$
   values, one per grid point.
2. **`smooth_cutoff`** ([`xc.F:913`](../src/xc/xc.F#L913)) multiplies $\varepsilon_{xc}$
   by a smooth taper that goes to 0 below `DENSITY_CUTOFF` so that noise in the
   near-vacuum tail of $\rho$ does not contaminate the integral.
3. **The integral is a plain Riemann sum**: $E_{xc} = \sum_{\text{points}}
   \varepsilon_{xc}(\mathbf r) \,\Delta V$, where `dvol` is the volume per grid point.
   `accurate_sum` is a compensated (Kahan-style) summation for reproducibility.
4. **MPI reduction**: when the grid is distributed across ranks, the partial sums are
   added together.

That scalar is `E_xc`.

### 3.4 The energy+potential variant — `xc_vxc_pw_create`

When the potential is also needed ([`xc.F:1144`](../src/xc/xc.F#L1144)), the functional
is evaluated to **order 1** (`deriv_order=1`,
[`xc.F:1201`](../src/xc/xc.F#L1201)). The energy is then formed at
[`xc.F:1401-1410`](../src/xc/xc.F#L1401):

```fortran
! 0-deriv -> value of exc
CALL xc_dset_recover_pw(deriv_set, [INTEGER::], v_drho_r, pw_grid)  ! eps_xc as a pw
CALL smooth_cutoff(pot=v_drho_r%array, rho=rho, ...)               ! same taper
exc = pw_integrate_function(v_drho_r)                               ! integrate
```

`pw_integrate_function` is the same `sum × dvol` (+ MPI) operation as in
`xc_exc_calc`, just wrapped for a plane-wave object — which is why the two routines are
required to agree. In the same call the **potential** $v_{xc}$ is recovered into
`vxc_rho` (the 1st-order derivatives, [`xc.F:1227-1233`](../src/xc/xc.F#L1227)); for
GGAs the gradient pieces are handled by `xc_pw_divergence`, and the meta-GGA $\tau$
potential into `vxc_tau` ([`xc.F:1426`](../src/xc/xc.F#L1426)).

### 3.5 What happens to the potential (not the energy, but for completeness)

`vxc_rho` flows back up to `qs_ks_build_kohn_sham_matrix` and is later integrated
against basis-function pairs inside `sum_up_and_integrate`
([`qs_ks_methods.F:1411`](../src/qs_ks_methods.F#L1411) →
[`qs_ks_utils.F`](../src/qs_ks_utils.F)), which calls `integrate_v_rspace` to produce
the XC block of the KS matrix. The energy `energy%exc`, however, was already finalized
in step 3.3/3.4 — it is **not** recomputed as `Tr[P · V_xc]`. (Note the contrast with
the Hartree term, which *is* a trace; XC is a direct grid integral of the energy
density.)

---

## 4. GAPW: grid (soft) + atomic (hard − soft) decomposition

GAPW splits every density into a **soft** part (smooth, representable on the plane-wave
grid) and a **hard** part (the sharp atomic region), using the PAW-like identity

$$
E_{xc}[\rho] \;=\; \underbrace{E_{xc}[\tilde\rho]}_{\text{grid, soft}}
\;+\; \sum_{a}\Big(\underbrace{E_{xc}[\rho_a^{1}]}_{\text{hard, atom }a}
- \underbrace{E_{xc}[\tilde\rho_a^{1}]}_{\text{soft, atom }a}\Big).
$$

In code this maps cleanly onto two energy slots:

- **`energy%exc`** = $E_{xc}[\tilde\rho]$ — the *soft* density on the plane-wave grid,
  computed by **exactly the same `qs_vxc_create` path as GPW** (Section 3).
- **`energy%exc1`** = $\sum_a (E_{xc}[\rho_a^1] - E_{xc}[\tilde\rho_a^1])$ — the
  per-atom hard-minus-soft correction, computed by
  [`calculate_vxc_atom`](../src/qs_vxc_atom.F#L83).

### 4.1 Preparing the GAPW densities

Before any XC work, [`qs_ks_methods.F:489-498`](../src/qs_ks_methods.F#L489) calls
[`prepare_gapw_den`](../src/qs_gapw_densities.F): this projects the density onto each
atom and builds the radial expansions of the hard density $\rho_a^1$ and soft density
$\tilde\rho_a^1$ stored in the `rho_atom_set`.

### 4.2 The soft (grid) part — same as GPW

The call at [`qs_ks_methods.F:792`](../src/qs_ks_methods.F#L792) computes
`energy%exc` from `rho_struct`, which in GAPW mode holds the **soft** total density.
Everything in Section 3 applies unchanged — GAPW reuses the GPW grid integral for the
soft term.

### 4.3 The atomic (hard − soft) part — `calculate_vxc_atom`

Immediately after, [`qs_ks_methods.F:809-810`](../src/qs_ks_methods.F#L809):

```fortran
IF (gapw .OR. gapw_xc) THEN
   CALL calculate_vxc_atom(qs_env, just_energy_xc, energy%exc1, xc_section_external=xc_section)
END IF
```

Inside [`calculate_vxc_atom`](../src/qs_vxc_atom.F#L83):

1. **Loop over atom kinds, then atoms** ([`qs_vxc_atom.F:203,268`](../src/qs_vxc_atom.F#L203));
   skip atoms that are not PAW atoms ([`:209`](../src/qs_vxc_atom.F#L209)).
2. **Build the angular density on the per-atom grid**
   ([`qs_vxc_atom.F:297-317`](../src/qs_vxc_atom.F#L297)): the radial coefficients
   `rho_rad_h`/`rho_rad_s` are expanded over the Lebedev angular grid (`na` angular ×
   `nr` radial points) into `rho_h` (hard) and `rho_s` (soft) via `calc_rho_angular`,
   then loaded into `rho_set_h` / `rho_set_s` by `fill_rho_set`.
3. **Evaluate the functional and integrate, twice** — once for hard, once for soft
   ([`qs_vxc_atom.F:322-335`](../src/qs_vxc_atom.F#L322)):

   ```fortran
   CALL vxc_of_r_new(..., rho_set_h, ..., exc_h, vxc_h, ...)   ! hard atomic density
   rho_atom%exc_h = rho_atom%exc_h + exc_h
   CALL vxc_of_r_new(..., rho_set_s, ..., exc_s, vxc_s, ...)   ! soft atomic density
   rho_atom%exc_s = rho_atom%exc_s + exc_s
   ```

4. **Form the correction** ([`qs_vxc_atom.F:339`](../src/qs_vxc_atom.F#L339)):

   ```fortran
   exc1 = exc1 + rho_atom%exc_h - rho_atom%exc_s
   ```

5. **MPI reduction over atoms** ([`qs_vxc_atom.F:370`](../src/qs_vxc_atom.F#L370)):
   `CALL para_env%sum(exc1)`.

### 4.4 The per-atom integral — `vxc_of_r_new`

[`src/xc/xc_atom.F:61`](../src/xc/xc_atom.F#L61). This is the GAPW analogue of
`xc_exc_calc`. The energy is formed at
[`xc_atom.F:123-134`](../src/xc/xc_atom.F#L123):

```fortran
CALL xc_functionals_eval(xc_fun_section, lsd, rho_set, deriv_set, deriv_order)
deriv_att => xc_dset_get_derivative(deriv_set, [INTEGER::])   ! the 0-th deriv
CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)      ! eps_xc(ia,ir)
exc = 0.0_dp
DO ir = 1, nr           ! radial points
   DO ia = 1, na        ! angular (Lebedev) points
      exc = exc + deriv_data(ia, ir, 1)*w(ia, ir)             ! eps_xc * weight
   END DO
END DO
```

So the atomic integral is again **$\sum \varepsilon_{xc} \times \text{weight}$**, but
the quadrature weights `w(ia,ir)` are the **spherical grid weights** (radial weight ×
Lebedev angular weight, from `grid_atom%weight`) rather than a uniform `dvol`. Same
functional `xc_functionals_eval`, same order-0 descriptor `[INTEGER::]`, different
quadrature.

### 4.5 Assembling the GAPW total

The two slots are summed independently into the total at
[`qs_ks_methods.F:1648-1649`](../src/qs_ks_methods.F#L1648):

```fortran
energy%total = ... + energy%exc + energy%exc1 + ...
```

and printed separately for transparency
([`qs_ks_utils.F:1136-1139`](../src/qs_ks_utils.F#L1136)):

```
GAPW| Exc from hard and soft atomic rho1:   <energy%exc1 + energy%exc1_aux_fit>
```

So the full GAPW XC energy is `energy%exc` (grid/soft) **+** `energy%exc1`
(atomic hard − soft). When comparing against GPW you must add both.

### 4.6 `gapw_xc`

`gapw_xc` ([`qs_ks_methods.F:771,809`](../src/qs_ks_methods.F#L771)) is a lighter
variant: it applies the GAPW soft/hard split **only to the XC term** (and not to the
Hartree term). The XC machinery is identical to Section 4.2–4.4 — `energy%exc` from the
soft grid density, `energy%exc1` from `calculate_vxc_atom`. The difference is purely in
how the rest of the density (Hartree) is treated, not in the XC integral itself.

---

## 5. The potential side (for the KS matrix)

For completeness — `vxc_of_r_new` also returns the potential `vxc_h`/`vxc_s` (and
gradient/tau pieces) when `energy_only = .FALSE.`
([`xc_atom.F:136`](../src/xc/xc_atom.F#L136)). These are turned into matrix elements by
`gaVxcgb_noGC` / `gaVxcgb_GC`
([`qs_vxc_atom.F:349-353`](../src/qs_vxc_atom.F#L349)), which integrate the potential
against pairs of GAPW_1C basis functions on the atomic grid and store the result in the
`ga_Vlocal_gb` blocks of `rho_atom`. Those blocks are later added to the KS matrix.
The **energy**, again, is the direct grid integral of $\varepsilon_{xc}$, not a
`Tr[P·V]`.

---

## 6. Summary table

| Quantity | GPW | GAPW |
|---|---|---|
| Soft/grid energy `energy%exc` | $\sum_{\mathbf r}\varepsilon_{xc}[\rho]\,\Delta V$ on the plane-wave grid | same, but on the **soft** density $\tilde\rho$ |
| Atomic energy `energy%exc1` | — (always 0) | $\sum_a(E_{xc}[\rho_a^1]-E_{xc}[\tilde\rho_a^1])$ on per-atom spherical grids |
| Functional evaluator | `xc_functionals_eval` via `xc_rho_set_and_dset_create` | `xc_functionals_eval` via `vxc_of_r_new` |
| Energy density slot | order-0 derivative `[INTEGER::]` (`e_0`) | order-0 derivative `[INTEGER::]` (`deriv_data`) |
| Quadrature weight | uniform `pw_grid%dvol` | `grid_atom%weight(ia,ir)` (radial × Lebedev) |
| Integrator | `accurate_sum` / `pw_integrate_function` | explicit double loop over `(ia,ir)` |
| Low-density taper | `smooth_cutoff` (`DENSITY_CUTOFF`) | functional-internal cutoffs (`density_cut`/`gradient_cut`) |
| Total XC energy | `energy%exc` | `energy%exc + energy%exc1` |

---

## 7. File / routine index

| File | Routine | Role |
|---|---|---|
| [`src/qs_ks_methods.F`](../src/qs_ks_methods.F#L792) | `qs_ks_build_kohn_sham_matrix` | launches XC, holds `energy%exc`/`energy%exc1`, sums total |
| [`src/qs_vxc.F`](../src/qs_vxc.F#L95) | `qs_vxc_create` | top-level grid XC: density prep, grid choice, NLCC, SIC scaling |
| [`src/xc/xc.F`](../src/xc/xc.F#L1459) | `xc_exc_calc` | grid **energy-only** integral |
| [`src/xc/xc.F`](../src/xc/xc.F#L1144) | `xc_vxc_pw_create` | grid **energy + potential** |
| [`src/xc/xc.F`](../src/xc/xc.F#L839) | `xc_rho_set_and_dset_create` | builds density arrays + functional derivative container |
| [`src/xc/xc.F`](../src/xc/xc.F#L913) | `smooth_cutoff` | low-density tapering of $\varepsilon_{xc}$ |
| [`src/qs_vxc_atom.F`](../src/qs_vxc_atom.F#L83) | `calculate_vxc_atom` | GAPW atomic hard−soft loop, fills `energy%exc1` |
| [`src/xc/xc_atom.F`](../src/xc/xc_atom.F#L61) | `vxc_of_r_new` | per-atom spherical-grid energy + potential |
| [`src/qs_gapw_densities.F`](../src/qs_gapw_densities.F) | `prepare_gapw_den` | builds hard/soft atomic radial densities |
| [`src/qs_energy_types.F`](../src/qs_energy_types.F#L41) | — | declares `exc`, `exc1`, `exc1_aux_fit` |

---

### One-paragraph takeaway

In both GPW and GAPW the XC energy is a **numerical quadrature of the XC energy
density**: evaluate the functional to order 0 to get $\varepsilon_{xc}$ at every grid
point, then sum $\varepsilon_{xc}\times(\text{quadrature weight})$ and reduce over MPI.
GPW does this once on the uniform plane-wave grid (`energy%exc`). GAPW does it three
times — once on the plane-wave grid for the soft density (`energy%exc`) and, per atom,
once for the hard and once for the soft atomic density, storing their difference in
`energy%exc1`. The total XC energy is `energy%exc` for GPW and
`energy%exc + energy%exc1` for GAPW. The potential is derived from the order-1
derivatives in the same evaluation but does **not** enter the energy via a `Tr[P·V]` —
the energy is the direct integral of the energy density.
