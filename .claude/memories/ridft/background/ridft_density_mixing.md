# RIDFT and SCF Density Mixing

**Status:** developer note for the `ri_dft` branch.
**Scope:** why the RIDFT SCF oscillates under Broyden/Pulay (g-space) mixing but converges
under `DIRECT_P_MIXING`, how this differs from GPW/GAPW, and the guard that enforces a
compatible mixer.

---

## 1. TL;DR

* GPW/GAPW build the Hartree potential from the **grid density** `rho_g`/`rho_r`.
  Under g-space mixing (Broyden, Pulay, …) it is exactly `rho_g`/`rho_r` that get mixed,
  so their Hartree feedback is damped and the SCF converges.
* RIDFT builds the Hartree potential from the **AO density matrix** `rho_ao(1,1)`
  (it needs a matrix for the RI 3-centre contraction) **and zeroes the grid Hartree**.
* The two mixing channels in CP2K are **mutually exclusive**:
  * `DIRECT_P_MIXING` mixes the **density matrix** `rho_ao`; the grid is *not* mixed.
  * `BROYDEN_MIXING`/`PULAY_MIXING`/… mix the **grid density** `rho_g`; the density
    matrix `rho_ao` is left as the raw, **un-mixed** output density `P_out`.
* Therefore RIDFT only sees a *mixed* density when `DIRECT_P_MIXING` is used. With a
  g-space mixer the quantity RIDFT depends on (`rho_ao`) is never mixed, the RI-Hartree
  feedback is undamped, and the SCF falls into a stable **period-2 limit cycle**.
* **"Can't RIDFT just use the same `rho` as GAPW?"** No — GAPW does not use a density
  *matrix* for its Hartree at all; it uses the mixed grid density. A *matrix*-based
  method has no mixed-matrix counterpart under g-space mixing. The only mixed density
  matrix that exists is the one produced by `DIRECT_P_MIXING`.
* **Fix applied:** a guard in `qs_ks_build_kohn_sham_matrix` aborts with a clear message
  if a g-space mixer is selected together with RIDFT (see §7).

---

## 2. The three representations of the density in CP2K

`qs_rho_type` carries the electron density in three coupled forms:

| Symbol | Object | Used by |
|---|---|---|
| `rho_ao` | DBCSR **density matrix** `P` in the AO basis | RI-Hartree (RIDFT), HFX, GAPW 1-centre terms, energy traces |
| `rho_r`  | density on the **real-space grid** | XC functional, integration of `V(r)` |
| `rho_g`  | density on the **reciprocal-space grid** | Poisson solve → grid Hartree (GPW/GAPW) |

`rho_r`/`rho_g` are obtained from `rho_ao` by `qs_rho_update_rho` (collocation onto the
grid). They are *derived* quantities — except that **density mixing can be applied to
either the matrix or the grid representation, and CP2K does exactly one of the two.**

---

## 3. The SCF density flow (diagonalisation + mixing)

One inner SCF iteration in `qs_scf.F` (`scf_env_do_scf` loop):

```
build KS  ──►  qs_ks_update_qs_env ──► qs_ks_build_kohn_sham_matrix   (RIDFT block lives here)
diagonalise ─► qs_scf_new_mos      (KS → MOs → new output density P_out, written into rho_ao)
mix (1)   ──►  qs_scf_density_mixing                          (src/qs_scf_loop_utils.F:483)
copy (2)  ──►  rho_ao_kp ← scf_env%p_mix_new                  (src/qs_scf.F:587-595)
mix (3)   ──►  qs_scf_rho_update(..., mix_rho = method>=gspace_mixing_nr)
                                                              (src/qs_scf.F:597, src/qs_scf_loop_utils.F:654)
```

### Step (1): `qs_scf_density_mixing` — sets `p_mix_new` (src/qs_scf_loop_utils.F:495-510)

```fortran
SELECT CASE (scf_env%mixing_method)
CASE (direct_mixing_nr)                       ! DIRECT_P_MIXING
   CALL scf_env_density_mixing(scf_env%p_mix_new, ... rho_ao_kp ...)   ! MIX THE MATRIX
CASE (gspace_mixing_nr, pulay_mixing_nr, broyden_mixing_nr, multisecant_mixing_nr)
   CALL self_consistency_check(rho_ao_kp, scf_env%p_delta, ..., scf_env%p_mix_new, ...)
   ! p_mix_new := P_out (raw), p_delta := P_out - P_in.  NO matrix mixing here.
CASE (no_mixing_nr)
END SELECT
```

* `DIRECT_P_MIXING`  ⇒  `p_mix_new` = **mixed** density matrix.
* g-space family      ⇒  `p_mix_new` = **un-mixed** output density matrix `P_out`.

### Step (2): copy back into the live density matrix (src/qs_scf.F:587-595)

```fortran
IF (scf_env%mixing_method > 0) THEN
   DO ic = ... ; DO ispin = ...
      CALL dbcsr_copy(rho_ao_kp(ispin, ic)%matrix, scf_env%p_mix_new(ispin, ic)%matrix, ...)
   END DO ; END DO
END IF
```

So after step (2), `rho_ao` holds:
* the **mixed matrix** for `DIRECT_P_MIXING`, or
* the **un-mixed `P_out`** for any g-space mixer.

### Step (3): grid recompute + optional grid mixing (src/qs_scf_loop_utils.F:654-675)

```fortran
SUBROUTINE qs_scf_rho_update(rho, qs_env, scf_env, ks_env, mix_rho)
   CALL qs_rho_update_rho(rho, qs_env=qs_env)         ! rho_r/rho_g recomputed FROM rho_ao
   IF (mix_rho) THEN                                  ! mix_rho = (method >= gspace_mixing_nr)
      CALL gspace_mixing(qs_env, scf_env%mixing_method, scf_env%mixing_store, rho, ...)
   END IF                                             ! mixes rho_g/rho_r ONLY
```

The header comment states the exclusivity explicitly:

> `! ** Density mixing through density matrix or on the reciprocal space grid (exclusive)`

Mixing-method constants (`src/qs_density_mixing_types.F:45-46`):

```
no_mixing_nr = 0,  direct_mixing_nr = 1,  gspace_mixing_nr = 2,
pulay_mixing_nr = 3,  broyden_mixing_nr = 4,  multisecant_mixing_nr = 5
```

### Net result per iteration

| mixing method | `rho_ao` (matrix) | `rho_g`/`rho_r` (grid) |
|---|---|---|
| `DIRECT_P_MIXING` (1) | **mixed** | recomputed from the mixed matrix → mixed, consistent |
| `BROYDEN`/`PULAY`/… (≥2) | **un-mixed `P_out`** | recomputed from `P_out`, then **mixed on the grid** |

This table is the whole story. The mixed density lives **in the matrix** for direct
mixing and **on the grid** for g-space mixing — never in both.

---

## 4. What GPW/GAPW use vs. what RIDFT uses

### GPW / GAPW (reference, converges)

The Hartree potential is obtained from the **total charge density on the reciprocal
grid** via the Poisson solver (`calc_rho_tot_gspace` + `pw_poisson_solve` in
`qs_ks_methods.F`), i.e. from `rho_g`. Under a g-space mixer, `rho_g` is precisely the
mixed quantity (step 3 above). XC likewise uses the mixed `rho_r`. The KS matrix is then
formed by integrating `V(r)` against the basis. **GAPW never needs a *mixed density
matrix*; its feedback density is the mixed grid density.** That is why GAPW (dir 12)
converges in 21 Broyden steps with this very basis.

### RIDFT (this branch)

RIDFT replaces the grid Poisson Hartree with a direct RI integral that is contracted
with the **AO density matrix**:

* `rho_ao_ridft` is sourced from `rho_ao(i,1)%matrix` — `src/qs_ks_methods.F:961`.
* the RI e–e Hartree `V_ee = J[rho_ao_ridft]` is built and added to the KS matrix.
* `build_core_erf_ridft` is passed `rho_ao` for the nuclear-attraction energy trace —
  `src/qs_ks_methods.F:1376`.
* **the grid Hartree is explicitly zeroed** so it cannot double count —
  `src/qs_ks_utils.F:1327` (`CALL pw_zero(v_rspace)` inside `sum_up_and_integrate`).

Because the grid Hartree is zeroed, **any mixing applied to `rho_g`/`rho_r` has no effect
on the RIDFT Hartree feedback.** RIDFT's feedback depends solely on the AO density matrix
`rho_ao`.

### Answer to "should that not be the correct `rho` matrix?"

`rho_ao(1,1)` is the *correct object* (it is the AO density matrix), but under a g-space
mixer it holds the **un-mixed** value `P_out`. GAPW avoids the problem not by using a
"more correct" matrix, but by **not using a matrix at all** for its Hartree — it uses the
mixed grid density. There is no mixed-matrix equivalent of the Broyden-mixed `rho_g`
(grid mixing destroys the one-to-one map to a single density matrix). Hence the only
self-consistent mixed density matrix available to a matrix-based method is the one from
`DIRECT_P_MIXING`.

---

## 5. Why `DIRECT_P_MIXING` works but `BROYDEN_MIXING` does not

Let `F` be the SCF map `P_in → P_out` (build KS, diagonalise, get new density). At a
solution `P* = F(P*)`. Convergence of plain iteration requires the dominant eigenvalue
`λ` of the Jacobian `dF/dP` to satisfy `|λ| < 1`. For systems with a strong Hartree
response (here: a tiny H2 in a large, near-linearly-dependent aug-cc-pVTZ basis) `λ` is
large and **negative** — successive densities over-/under-shoot, i.e. charge sloshing.
Damping (mixing) is what tames it: linear mixing with parameter `α` replaces `λ` by
`1 − α + α·λ`, which is `< 1` in magnitude for suitable `α`.

* **`DIRECT_P_MIXING`** applies that damping to the **density matrix** `rho_ao`
  (step 1, `direct_mixing_nr`). Since `rho_ao` is exactly what RIDFT reads, the RI-Hartree
  feedback is damped → the iteration contracts → **converges** (dir 15: 11 steps).

* **`BROYDEN_MIXING`** applies its (quasi-Newton) damping to the **grid density** `rho_g`
  (step 3). But RIDFT zeroes the grid Hartree and feeds on `rho_ao`, which under Broyden
  is the **raw `P_out`** (step 1 g-space branch). So RIDFT's feedback loop runs with
  **no damping at all** — effectively `α = 1` on the matrix. With `λ < −1` this is a
  stable period-2 limit cycle: `P_A → P_B → P_A → …`, neither collapsing to `P*`.
  Broyden's sophistication is wasted because it is acting on a quantity (`rho_g`) that is
  decoupled from RIDFT's actual feedback.

It is **not** that Broyden is "worse" than direct mixing in general — it is normally far
better. It simply mixes the *wrong representation* for a matrix-based, grid-Hartree-zeroed
method like RIDFT.

### Empirical confirmation (H2, aug-cc-pVTZ ORB, aug-cc-pVTZ-RIFIT, α = 0.4)

| run | method | mixing | result |
|---|---|---|---|
| dir 12 | GAPW | BROYDEN | converges, 21 steps |
| dir 14 | RIDFT | BROYDEN | **period-2 limit cycle**, never converges |
| dir 15 | RIDFT | DIRECT_P_MIXING | converges, 11 steps |

In dir 14 the Hartree energy alternates between two stable values every step:

```
energy%hartree:  2.5258  ↔  1.2435        (e_ee: 0.878 ↔ 2.953 ; e_nuc_erf: -1.449 ↔ -3.769)
```

Both branches stabilise (the step-to-step change shrinks to ~1e-5) — the hallmark of a
limit cycle rather than slow convergence. In dir 15 the same system converges to a single
density, and the RIDFT Hartree matches GAPW to ~2.7 mHa (RI-V fitting error; compare
RIDFT `energy%hartree` against GAPW `energy%hartree` **plus** `GAPW| local Eh = 1 center
integrals`).

---

## 6. Diagram

```
                         ┌───────────────────────────── SCF iteration ─────────────────────────────┐
  rho_ao (P_in) ─► build KS ─► diagonalise ─► P_out ─► [mix] ─► rho_ao(next) ─► recompute grid ─► [grid mix]
                      │                                  │                          │                   │
        GAPW Hartree ─┘ uses rho_g (grid) ◄──────────────┼──────────────────────────┘  mixed here for  │
                                                          │                             BROYDEN/PULAY ◄──┘
        RIDFT Hartree ── uses rho_ao (matrix) ◄───────────┘  mixed here ONLY for DIRECT_P_MIXING
        RIDFT grid Hartree = 0  (pw_zero, qs_ks_utils.F:1327)
```

* `DIRECT_P_MIXING`: the "[mix]" on the matrix is active → RIDFT sees a mixed `rho_ao`. ✅
* `BROYDEN`/`PULAY`: only the "[grid mix]" is active → RIDFT's `rho_ao` is `P_out`,
  un-mixed, and the grid mixing is invisible to RIDFT (grid Hartree zeroed). ❌

---

## 7. The fix (guard) — applied

`src/qs_ks_methods.F`, top of the RIDFT block in `qs_ks_build_kohn_sham_matrix`:

```fortran
CALL get_qs_env(qs_env, ..., scf_env=scf_env_ridft)

! RIDFT builds the Hartree potential directly from the density MATRIX rho_ao(1,1).
! Only matrix-level mixing (DIRECT_P_MIXING) produces a mixed density matrix; the
! g-space mixers (Broyden/Pulay) mix the grid density rho_g and leave rho_ao as the
! un-mixed output density, so the RIDFT Hartree feedback would be undamped and the
! SCF would charge-slosh into a period-2 limit cycle. Require a matrix mixer here.
IF (ASSOCIATED(scf_env_ridft)) THEN
   IF (scf_env_ridft%mixing_method >= gspace_mixing_nr) THEN
      CPABORT("RIDFT builds the Hartree from the density matrix and requires "// &
              "density-matrix mixing: set &MIXING METHOD DIRECT_P_MIXING. G-space "// &
              "mixers (BROYDEN/PULAY) only mix the grid density and leave rho_ao "// &
              "un-mixed, making the RIDFT SCF oscillate.")
   END IF
END IF
```

* Passes: `DIRECT_P_MIXING` (1), `no_mixing` (0/−1).
* Aborts: `gspace` (2), `pulay` (3), `broyden` (4), `multisecant` (5) — on SCF iteration 1.

### User guidance

For any RIDFT run use:

```
&SCF
  &MIXING
    METHOD DIRECT_P_MIXING
    ALPHA 0.4            ! lower (0.2-0.4) if the SCF is still stiff
  &END MIXING
&END SCF
```

---

## 8. Alternative not taken: internal matrix mixing (Option B)

A more general fix would let RIDFT work with *any* `&MIXING` method by mixing the density
matrix **inside** the RIDFT block:

```
P_ridft ← (1 − α)·P_ridft + α·rho_ao        ! α from scf_env%mixing_store%alpha
```

with `P_ridft` persisted across iterations (e.g. a new `p_mix_ridft` matrix set on
`scf_env`), and used for the RI-Hartree, the energy trace, and `build_core_erf_ridft`.

This was **not** implemented because:

* it duplicates mixing logic that `DIRECT_P_MIXING` already performs correctly;
* it must update **once per real SCF step only** — the RIDFT block currently also runs
  for `just_energy`/post-SCF evaluations, which would double-mix (see bug B5 in the
  developer bug list) and would have to be guarded first;
* RIDFT's internal density could then drift from the SCF's own convergence criterion.

It remains a viable future enhancement if Broyden/Pulay support for RIDFT is required.

---

## 9. Code reference index

| What | File:line |
|---|---|
| RIDFT block (sources `rho_ao_ridft` from `rho_ao(i,1)`) | `src/qs_ks_methods.F:961` |
| RIDFT passes `rho_ao` to `build_core_erf_ridft` | `src/qs_ks_methods.F:1376` |
| RIDFT zeroes the grid Hartree | `src/qs_ks_utils.F:1327` |
| Mixing guard (this fix) | `src/qs_ks_methods.F:~900` |
| `qs_scf_density_mixing` (matrix vs g-space branch) | `src/qs_scf_loop_utils.F:483-512` |
| `p_mix_new → rho_ao` copy (all mixers) | `src/qs_scf.F:587-595` |
| `qs_scf_rho_update` (grid recompute + grid mix) | `src/qs_scf_loop_utils.F:654-675` |
| Mixing-method constants | `src/qs_density_mixing_types.F:45-46` |
