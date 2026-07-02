# Point charges in RIDFT — the "fake H core" trick and why it misbehaves

**Scope:** what happens to the Hartree energy in `METHOD RIDFT` when the input adds a
charge-only "ghost" atom such as

```text
&KIND 0partial_charge          ! fake H atom, no wavefunctions in the QM calc
  ELEMENT H
  BASIS_SET     NONE
  BASIS_SET RI_AUX NONE
  CORE_CORRECTION <q>          ! partial charge, offset against the H core (Z-1)
  POTENTIAL GTH-PBE-q1
&END KIND
```

The short answer: **in the current RIDFT implementation this fake core is electrostatically
invisible to the QM electrons.** Its charge contributes *nothing* to `energy%hartree`, nothing
to the KS matrix, and nothing to the nucleus–nucleus repulsion. Below is the mechanism, why,
and two ways to do it correctly.

---

## 1. How the nuclear charge is *supposed* to reach the Hartree energy in RIDFT

RIDFT deliberately removes the plane-wave grid, so it cannot use GAPW/GPW's route
(`calculate_rho_core` → `rho_core` Gaussian on the grid → Poisson solve). Instead the nuclear
electrostatics is rebuilt **analytically in the AO basis** by `build_core_erf_ridft`
(`src/core_ae.F:539`), using the exact Coulomb partition

```
1/r = erf(√α·r)/r  +  erfc(√α·r)/r
```

- the **erfc** (sharp) half is the point-nucleus attraction → already in the core Hamiltonian
  (`build_core_ae` / `build_core_ppl`);
- the **erf** (smooth) half is the potential of a *Gaussian-smeared* nucleus of width
  `α = 1/(2 r_c²)` → this is exactly the `−rho_core` term GAPW Poisson-solves, and it is the
  piece RIDFT must add back.

`build_core_erf_ridft` therefore produces two contributions that flow into the Hartree energy:

1. **electron–nucleus** `−Z·⟨μ|erf(√α|r−R_A|)/|r−R_A| |ν⟩` → added to the KS matrix; its trace
   gives `e_nuc_erf` (negative, attractive);
2. **nucleus–nucleus** `Σ_{A<B} Z_A Z_B·erf(√α_{AB}·R)/R` → scalar `e_nn_erf` (positive);

and the caller assembles
```
energy%hartree += e_nuc_erf + e_nn_erf − energy%core_self
```
(`energy%core_self` is the Gaussian self-energy `−Σ_A ½ Z_A² √(2α_A/π)`).

So a charge can only be felt by the QM electrons in RIDFT if it appears in **`build_core_erf_ridft`**.

---

## 2. What the "fake H core" input actually creates

| input line | effect |
|---|---|
| `BASIS_SET NONE` | no orbital basis on this kind → it never appears as μ or ν in any AO integral, contributes no electrons, no RI density |
| `BASIS_SET RI_AUX NONE` | no auxiliary functions → does not participate in the RI fit of the e–e density |
| `CORE_CORRECTION <q>` | parsed as `zeff_correction` (`qs_kind_types.F:1854`); it **offsets** the potential's `zeff` and the electron-count bookkeeping (`nelectron = … − zeff_correction`, line 1190). It does *not* by itself create a standalone classical charge — it only shifts the charge carried by whatever `POTENTIAL` is attached. |
| `POTENTIAL GTH-PBE-q1` | allocates a **`gth_potential`** (not `all_potential`/`sgp_potential`). The net core charge `q = zeff − zeff_correction` is carried by the GTH local pseudopotential. |

So the intended physics is: a positive (or fractional) point-like charge sitting in the QM
region, with no orbitals of its own, whose field the real electrons feel through Hartree.

---

## 3. Why it fails in RIDFT — the GTH path is skipped

`build_core_erf_ridft` reads the core-charge parameters **only** from `all_potential` or
`sgp_potential`:

```fortran
! src/core_ae.F:687-700
DO kkind = 1, nkind
   CALL get_qs_kind(qs_kind_set(kkind), all_potential=all_potential, &
                    sgp_potential=sgp_potential)
   IF (ASSOCIATED(all_potential)) THEN
      CALL get_potential(all_potential, alpha_core_charge=alpha_c, zeff=zeta_c, ...)
   ELSE IF (ASSOCIATED(sgp_potential)) THEN
      CALL get_potential(sgp_potential, alpha_core_charge=alpha_c, zeff=zeta_c, ...)
   ELSE
      CYCLE          ! <-- a gth_potential kind lands here and is dropped
   END IF
   ...
```

The same restriction applies to the nucleus–nucleus loop: it pre-caches `z_nn(ikind)` only from
`all_potential`/`sgp_potential` (`core_ae.F:797-805`) and then does `IF (z_nn(ikind)==0.0) CYCLE`
(lines 810, 813). A GTH kind has `z_nn = 0` here, so it is excluded from `e_nn_erf` too.

Consequences for the fake `GTH-PBE-q1` core:

1. **No electron–nucleus term.** Its `q` never enters the `verfc` accumulation, so
   `e_nuc_erf` and the KS matrix contain *zero* attraction toward the fake charge. The QM
   electrons simply do not see it.
2. **No nucleus–nucleus term.** It is skipped in the `e_nn_erf` loop, so the repulsion between
   the fake charge and the real nuclei is missing.
3. **No grid fallback.** RIDFT zeroes `v_rspace` in `sum_up_and_integrate`
   (`qs_ks_utils.F`), so the grid `rho_core` path that *would* have carried a GTH core in
   GAPW/GPW is also gone.

Net: **`energy%hartree` is wrong by exactly the (e–n attraction + n–n repulsion) of the fake
charge, and the SCF feels no embedding field at all.** This is an unwanted effect, and it is
silent — there is no warning; the run completes and reports a Hartree energy that omits the
charge.

> **Important corollary (bigger than the fake atom):** the current RIDFT does **not** yet
> support *any* `gth_potential` atom's Hartree contribution. All the RIDFT validation so far
> used `POTENTIAL ALL` (all-electron H₂). The moment a real GTH-pseudopotential atom is used,
> its `rho_core`/erf-Gaussian Hartree term is missing for the same reason. The fake-core idea
> just exposes this gap in its purest form.

### Even in GAPW the trick is imperfect
For completeness: the fake-GTH-core trick is not a clean point charge even in GAPW/GPW. A GTH
local PP is `V_loc(r) = −Zeff·erf(r/√2 r_loc)/r + (C₁ + C₂(r/r_loc)² + …)·exp(−(r/r_loc)²/2)`.
- the `erf` part is a **Gaussian-smeared** charge of width `r_loc`, not a true point charge;
- the polynomial `C_i` part still produces `⟨μ|V_ppl|ν⟩` integrals on the *real* atoms' basis
  (`calculate_ppl_rspace`), so short-range pseudopotential structure leaks into the result.

So "fake H core = point charge" is an approximation regardless of method; RIDFT just makes it
*zero* instead of *approximate*.

---

## 4. How to simulate point charges correctly

There are two clean routes. Both keep RIDFT grid-free and reuse machinery that already exists.

### Option A — extend `build_core_erf_ridft` to also read GTH core charges *(small change, smeared charge)*

Add a `gth_potential` branch to the `kkind` loop and to the `z_nn` pre-cache, pulling
`zeff = zeff − zeff_correction` and `α = 1/(2 r_loc²)` from the GTH local part. This makes both
**real** GTH atoms and the **fake** GTH core work, and it is the minimal fix consistent with the
existing erf design.

Caveats: the resulting charge is Gaussian-smeared with width `r_loc` (≈0.4–0.6 bohr for typical
PBE PPs), and the GTH short-range `C_i` terms still enter through `build_core_ppl`. For a *clean*
classical charge you would supply a custom GTH entry with all `C_i = 0` and a chosen `r_loc`, and
accept the Gaussian width.

### Option B — a dedicated point-charge channel using `verfc`'s `1/r` output *(recommended, true point charge)*

`verfc` already computes the bare Coulomb integral `⟨μ|1/r|ν⟩` and returns it as its `vnuc`
output — which `build_core_erf_ridft` currently **discards** (`core_ae.F:736`, the `vnuc`
argument). A genuine point charge `q` at `R` contributes:

- electron side: `−q·⟨μ|1/|r−R| |ν⟩` (the `vnuc` integral, no erf, no Gaussian width);
- nucleus side: `Σ_A q·Z_A / |R−R_A|` (bare Coulomb) and `Σ q·q' / |R−R'|` for charge pairs;
- **no `core_self` term** — a point charge has no Gaussian self-energy, so the
  `− energy%core_self` correction is simply absent for these centers.

This is the cleanest fit for RIDFT: it is exact (a true `1/r` point charge), grid-free, adds no
self-energy bookkeeping, and reuses the `verfc` call that is already in the loop — you keep the
`vnuc` branch instead of throwing it away.

The charges should come from an **explicit input list** (charge + position), not from abusing a
`&KIND` with `BASIS_SET NONE`. A `&KIND`-based ghost has to be threaded through the whole
qs_kind / potential / neighbor-list machinery (`sac_ae` vs `sac_ppl`, `zeff`, electron counting)
where it keeps getting filtered out; an explicit charge list sidesteps all of it and is what a
"point charge" actually is.

> Note: `src/subsys/external_potential_types.F` is already on the `ri_dft` branch's modified
> list — an explicit external point-charge type is the natural home for Option B's input and
> data, feeding a new `build_point_charge_ridft` routine modelled on `build_core_erf_ridft` but
> contracting the `vnuc` (1/r) integral instead of `verf`.

### Recommendation

| goal | use |
|---|---|
| Make real GTH-pseudopotential atoms work in RIDFT at all | **Option A** (needed regardless) |
| Embed genuine classical point charges in the QM region | **Option B** (explicit list + `vnuc`/1·r), do **not** rely on a fake `&KIND` |

Do **not** ship the fake `GTH-PBE-q1` ghost as the point-charge mechanism for RIDFT: today it
contributes nothing, and even once Option A lands it is a smeared, PP-contaminated charge that
is hard to reason about.

---

## 5. Quick verification checklist (when implementing)

1. With a fake charge present but `XC NONE` and the real atoms held fixed, RIDFT `energy%hartree`
   must shift by exactly `e_q-nuc(en) + e_q-nuc(nn)` relative to the no-charge run — and that
   shift must match the analytic `−q·Tr[P·⟨μ|1/r|ν⟩] + Σ q·Z_A/R_A` (Option B) or the GAPW
   reference using the same smeared GTH core (Option A).
2. Confirm the charge appears in the KS matrix diagonal near the real atoms (non-zero
   `v_nuc`/`v_point` block), i.e. the electrons actually feel it.
3. Electron count: with `CORE_CORRECTION` the SCF electron number is `Σ(zeff − zeff_correction)`;
   make sure the fake center adds **no** electrons (it has no basis) and that the cell does not
   silently charge up (CP2K warns about this for `CORE_CORRECTION ≠ 0`).

---

## 6. Concrete implementation — Option A (GTH support in `build_core_erf_ridft`)

This is the change that makes both **real GTH-pseudopotential atoms** and the **fake GTH core**
work. It is small and reuses the existing `verfc` body unchanged.

### 6.0 Why it does not double-count (verified)

`build_core_ppl` (`src/core_ppl.F`) computes **only** the C-polynomial short-range local part
(`cexp_ppl`/`alpha_ppl`) + extended-local + nonlocal projectors. It contains no `zeff`, no
`erf`, no `alpha_core_charge`. The GTH local `−Zeff·erf(r/√2 r_loc)/r` long-range term is carried
by `rho_core` (`calculate_rho_core`), exactly as in the all-electron case. So the erf core lives
in the Hartree/grid channel, never in `H_core`. Adding it analytically in `build_core_erf_ridft`
correctly **replaces** the grid `rho_core` that RIDFT removed — it does not double-count
`build_core_ppl`.

```
all-electron:  build_core_ae  (erfc → H_core)  + rho_core (erf → Hartree)
GTH:           build_core_ppl (C-poly → H_core) + rho_core (erf → Hartree)
RIDFT swaps rho_core for the analytic erf integral in BOTH cases.
```

> Bonus fix: `calculate_ecore_self` (`qs_core_energies.F:401`) already includes GTH atoms in
> `energy%core_self` (it aggregates `zeff`/`alpha_core_charge` over any potential). So today, with
> a GTH atom, `energy%hartree += e_nuc_erf + e_nn_erf − energy%core_self` **over-subtracts**
> (core_self has the GTH term, the erf terms don't). Closing the gap balances this automatically.

### 6.1 Key fact — GTH cores live in `sac_ppl`, not `sac_ae`

`sac_ae` is built from orbital × (all/sgp)-present centers; GTH centers are in `sac_ppl`
(`qs_neighbor_lists.F:690-735`). `gth_potential` exposes the *same* accessors as `all_potential`
(`alpha_core_charge`, `ccore_charge`, `core_charge_radius`, `zeff`), with `α = 1/(2 r_loc²)`.
So the only structural change is: iterate a **second** neighbor list for GTH centers.

### 6.2 `src/core_ae.F` changes

**(a) Import** — add `gth_potential_type` to the `external_potential_types` USE list:
```fortran
USE external_potential_types,  ONLY: all_potential_type, get_potential, &
                                     gth_potential_type, sgp_potential_type
```

**(b) Signature** — add `sac_ppl`:
```fortran
SUBROUTINE build_core_erf_ridft(matrix_ks, rho_ao, e_nuc_erf, e_nn_erf, &
                                 qs_kind_set, atomic_kind_set, particle_set, &
                                 sab_orb, sac_ae, sac_ppl, nimages, cell_to_index, basis_type)
   ...
   TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb, sac_ae, sac_ppl
```

**(c) Declarations** — a second iterator, a GTH pointer, a selector pointer:
```fortran
TYPE(neighbor_list_iterator_p_type), DIMENSION(:), POINTER :: ap_iterator, ap_iterator_ppl, iter
TYPE(gth_potential_type), POINTER                          :: gth_potential
```

**(d) Iterator creation** — guard each list (a pure-GTH system has NULL `sac_ae`):
```fortran
NULLIFY (ap_iterator, ap_iterator_ppl)
IF (ASSOCIATED(sac_ae))  CALL neighbor_list_iterator_create(ap_iterator,     sac_ae,  search=.TRUE., nthread=1)
IF (ASSOCIATED(sac_ppl)) CALL neighbor_list_iterator_create(ap_iterator_ppl, sac_ppl, search=.TRUE., nthread=1)
```

**(e) kkind nuclear loop** — replace the `all/sgp` block (`core_ae.F:687-700`) so GTH selects the
`sac_ppl` iterator; everything inside the `DO WHILE` (the `verfc` call) is unchanged, just driven
by `iter`:
```fortran
DO kkind = 1, nkind
   NULLIFY (all_potential, sgp_potential, gth_potential)
   CALL get_qs_kind(qs_kind_set(kkind), all_potential=all_potential, &
                    sgp_potential=sgp_potential, gth_potential=gth_potential)
   IF (ASSOCIATED(all_potential)) THEN
      CALL get_potential(all_potential, alpha_core_charge=alpha_c, zeff=zeta_c, &
                         ccore_charge=core_charge, core_charge_radius=core_radius)
      iter => ap_iterator            ! all-electron core: sac_ae
   ELSE IF (ASSOCIATED(sgp_potential)) THEN
      CALL get_potential(sgp_potential, alpha_core_charge=alpha_c, zeff=zeta_c, &
                         ccore_charge=core_charge, core_charge_radius=core_radius)
      iter => ap_iterator            ! sgp erf core also lives in sac_ae
   ELSE IF (ASSOCIATED(gth_potential)) THEN
      IF (.NOT. ASSOCIATED(ap_iterator_ppl)) CYCLE
      CALL get_potential(gth_potential, alpha_core_charge=alpha_c, zeff=zeta_c, &
                         ccore_charge=core_charge, core_charge_radius=core_radius)
      iter => ap_iterator_ppl        ! GTH core: sac_ppl
   ELSE
      CYCLE
   END IF

   CALL nl_set_sub_iterator(iter, ikind, kkind, iatom)
   DO WHILE (nl_sub_iterate(iter) == 0)
      CALL get_iterator_info(iter, jatom=katom, r=rac)
      ...                            ! unchanged verfc body, screened by core_radius
   END DO
END DO
```

**(f) Nucleus–nucleus pre-cache** (`core_ae.F:796-805`) — collapse to the aggregating accessor
that already covers all/sgp/**gth** (same call `calculate_ecore_self` uses):
```fortran
DO ikind = 1, nkind
   CALL get_qs_kind(qs_kind_set(ikind), zeff=z_nn(ikind), alpha_core_charge=alpha_nn(ikind))
END DO
```
The rest of the `e_nn_erf` loop (the `z_nn==0 ⇒ CYCLE` guard for chargeless kinds and the
`erf(√α_AB R)/R` sum) is unchanged and now picks up GTH atoms automatically.

**(g) Cleanup** — release both iterators:
```fortran
IF (ASSOCIATED(ap_iterator))     CALL neighbor_list_iterator_release(ap_iterator)
IF (ASSOCIATED(ap_iterator_ppl)) CALL neighbor_list_iterator_release(ap_iterator_ppl)
```

### 6.3 `src/qs_ks_methods.F` caller changes (≈ line 1386)

```fortran
NULLIFY (sab_orb_ridft, sac_ae_ridft, sac_ppl_ridft, cell_to_index_ridft, kpoints_ridft)
CALL get_qs_env(qs_env, sab_orb=sab_orb_ridft, sac_ae=sac_ae_ridft, &
                sac_ppl=sac_ppl_ridft, kpoints=kpoints_ridft)
...
IF (ASSOCIATED(sac_ae_ridft) .OR. ASSOCIATED(sac_ppl_ridft)) THEN
   CALL build_core_erf_ridft(ks_matrix, rho_ao, e_nuc_erf, e_nn_erf, &
                             qs_kind_set, atomic_kind_set, particle_set, &
                             sab_orb_ridft, sac_ae_ridft, sac_ppl_ridft, &
                             nimages, cell_to_index_ridft, "ORB")
   energy%hartree = energy%hartree + e_nuc_erf + e_nn_erf - energy%core_self
END IF
```
(declare `sac_ppl_ridft` as `TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER`.)

### 6.4 Limitations of Option A
- The GTH core is a **Gaussian** of width `r_loc` (~0.4–0.6 bohr), *not* a true point charge.
- `sac_ppl`'s build range is `ppl_radius` (≥ the erf core's `core_charge_radius`), so the
  potential's own `core_radius` screen governs the cutoff and nothing is wrongly culled.
- For a **clean** classical point charge, prefer Option B (§4): contract `verfc`'s discarded
  `vnuc` (= `⟨μ|1/r|ν⟩`) for a bare `1/r` center driven from an explicit charge list, with no
  Gaussian width and no `core_self` term.

### 6.5 Regression / verification anchors
1. **POTENTIAL ALL H₂** (GTH branch inactive): `e_nuc_erf ≈ −3.83`, `e_nn_erf ≈ 0.715`,
   `energy%hartree ≈ 1.2785` must be **bit-unchanged** — the guard that the refactor is inert
   when no GTH atom is present.
2. **Real GTH H₂**: RIDFT `energy%hartree` should match GAPW `energy%hartree + hartree_1c` within
   the known RI-V fit residual (~mHa).
3. **Fake-core run**: `energy%hartree` shifts by exactly the fake charge's (e–n + n–n)
   electrostatics, and the KS diagonal near the real atoms gains the corresponding attraction.
```
