# RIDFT — Developer Overview & Working Memory

> **Purpose.** This is my orientation/working-memory doc for developing the new
> **RIDFT** (Resolution-of-Identity DFT) method option in CP2K. It maps the code that
> matters — `qs_ks_methods.F`, `core_ae.F`, and every supporting module/subroutine —
> so that a future session can pick up work on the `ridft_claude` branch without
> re-deriving the layout. Everything here is anchored to **real line numbers in this
> checkout** (branch `ridft_claude`, forked from `claude/ridft-method-overview-3qcgr7`,
> which itself sits on the local `ri_dft` branch). Line numbers in the older uploaded
> notes are ~1–7 lines off; trust this file for navigation.
>
> Deep physics/derivations live in [`background/`](background/) (5 notes). This file is
> the *code map*; those are the *theory*. Cross-references below point into them.

---

## 0. TL;DR — what RIDFT is, in one screen

RIDFT computes the **Hartree channel without a plane-wave grid**. Everything else
(kinetic, the sharp `erfc` core Hamiltonian, XC, core-overlap) is **bit-identical to
GAPW / all-electron**. Only `energy%hartree` and how the KS Hartree potential is built
change.

GAPW builds Hartree from one grid Poisson solve of `rho_elec + rho_core`. RIDFT replaces
that with three analytic pieces:

| Physical term | GAPW (grid) | RIDFT (analytic, no grid) | Routine |
|---|---|---|---|
| electron–electron `1/r` | Poisson of `rho_elec` | **RI-V fit** `J = (μν│P) V⁻¹ (Q│λσ)ρ` | `get_hartree_noenv` (`rt_bse_utils.F`) |
| electron–nucleus, smooth `erf` | Poisson of `rho_core`×`rho_elec` | `verf` AO integrals → `e_nuc_erf` | **`build_core_erf_ridft`** (`core_ae.F`) |
| nucleus–nucleus, smooth `erf` | Poisson of `rho_core`×`rho_core` | `Σ_{A<B} Z_AZ_B erf(√α_AB R)/R` → `e_nn_erf` | **`build_core_erf_ridft`** |
| grid self-energy correction | `+ core_self` line | rebalanced: `- energy%core_self` folded in | `qs_ks_methods.F` |

Final assembly (`qs_ks_methods.F:1134,1210,1374`):

```
energy%hartree = E_ee^RI  +  e_nuc_erf  +  e_nn_erf  -  energy%core_self
```

The KS matrix gets: the RI `J` matrix (`v_ao_dbcsr`, added at `qs_ks_methods.F:1359`)
+ `V_nuc_erf` (added inside `build_core_erf_ridft`, `core_ae.F:775`). The grid Hartree
potential `v_rspace` is **zeroed** in `sum_up_and_integrate` (`qs_ks_utils.F:1327`) so it
is not double-counted. XC still flows through `v_rspace_new` untouched.

See [`background/nuclear_contributions_core_ae.md`](background/nuclear_contributions_core_ae.md)
§5 for the term-by-term H₂ numbers and [`background/xc_energy_gpw_gapw.md`](background/xc_energy_gpw_gapw.md)
for the XC path.

---

## 1. Status of the branch — this is unfinished WIP

Be aware before touching anything:

- **Debug noise everywhere.** The RIDFT block in `qs_ks_methods.F` and the guards in
  `qs_ks_utils.F` are littered with `PRINT *` and `CALL dbcsr_print(...)` dumps
  (dozens). These must be stripped before this is remotely mergeable. Grep
  `PRINT *,` inside the `ridft` regions.
- **Hardcoded RI metric** (`qs_ks_methods.F:947-956`): `cutoff_radius = 10.0`,
  `omega = 0.0`, `potential_type = do_potential_coulomb`, `filename = "t_c_g.dat"`,
  `unit_nr = 1`. None of this is exposed through the input yet.
- **Dead experimental block** `IF (.FALSE.) THEN … END IF` at
  `qs_ks_methods.F:985-1079` — an *old* MO-based density reconstruction. The `×2.0`
  spin scaling (`cp_cfm_scale(2.0, …)`, lines 1047-1048) lives inside this dead block,
  so it is **inactive**. The live density is taken directly from `rho_ao` (lines
  962-971, 980-981). Delete the dead block eventually.
- **Only `DIRECT_P_MIXING` works.** A guard aborts g-space mixers (see §5).
- **No forces / no virial** for the RIDFT Hartree (RI + nuclear erf) terms.
- **`build_core_erf_ridft` is serial** (`nthread=1`, no OpenMP locks).
- **QM/MM + RIDFT aborts** (`qs_ks_methods.F:657`).
- Message typo to fix: `"RIDFT does not suppport g-space mixing methods."`
  (`qs_ks_methods.F:900`).

The reason RIDFT "never seemed to do XC" is mundane: every test input uses
`&XC_FUNCTIONAL NONE`, so `qs_vxc_create` returns early. The XC path itself is correct
and shared — see §7.

---

## 2. Control-flow map — every place `qs_control%ridft` branches

The method is a single boolean `dft_control%qs_control%ridft`, set from
`method_id == do_method_ridft` (=21). Grep `ridft` (case-insensitive) to find all sites.
The behaviourally meaningful gates, in execution order inside
`qs_ks_build_kohn_sham_matrix` (`qs_ks_methods.F`):

| Line | Gate | Effect |
|---|---|---|
| 454 | `ridft = dft_control%qs_control%ridft` | local flag |
| 508 | `IF (.NOT. ridft)` | skip creating `v_hartree_gspace`/`rho_tot_gspace` + `DETAILED_ENERGY` grid solve |
| 534 | `IF (ridft)` | create `rho_tot_gspace`=0 and `v_hartree_rspace`=0 (so downstream code has zeroed objects) |
| 547 | *(all methods)* | `calc_rho_tot_gspace` still runs — builds `rho_g` used later for XC/density |
| 572 | `IF (.NOT. ridft)` | **skip the grid `pw_poisson_solve`** for `energy%hartree` |
| 594, 602 | `IF (.NOT. ridft)` | skip DDAPC / periodic-image decoupling on the grid Hartree |
| 612 | `IF (.NOT. ridft)` | skip `give_back_pw(v_hartree_gspace)` (never created for RIDFT) |
| 657 | `IF (qmmm .AND. ridft)` | **`CPABORT`** — QM/MM not supported |
| 886 | `IF (ridft)` | **the main RI-Hartree block** (§4) |
| 1317 | `IF (ridft)` | point `my_rho` at `rho_ao_ridft_copy_dbcsr` for `sum_up_and_integrate` |
| 1354 | `IF (ridft)` | add RI `J` (`v_ao_dbcsr`) to KS + call `build_core_erf_ridft` + assemble `energy%hartree` |
| 1510 | `IF (gapw .OR. gapw_xc .OR. ridft)` | shared post-step rho handling (comment flags it "maybe unnecessary for RIDFT") |

And outside this driver:

| File:line | Effect |
|---|---|
| `qs_ks_utils.F:1325-1328` | `sum_up_and_integrate` zeroes grid Hartree `v_rspace` for RIDFT |
| `qs_neighbor_lists.F:425` | `ridft` grouped with `lrigpw`/`rigpw` so `sab_orb`/`sac_ae`/`sac_ppl` get built |
| `qs_environment.F:1376` | dispersion env allocated for RIDFT too |
| `topology_coordinate_util.F:346` | "ignore outside box" warning — intended for RIDFT + partial charges outside the cell |

---

## 3. `core_ae.F` — the nuclear erf Hartree term

`core_ae.F` builds nuclear contributions from one shared integral engine,
[`verfc`](../../../src/aobasis/ai_verfc.F) (Obara–Saika three-center integrals). `verfc`
returns three blocks from the exact partition `1/r = erf(√α r)/r + erfc(√α r)/r`:
`vnuc` (bare `1/r`), `verf` (smooth Gaussian `erf/r`), `vabc` (`-Z erfc/r`). A driver
picks its physics by choosing which block it keeps. Full theory:
[`background/nuclear_contributions_core_ae.md`](background/nuclear_contributions_core_ae.md).

`PUBLIC` routines (`core_ae.F:68`): `build_core_ae`, **`build_core_erf_ridft`**,
`build_erfc`, `verfc_force`.

| Routine | Keeps | Writes | Used by |
|---|---|---|---|
| `build_core_ae` (`:90`) | `vabc` (`-Z erfc/r`, sharp) | core Hamiltonian `matrix_h` | `qs_core_matrices.F` — **shared, RIDFT uses it too** |
| `build_erfc` (`:928`) | `vabc`, params hand-fed | scratch `e_mat` | energy-decomposition analysis |
| **`build_core_erf_ridft`** (`:540`) | **`verf`** (`+erf/r`, smooth) | **`matrix_ks`** + returns `e_nuc_erf`,`e_nn_erf` | `qs_ks_methods.F:1369` |
| `verfc_force` (`:844`) | derivative integrals | forces | `build_core_ae` only |

### 3.1 `build_core_erf_ridft` — anatomy (`core_ae.F:540-809`)

Signature (`:540`):
```fortran
build_core_erf_ridft(matrix_ks, rho_ao, e_nuc_erf, e_nn_erf, &
                     qs_kind_set, atomic_kind_set, particle_set, &
                     sab_orb, sac_ae, sac_ppl, nimages, cell_to_index, basis_type)
```

Flow:
1. **Setup** (`:590-616`): build `basis_set_list`, allocate `hab`/`verf`/`vnuc`
   scratch, and create iterators over **two** neighbor lists — `ap_iterator` over
   `sac_ae` (all-electron/sgp cores) and `ap_iterator_ppl` over `sac_ppl` (GTH cores).
   Each is guarded by `ASSOCIATED` so a pure-GTH or pure-AE system works.
2. **Scratch KS clones** (`:618-622`): `v_nuc_erf(img)` = a zeroed copy of
   `matrix_ks(1,img)` per image. All erf integrals accumulate here (kept separate so a
   single `dbcsr_dot` later handles the symmetry factor + MPI sum cleanly).
3. **Orbital-pair loop** over `sab_orb` (`:624-763`): for each (μ on iatom, ν on jatom),
   loop kinds `kkind` and select the potential:
   - `all_potential` → `iter => ap_iterator` (`:685`)
   - `sgp_potential` → `iter => ap_iterator` (`:691`)
   - **`gth_potential` → `iter => ap_iterator_ppl`** (`:697`) — **Option A already
     implemented here.** GTH cores live in `sac_ppl`, not `sac_ae`.
   - else `CYCLE`.
   Pull `alpha_core_charge`, `zeff`, `ccore_charge`, `core_charge_radius` via
   `get_potential`, then call `verfc(... alpha_c, core_radius, 0.0, -core_charge, …)`
   (`:730`). The `-core_charge` prefactor selects/scales the **`verf`** (smooth) block;
   `vnuc` is filled but discarded. Distance screening at `:717,721,724`.
4. **Contract + scatter** (`:740-761`): primitive → spherical (`sphi_a`/`sphi_b`
   matmuls), added into the stored triangle of `v_nuc_erf` (`iatom<=jatom` branch).
5. **Energy trace + KS update** (`:770-779`): per image
   `e_nuc_erf += dbcsr_dot(rho_ao, v_nuc_erf)` (symmetry factor & MPI reduction happen
   **inside** `dbcsr_dot` — caller must **not** reduce again), then
   `dbcsr_add(matrix_ks, v_nuc_erf)` folds `V_nuc_erf` into the KS matrix.
6. **Nucleus–nucleus erf** (`:781-805`): cache `z_nn`,`alpha_nn` per kind via
   `get_qs_kind(zeff=…, alpha_core_charge=…)` (aggregating accessor — covers all/sgp/gth
   automatically). Unscreened `O(N²)` pair sum
   `e_nn_erf = Σ_{A<B} Z_AZ_B erf(√α_AB R)/R`, `α_AB = α_Aα_B/(α_A+α_B)`. `z_nn==0`
   kinds are skipped, so ghost/basis-only centers drop out.

**Key correctness facts**
- `e_nuc_erf` is negative (attractive), `e_nn_erf` positive.
- No `core_self` term is produced here — it is rebalanced by the caller (§5.2).
- No forces; serial.
- **GTH support is live** (uploaded `ridft_point_charges_fake_core.md` describes this as
  "to do", but it's *already done* in this checkout). What is **still open** is a *true*
  point charge (Option B: contract `verfc`'s discarded `vnuc` = `⟨μ│1/r│ν⟩` from an
  explicit charge list). See that note §4/§6 and
  [`background/ridft_point_charges_fake_core.md`](background/ridft_point_charges_fake_core.md).

---

## 4. `qs_ks_methods.F` — the driver and the RI-Hartree block

`qs_ks_build_kohn_sham_matrix` (`:228`) is *the* place RIDFT lives. Relevant declarations
are `:239-343` (note `ridft`, `e_nuc_erf`, `e_nn_erf`, the `*_ridft` neighbor-list
pointers, `rho_ao_ridft`/`rho_ao_ridft_copy` cfm arrays, `scf_env_ridft`, and the RI
scratch dbcsr `v_ao_dbcsr`, `v_dbcsr`, `int_3c_ptr`, `rho_ao_calc`).

### 4.1 The main RI-Hartree block (`:886-1247`)

1. **Env + mixing guard** (`:890-902`): fetch `scf_control`, `qs_kind_set`,
   `particle_set`, `mos`, `atomic_kind_set`, `scf_env_ridft`. Then:
   ```fortran
   IF (scf_env_ridft%mixing_method >= gspace_mixing_nr) &
      CPABORT("RIDFT does not suppport g-space mixing methods.")
   ```
   RIDFT feeds on the density **matrix** `rho_ao`; only `DIRECT_P_MIXING` mixes the
   matrix. See [`background/ridft_density_mixing.md`](background/ridft_density_mixing.md).
2. **RI index bookkeeping** (`:904-945`): `basis_set_list_setup` for `"ORB"` and
   `"RI_AUX"`; build `sizes_RI(iatom)`, `i_RI_start_from_atom`, `i_RI_end_from_atom`,
   total `n_RI`. (Copied from the GW `get_RI_basis_and_basis_function_indices` pattern.)
3. **RI metric** (`:947-956`): hardcoded (see §1).
4. **Density → cfm** (`:958-983`): `cp_libint_static_init`; make `rho_fm_struct`
   (`n_ao × n_ao`). **Live path:** `rho_ao_ridft(i)` and `rho_ao_ridft_copy(i)` are
   filled directly from `rho_ao(i,1)%matrix` (`copy_dbcsr_to_fm` → `cp_fm_to_cfm`), and
   `rho_ao_ridft_copy_dbcsr` (a real dbcsr) likewise. `rho_ao_ridft` = density for the RI
   `J` build; `rho_ao_ridft_copy_dbcsr` = density used both for the energy trace and, via
   `my_rho`, for grid-potential integration.
5. **Dead block** (`:985-1079`): `IF (.FALSE.)` — ignore/delete (see §1).
6. **3-center integrals** (`:1089`): `create_hartree_ri_3c(... int_3c=int_3c_ptr, n_ao,
   n_RI, basis_set_AO, basis_set_RI, i_RI_start_from_atom, ri_metric, qs_env, unit_nr)`
   → `(μν│P)`.
7. **2-center metric** (`:1119`): `create_2c_coulomb_rep_v(v_dbcsr, basis_set_RI, qs_env,
   n_RI, blacs_env, para_env, ri_metric)` → the inverted Coulomb metric `V_PQ`.
8. **RI Hartree per spin** (`:1130-1226`):
   - `energy%hartree = 0` (`:1134`) — **reset**, RIDFT owns this term.
   - `get_hartree_noenv(v_fm(i), rho_ao_ridft(i), int_3c_ptr, v_dbcsr, n_RI, sizes_RI,
     para_env, rho_ao_calc, v_ao_dbcsr)` (`:1141`) → `v_fm(i)` = RI `J` matrix in AO basis.
   - Copy `v_fm → v_ao_dbcsr`; transpose + `dbcsr_desymmetrize` → `v_hartree_T_unsym`.
   - `dbcsr_dot(rho_ao_ridft_copy_dbcsr, v_hartree_T_unsym, e_hartree_ri_local)` (`:1203`);
     `× 0.5`; `energy%hartree += e_hartree_ri_local` (`:1207-1210`). This is `E_ee^RI`.
9. **Cleanup** (`:1228-1245`): release transpose temporaries, deallocate RI arrays,
   `cp_libint_static_cleanup`.

### 4.2 KS assembly + final Hartree (`:1260-1402`)

- Core H init: `ks_matrix = matrix_h` (`:1289`) — shared with all methods.
- **`my_rho` swap** (`:1317-1326`): RIDFT → `rho_ao_ridft_copy_dbcsr`; else `rho_ao`.
  `my_rho` is what `sum_up_and_integrate` integrates the grid potentials (XC, …) against.
- **Add RI `J`** (`:1359`): `dbcsr_add(ks_matrix, v_ao_dbcsr)`.
- **Nuclear erf** (`:1362-1379`): fetch `sab_orb`/`sac_ae`/`sac_ppl` (+ `cell_to_index`
  for k-points), then
  ```fortran
  CALL build_core_erf_ridft(ks_matrix, rho_ao, e_nuc_erf, e_nn_erf, …, "ORB")
  energy%hartree = energy%hartree + e_nuc_erf + e_nn_erf - energy%core_self
  ```
- **`sum_up_and_integrate`** (`:1387`): integrates XC etc. into KS; internally zeroes the
  grid Hartree for RIDFT (`qs_ks_utils.F:1327`).

---

## 5. Density mixing — why only `DIRECT_P_MIXING`

Full analysis: [`background/ridft_density_mixing.md`](background/ridft_density_mixing.md).

- CP2K mixes **either** the density matrix `rho_ao` (`DIRECT_P_MIXING`) **or** the grid
  density `rho_g` (Broyden/Pulay/…), never both.
- RIDFT's Hartree feedback depends **only** on `rho_ao` and it **zeroes the grid
  Hartree**. Under a g-space mixer, `rho_ao` is the raw un-mixed `P_out` → undamped
  feedback → a **period-2 limit cycle** (charge sloshing) for stiff systems.
- Guard: `qs_ks_methods.F:898-900` aborts if `mixing_method >= gspace_mixing_nr`
  (constants in `qs_density_mixing_types.F`: direct=1, gspace=2, pulay=3, broyden=4,
  multisecant=5). Passes: `DIRECT_P_MIXING`, `no_mixing`.
- Recommended input: `&MIXING METHOD DIRECT_P_MIXING / ALPHA 0.4 &END`.
- "Option B" (mix a persisted RIDFT matrix internally so any mixer works) is documented
  but **not** implemented.

### 5.2 The `core_self` rebalance

Total energy is a plain sum `core + core_overlap + core_self + exc + hartree`, with
`core_self` (negative) already a line item. GAPW *needs* it (its grid Poisson
over-counts Gaussian self-interaction). RIDFT's `e_nn_erf` is a clean `Σ_{A<B}` pair sum
with **no** self-interaction, so the `- energy%core_self` folded into `energy%hartree`
(`:1374`) cancels the `+ core_self` the total will add back. Net: RIDFT carries the true
nuclear repulsion as the explicit analytic `e_nn_erf`.

---

## 6. Supporting modules & subroutines

### 6.1 `rt_bse_utils.F` — the RI-V machinery (borrowed from RT-BSE / GW, author S. Marek)

RIDFT reuses three routines, made env-free so they don't need `bs_env`:

| Routine (line) | Role |
|---|---|
| `get_hartree_noenv` (`:109`) | Builds RI `J` in AO basis: `Q_Q = (Q│λσ)ρ`, then `J_μν = (μν│P) V⁻¹ Q`. Env-free twin of `get_hartree`. |
| `create_hartree_ri_3c` (`:265`) | Allocates + fills the 3-center integrals `(μν│P)` as a plain array (no DBT). |
| `create_2c_coulomb_rep_v` (`:347`) | Builds + inverts the 2-center Coulomb metric `V_PQ` (`= init_hartree` without `bs_env`). |

RI-V is variational (Coulomb metric) → a strict **lower bound**, which is why RIDFT's
`E_ee^RI` sits *below* GAPW's grid value (the residual ~mHa is intrinsic RI fit error,
not a bug). Theory: `nuclear_contributions_core_ae.md` §5.6.

### 6.2 XC — unchanged, shared with GPW/GAPW

RIDFT calls the **same** `qs_vxc_create` (`qs_ks_methods.F:792`) on the same auxbas/xc
grid; `energy%exc` is filled identically to GPW. The density grids `rho_r`/`rho_g` that XC
needs **are** built for RIDFT (`qs_rho_update_rho_low` takes the final `ELSE` branch, full
valence density collocated). The GAPW atomic XC (`calculate_vxc_atom`, `energy%exc1`) is
gated by `IF (gapw .OR. gapw_xc)` so it is correctly **skipped**. **No RIDFT-specific XC
code exists or is wanted** — the faithful implementation is the empty diff. Details:
[`background/ridft_exchange_correlation.md`](background/ridft_exchange_correlation.md) and
[`background/xc_energy_gpw_gapw.md`](background/xc_energy_gpw_gapw.md).

### 6.3 `qs_ks_utils.F` — grid Hartree zeroing

`sum_up_and_integrate` (`:1318+`): for RIDFT, `CALL pw_zero(v_rspace)` (`:1327`) — this is
`v_hartree_rspace`, **not** the XC potential `v_rspace_new`. Prevents double-counting the
Hartree the RI path already added to KS. (Also carries a debug `PRINT`.)

### 6.4 Input & control plumbing

| File:line | Role |
|---|---|
| `input_constants.F:271` | `do_method_ridft = 21` |
| `input_cp2k_qs.F:441,459` | registers `"RIDFT"` keyword → `do_method_ridft` in `&QS METHOD` enum |
| `cp_control_types.F:397` | `LOGICAL :: ridft = .FALSE.` in `qs_control_type` |
| `cp_control_utils.F:1013,1031-1032` | init `ridft=.FALSE.`; `CASE (do_method_ridft) → ridft=.TRUE.` |
| `qs_environment.F:1376` | include `do_method_ridft` in dispersion-env setup |
| `qs_neighbor_lists.F:51,425` | `ridft` grouped with lri/rigpw so orbital/`sac_ae`/`sac_ppl` lists exist |
| `topology_coordinate_util.F:346` | ignore-outside-cell warning, for RIDFT + external partial charges |
| `input_cp2k_print_dft.F:64` | imports `do_method_ridft` for method/printing |

---

## 7. Energy bookkeeping cheat-sheet (H₂ reference)

From `nuclear_contributions_core_ae.md` §5.3 (H₂, `XC NONE`), the RIDFT Hartree
decomposition:

| Sub-term | value (Ha) | source |
|---|---|---|
| `E_ee^RI` | +1.5731 | `get_hartree_noenv` |
| `e_nuc_erf` | −3.8306 | `build_core_erf_ridft`, `verf` |
| `e_nn_erf` | +0.7151 | `build_core_erf_ridft`, pair loop |
| `− core_self` | +2.8209 | analytic Gaussian self-energy |
| **RIDFT `energy%hartree`** | **1.2785** | sum |
| GAPW `energy%hartree` | 1.2932 | grid Poisson |

Residual **14.7 mHa** is RI fit error in `E_ee^RI` only (nuclear pieces are analytic/exact).
**When comparing to GAPW, compare RIDFT `energy%hartree` against GAPW
`energy%hartree + energy%hartree_1c`** (the 1-center term), not `energy%hartree` alone.

---

## 8. Verification anchors

1. **XC-off regression:** existing RIDFT test (`&XC_FUNCTIONAL NONE`) — source unchanged
   ⇒ bit-identical.
2. **XC-on:** set a real functional in matching RIDFT and GPW inputs (same cell/cutoff/
   basis); `energy%exc` must be equal (same routine, same grid).
3. **`POTENTIAL ALL` H₂** (GTH branch inactive): `e_nuc_erf ≈ −3.83`, `e_nn_erf ≈ 0.715`,
   `energy%hartree ≈ 1.2785`.
4. **Real GTH atom:** RIDFT `energy%hartree` ≈ GAPW `energy%hartree + hartree_1c` within
   the RI-V residual.
5. **Convergence:** `DIRECT_P_MIXING` converges (~11 steps for the H₂/aug-cc-pVTZ case);
   Broyden must abort at the guard.

---

## 9. Cleanup / TODO backlog (before this is mergeable)

- [ ] Strip all `PRINT *` / `dbcsr_print` debug from the `ridft` regions of
      `qs_ks_methods.F` and `qs_ks_utils.F`.
- [ ] Delete the dead `IF (.FALSE.)` block (`qs_ks_methods.F:985-1079`).
- [ ] Expose the RI metric via input instead of hardcoding
      (`cutoff_radius`, potential type, `t_c_g.dat`, `unit_nr`).
- [ ] Fix the abort-message typo `"suppport"` (`:900`).
- [ ] Implement **forces / virial** for the RI-Hartree and nuclear-erf terms
      (currently missing → geometry opt / MD impossible under RIDFT).
- [ ] Thread `build_core_erf_ridft` (OpenMP, like `build_core_ae`).
- [ ] Validate/clean **spin** handling (`spin_degeneracy`, open-shell) — currently
      closed-shell-centric; the `×2.0` scaling logic is buried in the dead block.
- [ ] Test/verify the **k-points** path (`dokp`, `cell_to_index`) — plumbed but untested.
- [ ] Decide on **true point charges** (Option B: `vnuc` `1/r` from an explicit list,
      likely homed in `external_potential_types.F`); GTH-core (Option A) is already done.
- [ ] Consider Option B mixing (internal matrix mix) if Broyden/Pulay support is needed.
- [ ] Revisit the `IF (gapw .OR. gapw_xc .OR. ridft)` post-step at `:1510`
      ("maybe unnecessary for RIDFT").

---

## 10. File index (quick jump)

| File | What to look at |
|---|---|
| `src/qs_ks_methods.F` | driver; RI-Hartree block `:886-1247`; KS assembly `:1260-1402`; gates in §2 |
| `src/core_ae.F` | `build_core_erf_ridft` `:540-809`; `build_core_ae` `:90`; `verfc_force` `:844` |
| `src/aobasis/ai_verfc.F` | `verfc` integral engine (`vnuc`/`verf`/`vabc` blocks) |
| `src/rt_bse_utils.F` | `get_hartree_noenv` `:109`; `create_hartree_ri_3c` `:265`; `create_2c_coulomb_rep_v` `:347` |
| `src/qs_ks_utils.F` | `sum_up_and_integrate` grid-Hartree zero `:1325-1328` |
| `src/qs_vxc.F`, `src/xc/*` | shared XC path (unchanged) |
| `src/qs_core_energies.F` | `calculate_ecore_self` (`core_self`) |
| input/control | `input_constants.F`, `input_cp2k_qs.F`, `cp_control_types.F`, `cp_control_utils.F`, `qs_environment.F`, `qs_neighbor_lists.F` |
| `.claude/memories/ridft/background/` | 5 deep-dive notes (physics + derivations) |

---

*Maintainer note: keep this file's line numbers in sync when editing the RIDFT regions.
If numbers drift, re-anchor with `grep -n ridft src/qs_ks_methods.F` and
`grep -n build_core_erf_ridft src/core_ae.F`.*
