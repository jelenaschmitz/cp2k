# Nuclear Contributions to the Core Hamiltonian in CP2K

### A guided tour of `core_ae.F`, the `verfc` integral engine, and how RIDFT fits in

This document explains, from the physics down to the Fortran, how CP2K includes the
**electron–nucleus** and **nucleus–nucleus** Coulomb interactions in the Kohn–Sham (KS)
matrix. It focuses on the three `build_…` driver routines in
[`src/core_ae.F`](../src/core_ae.F) and the shared integral kernel
[`verfc`](../src/aobasis/ai_verfc.F#L138) that all of them call.

It is written to be read top-to-bottom by someone new to the code and to be
presentable to others. Every non-trivial claim points either to a source line or a
paper.

---

## 1. The physical problem

In an electronic-structure calculation the nucleus enters the Hamiltonian through the
bare attractive Coulomb operator

$$
\hat V_{\text{en}}(\mathbf r) \;=\; -\sum_A \frac{Z_A}{|\mathbf r - \mathbf R_A|},
$$

where $Z_A$ is the charge of nucleus $A$ at position $\mathbf R_A$. In a Gaussian
basis we need the matrix elements

$$
V_{\mu\nu} \;=\; \Big\langle \mu \,\Big|\, -\sum_A \frac{Z_A}{|\mathbf r-\mathbf R_A|} \,\Big|\, \nu \Big\rangle .
$$

The difficulty is the $1/r$ singularity at the nucleus. CP2K's mixed
**Gaussian-and-plane-waves (GPW/GAPW)** scheme cannot put a sharp $1/r$ cusp on a
finite plane-wave grid, so it **splits the Coulomb kernel into a smooth piece and a
sharp piece** and treats each with the tool best suited to it.

### 1.1 The exact Coulomb partition (the central idea)

For any positive screening parameter $\alpha$ (units of length$^{-2}$):

$$
\boxed{\;\frac{1}{r} \;=\; \underbrace{\frac{\operatorname{erf}(\sqrt{\alpha}\,r)}{r}}_{\text{smooth, long range}} \;+\; \underbrace{\frac{\operatorname{erfc}(\sqrt{\alpha}\,r)}{r}}_{\text{sharp, short range}}\;}
$$

This identity is exact (`erf + erfc = 1`). Its two halves have complementary behaviour:

| piece | at $r\to 0$ | at $r\to\infty$ | physical picture |
|---|---|---|---|
| $\operatorname{erf}(\sqrt\alpha r)/r$ | finite (smooth) | $\to 1/r$ | potential of a **Gaussian charge cloud** |
| $\operatorname{erfc}(\sqrt\alpha r)/r$ | $\to 1/r$ (singular) | decays $\propto e^{-\alpha r^2}$ | short-range **point-minus-Gaussian** correction |

The key fact that makes this useful: the electrostatic potential of a **normalized
Gaussian charge** of total charge $Z$ and exponent $\alpha$ is *exactly* the erf half,

$$
\rho_{\text{core},A}(\mathbf r) = Z_A\Big(\tfrac{\alpha_A}{\pi}\Big)^{3/2} e^{-\alpha_A|\mathbf r-\mathbf R_A|^2}
\;\;\Longrightarrow\;\;
V_{\text{Gauss},A}(r) = \frac{Z_A\,\operatorname{erf}(\sqrt{\alpha_A}\,r)}{r}.
$$

So CP2K models each nucleus as a Gaussian charge cloud `rho_core` whose potential is
the **smooth erf part**, and adds back the **sharp erfc part** analytically. This is
the heart of GAPW [Lippert 1999; VandeVondele 2005].

### 1.2 Who handles which half

| Coulomb half | how it is evaluated | routine |
|---|---|---|
| **erf** (smooth, Gaussian) | GAPW: put `rho_core` on the grid, solve Poisson | grid / `qs_ks_utils.F` |
| **erf** (smooth, Gaussian) | RIDFT: analytic AO integrals (no grid) | **`build_core_erf_ridft`** |
| **erfc** (sharp, short range) | analytic AO integrals | **`build_core_ae`** |

This is the whole story in one line: **the three `build_…` routines and the grid each
take responsibility for a different slice of the same $-Z/r$ operator, chosen so that
no slice is hard to integrate.**

### 1.3 The Gaussian-nucleus model and its parameters

The whole construction rests on replacing the point nucleus by a **normalized Gaussian
charge cloud**. It is worth seeing exactly what that cloud is and where its parameters
come from in the code, because those same numbers (`alpha_c`, `core_charge`,
`core_radius`, `zeta_c`) are the arguments handed to `verfc`.

A point charge $Z_A\,\delta(\mathbf r-\mathbf R_A)$ is smeared into

$$
\rho_{\text{core},A}(\mathbf r) \;=\; N_A\,e^{-\alpha_A|\mathbf r-\mathbf R_A|^2},
\qquad N_A = Z_A\Big(\frac{\alpha_A}{\pi}\Big)^{3/2},
$$

where $N_A$ is fixed by charge conservation $\int\rho_{\text{core},A}\,d^3r = Z_A$. Its
electrostatic potential follows from Poisson's equation for a spherical Gaussian and is
the **erf** half of the Coulomb kernel exactly:

$$
V_{\text{core},A}(r) \;=\; \int \frac{\rho_{\text{core},A}(\mathbf r')}{|\mathbf r-\mathbf r'|}\,d^3r'
\;=\; \frac{Z_A\,\operatorname{erf}(\sqrt{\alpha_A}\,r)}{r}.
$$

This single identity is *why* `verf` (the matrix elements of the Gaussian potential) and
`rho_core` (the Gaussian charge on the grid) describe the same physics — one in the AO
basis, one on the plane-wave grid.

**Where the numbers come from** (CP2K, all-electron `POTENTIAL ALL`,
[`set_default_all_potential`, external_potential_types.F:2377-2385](../src/subsys/external_potential_types.F#L2377)):

```fortran
r     = clamp( 0.5 * covalent_radius(Z),  0.2, 1.0 )   ! core_charge_radius  (bohr)
alpha = 1 / (2 r^2)                                     ! alpha_core_charge   = α
ccore = zeff * sqrt( (alpha/pi)^3 )                     ! ccore_charge        = N
```

so in the `verfc` argument list (§2): `alpha_c` $=\alpha_A$ (cloud width),
`zeta_c` $=Z_A$ (`zeff`, the effective nuclear charge), `core_charge` $=N_A$ (the
normalization that turns `verf` into $Z_A\operatorname{erf}/r$), and `core_radius`
$=r$ (a screening cutoff). A *smaller* $r$ (sharper cloud, larger $\alpha$) means the
erf potential follows $1/r$ in closer to the nucleus — i.e. less of the operator is left
for the analytic erfc remainder. The choice $r\sim$ half the covalent radius is a
pragmatic compromise: tight enough to be physically reasonable, soft enough that
`rho_core` is representable on a normal plane-wave grid.

**Connection to range separation / Ewald.** This is the same mathematics as Ewald
summation and range-separated (screened) hybrid functionals: a long-range smooth part
($\operatorname{erf}$) handled in reciprocal space / on a grid, and a short-range part
($\operatorname{erfc}$) handled analytically in real space. CP2K applies it to the
electron–nucleus attraction rather than to the electron–electron interaction, but the
splitting parameter $\alpha$ plays the identical role of "where do we hand the problem
from real space to reciprocal space."

---

## 2. The integral engine: `verfc`

All three `build_…` routines delegate the actual integral evaluation to one shared
subroutine, [`verfc`](../src/aobasis/ai_verfc.F#L138). Understanding its three output
blocks is the key to understanding the differences between the drivers.

### 2.1 What the header says

The module header of [`ai_verfc.F`](../src/aobasis/ai_verfc.F#L8-L39) states the goal
directly: build $\langle a|-Z\,\operatorname{erfc}(r)/r|b\rangle$ by writing it as

$$
-Z\Big\langle a\Big|\frac{\operatorname{erfc}(r)}{r}\Big|b\Big\rangle
= -Z\Big\langle a\Big|\frac1r\Big|b\Big\rangle
\;+\; Z\,N\,\langle ab\,\|\,c\rangle ,
$$

i.e. the erfc nuclear attraction is rewritten as a **bare $1/r$ nuclear-attraction
integral** *minus* a **three-center Coulomb integral** $\langle ab\|c\rangle$ in which
$c$ is a single normalized **s-type Gaussian sitting on the nucleus** (multiplied by a
normalization factor $N$). That three-center Coulomb integral *is* the erf (Gaussian)
half. So a single routine naturally produces all three kernels at once.

> **Cited literature (in the header):** the recursion used to climb from
> $\langle s|\cdots|s\rangle$ up to arbitrary angular momentum is the
> **Obara–Saika** scheme: *S. Obara and A. Saika, J. Chem. Phys. **84**, 3963 (1986).*

### 2.2 The three returned blocks

`verfc` fills three Cartesian-primitive arrays. They are the physically meaningful
output and the reason the drivers differ:

| array | symbol | operator represented | line |
|---|---|---|---|
| `vnuc` | $\langle a\lvert 1/r_C\rvert b\rangle$ | **bare point-nucleus** $1/r$ (charge factored out) | [ai_verfc.F:287](../src/aobasis/ai_verfc.F#L287) |
| `verf` | $\langle a\lvert\, \text{erf-Coulomb of s-Gaussian } c\,\rvert b\rangle$ | **smeared / Gaussian** $\operatorname{erf}(\sqrt\alpha r)/r$ | [ai_verfc.F:299](../src/aobasis/ai_verfc.F#L299) |
| `vabc` | combination below | **screened** $-Z\,\operatorname{erfc}(\sqrt\alpha r)/r$ | [ai_verfc.F:911](../src/aobasis/ai_verfc.F#L911) |

The final assembly at [ai_verfc.F:909-914](../src/aobasis/ai_verfc.F#L909) is literally
the partition of §1.1 in code:

```fortran
vabc(na+i, nb+j) = vabc(na+i, nb+j) - zc*vnuc(i,j,1) + cerf*verf(i,j,1)
!                                      \___________/   \____________/
!                                       -Z * (1/r)      +Z * erf/r
!                                      = -Z*(1/r) + Z*erf/r = -Z*erfc/r
```

So **a driver picks its physics by choosing which block it keeps**:

- keep `vabc` → you get the **short-range erfc** nuclear attraction;
- keep `verf` (scaled by the nuclear charge) → you get the **smooth erf / Gaussian**
  nuclear attraction;
- `vnuc` is the bare $1/r$ building block, usually only an intermediate.

### 2.3 The meaning of the prefactors (where the physics lives)

Inside `verfc`, for each pair of primitive Gaussians $a$ (exponent $\zeta_a$) and $b$
(exponent $\zeta_b$), and an s-Gaussian $c$ on the nucleus (exponent
$\zeta_c=\alpha$), the seed integrals are built from
[ai_verfc.F:259-300](../src/aobasis/ai_verfc.F#L259):

- `zetp = 1/(ζa+ζb)`, `zetw = 1/(ζa+ζb+ζc)` — combined exponents from the **Gaussian
  product theorem**.
- `f0 = exp(-ζa·ζb/(ζa+ζb)·R_ab²)` — the **overlap prefactor** of the bra–ket
  Gaussian product (how much $a$ and $b$ overlap).
- `fnuc = 2π·zetp·f0` seeds the bare nuclear attraction
  $\langle s|1/r|s\rangle = \texttt{fnuc}\cdot F_0(t)$, with $t=R_{CP}^2/\texttt{zetp}$.
- `ferf = 2\sqrt{\pi^5\,\texttt{zetw}}\cdot\texttt{zetp}\cdot\texttt{zetq}\cdot f0$
  seeds the three-center Coulomb $\langle ss\|s\rangle$.
- `f(n) = F_n(t)` is the **Boys function** (incomplete Gamma function for Gaussian
  integrals), evaluated by [`fgamma`](../src/common/gamma.F#L153). It is the universal
  special function that resolves the $1/r$ singularity for Gaussian charge
  distributions [Boys 1950; McMurchie–Davidson 1978, cited in `gamma.F`].

From these s-function seeds, the long blocks of `verfc`
([lines 305-887](../src/aobasis/ai_verfc.F#L305)) are nothing but the **Obara–Saika
vertical and horizontal recurrences** that raise angular momentum:
$[s|\cdots|s]\to[a|\cdots|s]\to[a|\cdots|b]$. Each `coset(...)` index is a Cartesian
component $(l_x,l_y,l_z)$; the recurrences shift one unit of angular momentum at a time.
This is mechanical bookkeeping of one formula — Obara–Saika eq. for three-center
nuclear-attraction / Coulomb integrals.

### 2.4 The optional relativistic block (DKH `pVp`)

`verfc` has optional arguments (`pVp`, `pVp_sum`, `vnabc`, `dkh_erfc`) used only for
**scalar-relativistic Douglas–Kroll–Hess (DKH)** all-electron calculations
([ai_verfc.F:927-1012](../src/aobasis/ai_verfc.F#L927)). When `do_dkh` is set it also
assembles $\langle \mathbf\nabla a\,|V|\,\mathbf\nabla b\rangle$ matrices
($\langle p|V|p\rangle$), needed by the DKH kinetic-energy transformation. The
identity used,
$[\,p_a|V|p_b\,]=4\zeta_a\zeta_b[a{+}1|V|b{+}1]-\dots$, is the standard rule for the
derivative of a Gaussian (raising/lowering $l$). The RIDFT and standard all-electron
paths do **not** activate this block; mention it only so the extra arguments are not
mysterious.

### 2.5 The Boys function: how Gaussians tame the $1/r$ singularity

The deep reason all of this works analytically is that the Coulomb integral of two
Gaussians has a **closed form** built on one special function. Two ingredients:

**(i) The Gaussian product theorem.** A product of two Gaussians on centers $A$ and $B$
is a single Gaussian on a new center $P$ (the exponent-weighted midpoint):

$$
e^{-\zeta_a|\mathbf r-\mathbf A|^2}\,e^{-\zeta_b|\mathbf r-\mathbf B|^2}
= \underbrace{e^{-\frac{\zeta_a\zeta_b}{\zeta_a+\zeta_b}|\mathbf A-\mathbf B|^2}}_{=\,f0\ \text{in the code}}\;
e^{-(\zeta_a+\zeta_b)|\mathbf r-\mathbf P|^2}.
$$

The exponential prefactor is exactly `f0` ([ai_verfc.F:267](../src/aobasis/ai_verfc.F#L267)):
it measures the bra–ket overlap and makes integrals between far-apart Gaussians
negligible — the basis of the **distance screening** that lets the neighbor-list loops
skip most atom pairs.

**(ii) The Boys function.** Writing $1/r$ as a Gaussian integral
$\frac1r = \frac{2}{\sqrt\pi}\int_0^\infty e^{-r^2 u^2}\,du$ turns the singular Coulomb
integral over Gaussians into an ordinary Gaussian integral, leaving behind the
**Boys function**

$$
F_n(t) = \int_0^1 s^{2n}\,e^{-t s^2}\,ds ,
$$

evaluated by [`fgamma`](../src/common/gamma.F#L153) and returned in the array `f`. Its
argument $t = R_{CP}^2/\texttt{zetp}$ is the (scaled) squared distance from the Gaussian
product center $P$ to the nucleus $C$. Two limits make its role intuitive:

- $t\to 0$ (nucleus sits on the charge): $F_n(0)=1/(2n+1)$ — **finite**, so the
  smeared-charge integral has no singularity to diverge.
- $t\to\infty$ (nucleus far away): $F_n(t)\sim \tfrac12\Gamma(n+\tfrac12)\,t^{-(n+1/2)}$ —
  decays, recovering the expected $1/R$ falloff of a distant charge.

The *same* $F_n$ supplies both kernels in `verfc`; only its argument differs — `vnuc`
uses $t=R_{CP}^2/\texttt{zetp}$ (point nucleus), `verf` uses
$t=-f4\cdot R_{CP}^2/\texttt{zetp}$ with $f4=-\zeta_c\,\texttt{zetw}$, which folds in the
finite width $\zeta_c=\alpha$ of the nuclear Gaussian. That one substitution is the
entire difference between "point charge" and "Gaussian charge" at the integral level.
The higher orders $F_n$, $n>0$, are precisely what the Obara–Saika recurrences consume
to build up angular momentum (each step trades one unit of $l$ for one order of $n$).

---

## 3. The three `build_…` drivers compared

All three share the same skeleton (it is worth recognizing the pattern once):

1. Build a per-kind list of orbital basis sets.
2. Loop over orbital–orbital neighbor pairs `sab_orb` (gives $\mu$ on atom A, $\nu$ on
   atom B).
3. For each pair, loop over **all nuclei** via the orbital–nucleus neighbor list
   `sac_ae`, call `verfc`, and accumulate primitive integrals into `hab`.
4. Contract primitives → spherical and scatter into the target DBCSR matrix.

They differ in **(a) which `verfc` block they keep, (b) what target matrix they write,
(c) whether they compute forces/energies, and (d) threading.** The table summarizes;
the subsections give detail.

| | `build_core_ae` | `build_erfc` | `build_core_erf_ridft` |
|---|---|---|---|
| File / line | [core_ae.F:90](../src/core_ae.F#L90) | [core_ae.F:928](../src/core_ae.F#L928) | [core_ae.F:539](../src/core_ae.F#L539) |
| Block kept | `vabc` = $-Z\,\operatorname{erfc}/r$ | `vabc` = $-Z\,\operatorname{erfc}/r$ | `verf` = $+\,\text{erf}/r$ |
| Physics | short-range (sharp) nuclear attraction | same, for analysis | smooth Gaussian nuclear attraction |
| Writes into | core Hamiltonian `matrix_h` | a scratch energy matrix `e_mat` | KS matrix `matrix_ks` |
| $\alpha,Z,N$ from | the atomic potential object | **caller-supplied arrays** `calpha,ccore` | the atomic potential object |
| Forces / virial | **yes** (`verfc_force`) | no | no |
| Energy returned | atomic decomposition `atcore` (optional) | `atcore` (optional) | $E_{\text{nuc-erf}}$ **and** $E_{\text{nn-erf}}$ |
| OpenMP | yes (locks) | yes (locks) | no (serial) |
| Caller | `qs_core_matrices.F:295` | `ed_analysis.F:1015` | `qs_ks_methods.F:1376` |
| Role | standard GAPW/all-electron core H | energy-decomposition analysis | **RIDFT grid-free Hartree** |

### 3.1 `build_core_ae` — the standard short-range nuclear attraction

[core_ae.F:90](../src/core_ae.F#L90). This is the production routine that puts the
**sharp erfc** part of $-Z/r$ into the core Hamiltonian for GAPW / all-electron runs.
It is called from [`qs_core_matrices.F:295`](../src/qs_core_matrices.F#L295).

- Keeps `vabc` (the `hab` array it passes to `verfc` *is* `vabc`), then contracts and
  adds it to `matrix_h` ([core_ae.F:444-455](../src/core_ae.F#L444)).
- The screening parameters per nucleus — `alpha_c` ($\alpha$), `zeta_c` ($Z$),
  `core_charge` ($N$), `core_radius` — come from the atom's potential object via
  `get_potential` ([core_ae.F:334-340](../src/core_ae.F#L334)).
- **Forces & virial:** when `calculate_forces` is set it calls `verfc` with one extra
  unit of angular momentum (`la_max+nder`) to obtain the derivative integrals
  `habd` (= `vabc_plus`), then contracts them with the density in
  [`verfc_force`](../src/core_ae.F#L844). Translational invariance gives the force on
  the nucleus $C$ as $-(\mathbf f_a+\mathbf f_b)$
  ([core_ae.F:410-412](../src/core_ae.F#L410)).
- **Atomic energy decomposition:** if `atcore` is present it accumulates
  $\operatorname{Tr}[\rho\,V]$ per atom ([core_ae.F:426-429](../src/core_ae.F#L426)).
- **Threaded** with per-block OpenMP locks.

> The complementary smooth half for `build_core_ae` is supplied elsewhere: in GAPW it
> is the Poisson solution of `rho_core` on the grid.

### 3.2 `build_erfc` — the same integral, repackaged for analysis

[core_ae.F:928](../src/core_ae.F#L928). Functionally a **stripped-down twin** of
`build_core_ae`: it also keeps `vabc` = $-Z\,\operatorname{erfc}(\alpha r)/r$, but

- the screening parameters are passed in explicitly as arrays `calpha`, `ccore`
  ([core_ae.F:936](../src/core_ae.F#L936)) instead of read from the potential object,
  and `core_radius` is derived locally from `exp_radius`
  ([core_ae.F:1148](../src/core_ae.F#L1148));
- it has **no forces, no virial** — only the optional per-atom energy `atcore`;
- it writes into a caller-provided matrix `e_mat`, not the core Hamiltonian.

It is used by the **energy-decomposition analysis** in
[`ed_analysis.F:1015`](../src/ed_analysis.F#L1015), where one wants the erfc nuclear
term evaluated with custom $\alpha/Z$ values for a chosen partitioning. Think of it as
“`build_core_ae` minus forces, with hand-fed parameters.”

### 3.3 `build_core_erf_ridft` — the RIDFT grid-free Hartree term

[core_ae.F:539](../src/core_ae.F#L539). This is the **new** routine and the mirror
image of `build_core_ae`: it keeps the **other** half of the partition, `verf`, i.e.
the **smooth Gaussian** nuclear attraction $-Z\,\operatorname{erf}(\sqrt\alpha r)/r$.
Called from [`qs_ks_methods.F:1376`](../src/qs_ks_methods.F#L1376).

Why this exists: RIDFT deliberately **skips the plane-wave grid**, so it never solves
the Poisson equation for `rho_core`. That smooth erf contribution would simply be
missing. `build_core_erf_ridft` reconstructs it analytically in the AO basis. It
produces **two** quantities:

1. **Electron–nucleus erf potential → KS matrix.** It accumulates
   `hab -= core_charge * verf` ([core_ae.F:735](../src/core_ae.F#L735)) into a scratch
   DBCSR matrix `v_nuc_erf` cloned from the KS matrix
   ([core_ae.F:620-624](../src/core_ae.F#L620)). After the loop it does, per image,
   ([core_ae.F:777-784](../src/core_ae.F#L777)):
   - `dbcsr_dot(rho_ao, v_nuc_erf)` → the energy
     $E_{\text{nuc-erf}}=\operatorname{Tr}[\rho\,V_{\text{nuc-erf}}]$ (the symmetric
     off-diagonal factor of 2 and the MPI reduction are applied **inside** `dbcsr_dot`,
     so the caller must not reduce again — see [[project_build_core_erf_ridft]]);
   - `dbcsr_add(matrix_ks, v_nuc_erf)` → fold $V_{\text{nuc-erf}}$ into the KS matrix.

   Keeping a separate scratch matrix (rather than the per-block manual trace used by
   `build_core_ae`) is what lets one `dbcsr_dot` handle the symmetry factor and the MPI
   sum cleanly.

2. **Nucleus–nucleus erf repulsion (scalar).** The Gaussian nuclei also repel each
   other; this is the nuclear part of the GAPW grid Hartree energy. It is computed in a
   plain $O(N^2)$ atom-pair loop ([core_ae.F:807-820](../src/core_ae.F#L807)):

   $$
   E_{\text{nn-erf}}=\sum_{A<B}\frac{Z_A Z_B\,\operatorname{erf}(\sqrt{\alpha_{AB}}\,R_{AB})}{R_{AB}},
   \qquad \alpha_{AB}=\frac{\alpha_A\alpha_B}{\alpha_A+\alpha_B}.
   $$

   This loop is deliberately **not screened** — `erf(√α R)/R` does not decay, so every
   pair contributes. The reduced exponent $\alpha_{AB}$ is the standard result for the
   Coulomb interaction of two Gaussians.

Differences from `build_core_ae` to note explicitly:

- **Opposite block:** `verf` (smooth) vs `vabc` (sharp). The `vabc_scratch` array here
  exists only to absorb and discard `verfc`'s erfc output
  ([core_ae.F:727](../src/core_ae.F#L727)).
- **No forces** (energy/SCF path only), and currently **serial** (`nthread=1`,
  [core_ae.F:607](../src/core_ae.F#L607)).
- Returns explicit scalar energies `e_nuc_erf`, `e_nn_erf` rather than writing the core
  Hamiltonian.

### 3.4 `verfc_force` — the gradient helper

[core_ae.F:844](../src/core_ae.F#L844). Not a driver but the companion of
`build_core_ae`. Given the “plus” derivative integrals `habd` (angular momentum
$l{+}1$) and the density block `pab`, it assembles the Hellmann–Feynman-type gradients
$\mathbf f_a,\mathbf f_b$ using the Gaussian derivative relations
$\partial_i e^{-\zeta r^2}\to$ raise/lower $l_i$ with weight $2\zeta$ or $l_i$
([core_ae.F:890-905](../src/core_ae.F#L890)). The force on the nucleus follows from
translational invariance in the caller.

---

## 4. From a primitive integral to the Kohn–Sham matrix

Sections 2–3 covered the *integrals*. This section covers the *implementation machinery*
that turns those integrals into actual matrix elements inside CP2K's data structures —
the part every reader of `core_ae.F` has to recognize, and the context that defines what
"core contribution" means in the code.

### 4.1 The shared assembly pipeline

Every `build_…` routine runs the same four-stage pipeline. Knowing it once lets you read
all of them:

1. **Two nested neighbor lists.** The outer loop walks `sab_orb` — pairs of *orbital*
   basis functions $(\mu\,\text{on atom }A,\ \nu\,\text{on atom }B)$ that overlap. The
   inner loop walks `sac_ae` — for the current bra atom, the list of *nuclei* $C$ close
   enough to contribute. So each `verfc` call evaluates a genuine **three-center**
   quantity $(\mu_A,\nu_B \mid C)$. The neighbor lists are pre-screened by radius, so
   distant, negligible triples are never visited
   ([core_ae.F:359-372](../src/core_ae.F#L359)).

2. **Distance screening.** Inside the loops, `set_radius_a + core_radius < dac`-type
   tests skip pairs whose Gaussian overlap (the `f0` factor of §2.5) is below threshold.
   This is what makes the formally $O(N_{\text{atom}}^3)$ nuclear-attraction build scale
   *linearly* for large systems.

3. **Contraction: primitive → contracted → spherical.** `verfc` returns integrals over
   **primitive Cartesian** Gaussians. They are contracted to the **spherical** AO basis
   actually used by the SCF, by two matrix multiplies with the contraction coefficients
   `sphi_a`, `sphi_b` ([core_ae.F:444-455](../src/core_ae.F#L444)):
   $H_{\text{sph}} = \texttt{sphi\_a}^{\!\top}\,(\text{primitive } H)\,\texttt{sphi\_b}$.
   This is where the cheap, dense primitive work is folded into the compact basis.

4. **Scatter into DBCSR with symmetric storage.** The result is added into a sparse
   distributed block-matrix (DBCSR). Because $H$ is symmetric, only the
   $\text{irow}\le\text{icol}$ triangle is stored; the `IF (iatom <= jatom)` /
   transpose branch present in all three routines
   ([core_ae.F:447-455](../src/core_ae.F#L447)) writes into the canonical triangle.
   OpenMP per-block locks guard the concurrent writes in `build_core_ae`/`build_erfc`.

The `f0 = 2.0` factor for off-diagonal atom pairs ([core_ae.F:286](../src/core_ae.F#L286))
and the `dbcsr_dot` energy traces (§5.3) both exist to account for this stored-triangle
convention: an off-diagonal block represents *two* symmetric matrix elements.

### 4.2 What "core" means: the core-Hamiltonian family

In CP2K the **core Hamiltonian** $H_{\text{core}}$ is the density-independent part of
the Kohn–Sham operator — everything that does *not* change as the SCF density updates:

$$
H_{\text{core}} \;=\; \underbrace{T}_{\text{kinetic}} \;+\; \underbrace{V_{\text{en}}}_{\text{electron–nucleus}} \;(+\,V_{\text{nl}}\ \text{for pseudopotentials}).
$$

It is assembled **once per geometry** in
[`qs_core_matrices.F`](../src/qs_core_matrices.F) and then reused unchanged in every SCF
iteration (only the kinetic and overlap matrices and these nuclear terms — never the
Hartree or XC potentials, which are rebuilt each step). `build_core_ae` is one member of
a small family that fills the nuclear part, selected by how the nucleus is modeled:

| nucleus model | routine | what it adds to $H_{\text{core}}$ | line |
|---|---|---|---|
| **all-electron** (`POTENTIAL ALL`) | **`build_core_ae`** | $-Z\,\operatorname{erfc}(\sqrt\alpha r)/r$ (sharp nuclear attraction) | [qs_core_matrices.F:295](../src/qs_core_matrices.F#L295) |
| pseudopotential, **local** part | `build_core_ppl` | GTH local PP: erf core + short-range Gaussian×polynomial | [qs_core_matrices.F:328](../src/qs_core_matrices.F#L328) |
| pseudopotential, **non-local** part | `build_core_ppnl` | separable projector $\sum |p\rangle h\langle p|$ | [qs_core_matrices.F:359](../src/qs_core_matrices.F#L359) |

This is the precise sense in which `build_core_ae` is a "core contribution": it is the
all-electron nuclear-attraction entry of the core Hamiltonian, the analogue of the GTH
pseudopotential routines for the case where no pseudopotential is used. Its *smooth*
counterpart (the erf half) is deliberately **not** here — it lives in the
density-dependent Hartree term (next), because on the grid it is naturally lumped with
the electron density.

### 4.3 Two different places in the SCF cycle

The split of the nuclear attraction across the static $H_{\text{core}}$ and the dynamic
Hartree term shows up as two different *call sites* in the code:

- **Built once, before SCF** — `build_core_ae` (the erfc part), from
  `qs_core_matrices.F`. Density-independent; the matrix is stored and added into the KS
  matrix every iteration for free.
- **Built every SCF step, inside the KS build** — `build_core_erf_ridft` (the erf part,
  RIDFT only), from [`qs_ks_methods.F:1376`](../src/qs_ks_methods.F#L1376), right next to
  the RI electron–electron Hartree build. It *is* density-dependent (the energy trace
  $\operatorname{Tr}[\rho\,V_{\text{nuc-erf}}]$ uses the current $\rho$), so it cannot be
  precomputed.

In GAPW the erf part is instead produced inside `sum_up_and_integrate` /
`qs_ks_utils.F`, where the `rho_core` Gaussian is added to the total density and the grid
Poisson solve is performed. RIDFT replaces that grid step with the analytic
`build_core_erf_ridft` call — the v_rspace Hartree potential is **zeroed for RIDFT** so
the grid contribution is not double-counted ([CLAUDE.md project notes]). This is the
single concrete code change that distinguishes the two methods' Hartree channel.

---

## 5. Why these contributions exist: the SCF energy framework

The previous section showed *how* the integrals are split. This section explains
*why* — what each piece means inside the Kohn–Sham SCF cycle, which energy term it
feeds, and exactly how the bookkeeping differs between GAPW and RIDFT. Concrete numbers
are the converged H₂ test (`s_one_exponent` orbital basis, `aug-cc-pVTZ-RIFIT` RI basis,
`XC NONE`) from [[project_ridft_larger_ri_basis_comparison]].

### 5.1 The Kohn–Sham energy and where the nucleus enters

The total DFT energy that CP2K minimizes is

$$
E_{\text{tot}}[\rho] \;=\; \underbrace{T_s[\rho] + E_{\text{en}}[\rho]}_{\text{core Hamiltonian}} \;+\; \underbrace{E_{\text{H}}[\rho]}_{\text{Hartree}} \;+\; E_{\text{xc}}[\rho] \;+\; E_{\text{nn}},
$$

with $T_s$ the kinetic energy, $E_{\text{en}}$ the electron–nucleus attraction,
$E_{\text{H}}$ the electron–electron Coulomb (Hartree) energy, $E_{\text{xc}}$ the
exchange–correlation energy, and $E_{\text{nn}}$ the nucleus–nucleus repulsion.

The SCF loop builds the Kohn–Sham matrix
$F = H_{\text{core}} + J[\rho] + V_{\text{xc}}[\rho]$, diagonalizes it for a new density
$\rho$, and repeats until self-consistent. The two matrix pieces have very different
character, and this is the reason the nuclear attraction is **deliberately split across
both of them**:

| matrix piece | density-dependent? | rebuilt each SCF step? | nuclear part it carries |
|---|---|---|---|
| **$H_{\text{core}}$** (core Hamiltonian) | no | no — built **once** | the **sharp erfc** electron–nucleus attraction (`build_core_ae`) |
| **$J[\rho]$** (Hartree / Coulomb) | yes | yes — **every** step | the **smooth erf** electron–nucleus attraction + e–e + n–n |

So the single physical operator $-Z/r$ is divided **not for physics but for numerics**:
its smooth half travels with the (equally smooth) electron density into the
density-dependent Hartree term, where a grid or an RI fit can handle them together;
its sharp half is frozen into the static core Hamiltonian, where it is integrated
analytically once and never revisited. The erf/erfc partition of §1.1 is precisely the
seam along which the operator is cut to make each half live comfortably in its matrix.

### 5.2 What energy goes into which term

Mapping the physical interactions onto CP2K's printed energy components:

| physical interaction | kernel | CP2K energy term | how evaluated (GAPW) | how evaluated (RIDFT) |
|---|---|---|---|---|
| electron kinetic | — | `energy%core` (part) | analytic | analytic *(identical)* |
| electron–nucleus, **short range** | $-Z\,\operatorname{erfc}/r$ | `energy%core` (part) | `build_core_ae` | `build_core_ae` *(identical)* |
| electron–nucleus, **smooth** | $-Z\,\operatorname{erf}/r$ | `energy%hartree` | grid: Poisson of `rho_core` | **`build_core_erf_ridft`** → `e_nuc_erf` |
| electron–electron | $1/r$ | `energy%hartree` | grid: Poisson of `rho_elec` | **RI-V** → `E_ee^RI` |
| nucleus–nucleus, **smooth** | $\operatorname{erf}/r$ | `energy%hartree` | grid: `rho_core`–`rho_core` | **`build_core_erf_ridft`** → `e_nn_erf` |
| nucleus–nucleus, **short range** | $\operatorname{erfc}/r$ | `energy%core_overlap` | analytic pair sum | analytic pair sum *(identical)* |
| Gaussian-core self-energy correction | — | `energy%core_self` | analytic, removes grid self-interaction | folded into `energy%hartree` (see §5.4) |
| exchange–correlation | — | `energy%exc` | grid | grid *(identical)* |

The crucial observation: **everything outside `energy%hartree` is bit-identical between
GAPW and RIDFT** — same kinetic, same erfc core Hamiltonian (`build_core_ae` is used by
both), same short-range core overlap, same XC. RIDFT changes **only** how the Hartree
term is assembled. That is the entire surface of the method.

### 5.3 The Hartree term, side by side

**GAPW** computes the whole Hartree energy in one grid Poisson solve of the *total*
charge $\rho_{\text{tot}} = \rho_{\text{elec}} + \rho_{\text{core}}$ (electrons plus
Gaussian nuclei). Expanding that single $\tfrac12\langle\rho_{\text{tot}}|V_{\text{tot}}\rangle$:

$$
E_{\text{H}}^{\text{GAPW}} = \underbrace{E_{ee}}_{\text{e–e}} + \underbrace{E_{en}^{\text{erf}}}_{\text{e–n, smooth}} + \underbrace{E_{nn}^{\text{erf}}}_{\text{n–n, smooth}} - \underbrace{E_{\text{self}}}_{\text{remove spurious self-interaction}} .
$$

Because the Poisson solve treats each Gaussian nucleus as a charge distribution, it
spuriously includes each nucleus interacting with *itself*; `energy%core_self`
($-\sum_A \tfrac12 Z_A^2\sqrt{2\alpha_A/\pi}$) is the analytic correction that removes it.

**RIDFT** assembles the *same four terms* without ever touching the grid
([qs_ks_methods.F:1382](../src/qs_ks_methods.F#L1382)):

```fortran
energy%hartree = energy%hartree        & ! E_ee^RI   : electron-electron, from the RI-V fit
               + e_nuc_erf             & ! E_en^erf  : Tr[P*V_nuc_erf], from build_core_erf_ridft
               + e_nn_erf              & ! E_nn^erf  : Z_A*Z_B*erf(...)/R, the non-screened pair sum
               - energy%core_self        ! cancels the bookkeeping self-energy (see 4.4)
```

For the H₂ test the two reconstructions agree to 14.7 mHa
([[project_ridft_larger_ri_basis_comparison]]):

| Hartree sub-term | value (Ha) | source |
|---|---|---|
| $E_{ee}^{\text{RI}}$ | $+1.5731$ | RI-V fit (`get_hartree_noenv`) |
| $E_{en}^{\text{erf}}$ = `e_nuc_erf` | $-3.8306$ | `build_core_erf_ridft`, `verf` block |
| $E_{nn}^{\text{erf}}$ = `e_nn_erf` | $+0.7151$ | `build_core_erf_ridft`, atom-pair loop |
| $-E_{\text{self}}$ = `-core_self` | $+2.8209$ | analytic Gaussian self-energy |
| **RIDFT `energy%hartree`** | **$1.2785$** | sum of the above |
| **GAPW `energy%hartree`** | **$1.2932$** | one grid Poisson solve |

The residual 14.7 mHa is **not** a nuclear-term error — `e_nuc_erf` and `e_nn_erf` are
analytic and effectively exact. It sits entirely in $E_{ee}^{\text{RI}}$: the RI-V fit
($1.5731$) underestimates the grid e–e energy ($\approx1.588$) because the auxiliary
basis cannot fully represent the H–H bond-midpoint density. Coulomb-metric RI-V is
variational, so it is a strict *lower bound* — consistent with RIDFT's Hartree sitting
below GAPW's. This is intrinsic RI approximation error, not a bug.

### 5.4 The self-energy subtraction — a subtle but exact cancellation

CP2K's total energy is the plain sum of its components
([qs_ks_methods.F](../src/qs_ks_methods.F)):

$$
E_{\text{tot}} = \texttt{core} + \texttt{core\_overlap} + \texttt{core\_self} + \texttt{exc} + \texttt{hartree},
$$

where `core_self` is already negative. For the H₂ test (GAPW):

$$
1.42447 + 0.0000005 + (-2.82095) + 0 + 1.29316 = -0.10332\ \text{Ha} \;\checkmark
$$

GAPW *needs* the `core_self` line because its grid Poisson over-counts the Gaussian
self-interaction. **RIDFT's nuclear–nuclear term `e_nn_erf` is a clean $\sum_{A<B}$ pair
sum that contains no self-interaction at all** — so the `core_self` correction must be
removed. That is exactly what the `- energy%core_self` inside the RIDFT Hartree
assignment does: it cancels, term-for-term, the `+ core_self` that the total-energy sum
will add back. Net effect — RIDFT carries the nuclear repulsion as the explicit analytic
`e_nn_erf`, while GAPW carries it implicitly through `rho_core`-Poisson plus the
`core_self` correction. Both reach the same total
([qs_ks_methods.F:1382](../src/qs_ks_methods.F#L1382)):

$$
\text{RIDFT:}\quad 1.42447 + 0.0000005 + (-2.82095) + 0 + 1.27850 = -0.11798\ \text{Ha} \;\checkmark
$$

(the $-2.82095$ in the total and the $+2.82095$ folded into `hartree` cancel, leaving
the analytic `e_nn_erf` as the true nuclear repulsion).

### 5.5 Summary of the comparison

- **Identical by construction:** kinetic energy, the erfc core-Hamiltonian nuclear
  attraction (`build_core_ae`), short-range core overlap, XC, and the converged density
  (for a symmetry-fixed system). RIDFT touches **none** of these.
- **The only difference is the Hartree channel.** GAPW does one grid Poisson solve over
  electrons + Gaussian nuclei; RIDFT replaces it with (i) an **RI-V fit** for the
  electron–electron part and (ii) **`build_core_erf_ridft`** for the smooth
  electron–nucleus and nucleus–nucleus parts, with `core_self` rebalanced.
- **The motivation:** removing the plane-wave grid (and its FFT, cutoff, and global
  communication) from the Hartree build. The price is the RI fitting error in the e–e
  term, bounded below and controllable by the auxiliary basis quality.
- **The verification target** (CLAUDE.md) — RIDFT KS ≡ GAPW KS — is therefore reachable
  only *within RI accuracy*: the nuclear pieces match analytically; the e–e piece matches
  only as well as the RI auxiliary basis spans the orbital-product density.

### 5.6 GAPW's grid machinery vs RIDFT's RI machinery — in detail

The summary above says "GAPW does a Poisson solve, RIDFT does RI." This subsection
unpacks both sides, because the contrast is the whole point of the method.

**What GAPW actually does (the soft/hard decomposition).** GAPW
[Lippert 1999; Krack–Parrinello] does not put the true, cuspy density on the grid. It
writes the density as a *soft* part plus, per atom, a *hard* local correction:

$$
\rho(\mathbf r) = \tilde\rho(\mathbf r) \;+\; \sum_A\big(\rho^1_A(\mathbf r) - \tilde\rho^1_A(\mathbf r)\big),
$$

where $\tilde\rho$ is smooth (grid-representable) and the on-atom pair $\rho^1_A,
\tilde\rho^1_A$ (hard/soft one-center densities) live on dense **atomic radial grids**.
The Hartree energy then splits the same way:

- the **soft** part $\tilde\rho + \rho_{\text{core}}$ goes through the plane-wave
  **Poisson solve** → `energy%hartree`;
- the **one-center** parts go through analytic atomic integrals
  (`Vh_1c_gg_integrals`) → `energy%hartree_1c`.

So a *full* GAPW Hartree energy is `energy%hartree + energy%hartree_1c`. RIDFT, by
contrast, never forms the one-center correction: it computes one analytic
electron–electron Hartree from the RI fit and leaves `energy%hartree_1c = 0`. **When
comparing the two methods you must compare RIDFT's `energy%hartree` against GAPW's
`energy%hartree + energy%hartree_1c`**, not against `energy%hartree` alone — see the
CLAUDE.md project notes and [[project_ridft_hartree_energy]].

**What RIDFT does instead (resolution of identity, RI-V).** RIDFT expands every orbital
*product density* in an atom-centered auxiliary (RI) basis $\{P\}$:

$$
\rho_{\mu\nu}(\mathbf r) = \phi_\mu(\mathbf r)\phi_\nu(\mathbf r) \;\approx\; \sum_P d^{\,\mu\nu}_P\,\chi_P(\mathbf r),
$$

and obtains the fit coefficients by minimizing the **Coulomb-metric** error, which gives
the symmetric "RI-V" expression for the Hartree (Coulomb $J$) matrix
([rt_bse_utils.F](../src/rt_bse_utils.F), `get_hartree_noenv`):

$$
J_{\mu\nu} = \sum_{PQ} (\mu\nu\,|\,P)\,[\mathbf V^{-1}]_{PQ}\,Q_Q,
\qquad Q_Q = \sum_{\lambda\sigma}(Q\,|\,\lambda\sigma)\,\rho_{\lambda\sigma},
\qquad V_{PQ} = (P\,|\,Q),
$$

where $(\mu\nu|P)$ are the **three-center** Coulomb integrals
(`create_hartree_ri_3c`) and $V_{PQ}$ is the **two-center** Coulomb metric, inverted by
`create_2c_coulomb_rep_v`. The electron–electron Hartree energy is
$E_{ee}^{\text{RI}} = \tfrac12\operatorname{Tr}[\rho\,J]$. This is the standard density
fitting of Whitten / Dunlap / Vahtras — exact in a complete aux basis, variationally a
**lower bound** in a finite one (which is why RIDFT's Hartree sits *below* GAPW's, §5.3).

**The clean division of labor in RIDFT's Hartree.** Putting §3.3 and the RI fit together,
RIDFT builds the *entire* GAPW Hartree contribution from three analytic pieces with no
grid:

| GAPW (grid) | RIDFT (analytic) | routine |
|---|---|---|
| Poisson of $\tilde\rho$ (electron–electron) | RI-V $J$ matrix | `get_hartree_noenv` |
| Poisson of $\rho_{\text{core}}$ × $\tilde\rho$ (electron–nucleus, smooth) | `verf` integrals $\to E_{en}^{\text{erf}}$ | **`build_core_erf_ridft`** |
| Poisson of $\rho_{\text{core}}$ × $\rho_{\text{core}}$ (nucleus–nucleus, smooth) | $\sum_{A<B}Z_AZ_B\operatorname{erf}/R$ | **`build_core_erf_ridft`** |
| `core_self` correction (grid self-interaction) | rebalanced (§5.4) | `qs_ks_methods.F` |

Read top-to-bottom, that table *is* the RIDFT method: every row that GAPW evaluates by
laying charge on a plane-wave grid and solving Poisson, RIDFT evaluates by an analytic AO
or RI integral. The replacement is exact for the two nuclear rows and approximate (RI
fitting) for the electron–electron row — which is exactly where the residual 14.7 mHa of
§5.3 lives.

### 5.7 Why "differently" means only the smooth half — and why the split exists at all

Two questions naturally arise once both methods are all-electron: *if the nucleus is the
same point charge $-Z/r$ in both, why is its core contribution computed differently? And
why introduce the smooth/sharp separation in the first place?* The answer to both is
**numerics, not physics**, and it is worth stating crisply because it is easy to
over-state the difference between the methods.

**The methods agree on the sharp half exactly.** GAPW and RIDFT do *not* compute all of
the nuclear contribution differently — they compute only the **smooth (erf) half**
differently. The **sharp (erfc) half is bit-identical**: both call the same
[`build_core_ae`](../src/core_ae.F#L90) and freeze $-Z\,\operatorname{erfc}(\sqrt\alpha
r)/r$ into the static core Hamiltonian once per geometry (§4.2). So the real question is
narrow: why is the *smooth electron–nucleus* piece evaluated one way in GAPW and another
in RIDFT?

**Why the split exists at all.** The all-electron nucleus enters through one operator
with a hard singularity, $-\sum_A Z_A/|\mathbf r-\mathbf R_A|$. That $1/r$ cusp is the
entire difficulty: CP2K's plane-wave grid represents *smooth* functions, and **a sharp
$1/r$ cusp cannot be put on a finite grid** (it would demand an infinite plane-wave
cutoff). The erf/erfc identity of §1.1 cuts the operator along exactly the seam that
makes each half tractable — the smooth erf half looks like the potential of a Gaussian
cloud `rho_core` and rides with the (equally smooth) electron density into the
density-dependent Hartree term; the sharp erfc half is short-ranged and is integrated
analytically into the static core Hamiltonian. The split is therefore **the line along
which the single operator $-Z/r$ is cut so that neither half is hard to integrate**, not
a statement that different electrons feel different physics. (This is the same
range-separation mathematics as Ewald sums and screened hybrids — see §1.3.)

**Why the smooth half is then handled differently in the two methods.** Both methods put
the smooth half in the same slot — the Hartree channel — but they *build that channel
with different machinery.* RIDFT's entire purpose is to **delete the plane-wave grid**
(no FFT, no cutoff, no global communication in the Hartree build). The instant the grid
is gone, the Poisson solve is gone, and the smooth erf nuclear potential that GAPW
obtained *for free* as a by-product of that Poisson solve simply **vanishes** — electrons
would never feel the nuclei through the Hartree channel. RIDFT therefore has no choice
but to **reconstruct the identical smooth erf contribution by a different route**:
analytically, in the AO basis, from the *other* output block of the very same `verfc`
engine. GAPW keeps `verfc`'s `vabc` (erfc) block and lets the grid supply the erf half;
RIDFT keeps `verfc`'s `verf` (erf) block and computes the erf half explicitly
([`build_core_erf_ridft`](../src/core_ae.F#L539)). They are **mirror images** picking
opposite halves of one exact partition (§2.2, §3.3).

**The one-sentence summary.** *All-electron* is what forces the smooth/sharp split to
exist (the cusp must be cut to be integrable); *removing the grid* is what forces the
smooth half to be computed two different ways (Poisson in GAPW vs analytic `verf` in
RIDFT). Everything downstream of the difference is downstream of "no grid," not of
"all-electron" — had RIDFT kept the grid, `build_core_erf_ridft` would not need to exist.
The verification target is that the two routes give the same number, which they do up to
the RI fitting error in the *electron–electron* term (≈14.7 mHa for H₂, §5.3); the
nuclear pieces themselves match analytically.

---

## 6. How this fits into nuclear-contribution handling across CP2K

Stepping back, the `core_ae.F` routines are one corner of a consistent strategy. The
common thread is always the **erf/erfc Coulomb partition** of §1.1, applied differently
per method:

- **GPW with GTH pseudopotentials** (the common production setup): the nucleus + core
  electrons are replaced by a pseudopotential. Its **local long-range part** is exactly
  an erf (Gaussian `rho_core`) put on the grid via Poisson; its **local short-range
  part** (Gaussian × polynomial) and **non-local projectors** are analytic
  (`core_ppnl`). No bare $1/r$ ever appears. [Lippert 1997; VandeVondele 2005;
  Goedecker–Teter–Hutter pseudopotentials.]

- **GAPW / all-electron:** the full $-Z/r$ is kept and **split**:
  - smooth erf half → Gaussian `rho_core` on the grid → Poisson;
  - sharp erfc half → **`build_core_ae`** analytic AO integrals.
  [Lippert 1999; Krack & Parrinello.]

- **RIDFT (this branch):** removes the grid entirely. The smooth erf half that GAPW got
  from Poisson is now produced analytically by **`build_core_erf_ridft`**, while the
  sharp erfc half still comes from `build_core_ae`. Together they reconstruct the full
  GAPW nuclear contribution without any FFT. See [[concept_erf_erfc_partition]] and
  [[project_ridft_larger_ri_basis_comparison]] for the term-by-term H₂ verification.

- **Energy decomposition analysis** uses **`build_erfc`** to evaluate the erfc nuclear
  term with user-chosen screening, for interpreting interaction energies.

- **Scalar-relativistic (DKH) all-electron** runs reuse the very same `verfc` kernel,
  switching on its optional `pVp` block (§2.4).

The pleasing symmetry to take away: **every nuclear contribution in CP2K is the same
two integral families — bare $1/r$ (`vnuc`), Gaussian erf (`verf`), and their
difference erfc (`vabc`) — computed once by `verfc` and routed by each driver to
whichever part of the algorithm (grid or analytic, Hamiltonian or energy) needs it.**

---

## References

1. **S. Obara and A. Saika**, *Efficient recursive computation of molecular integrals
   over Cartesian Gaussian functions*, **J. Chem. Phys. 84, 3963 (1986)**.
   — The recursion scheme implemented in `verfc` (cited in its header,
   [ai_verfc.F:38-39](../src/aobasis/ai_verfc.F#L38)).
2. **L. E. McMurchie and E. R. Davidson**, *One- and two-electron integrals over
   Cartesian Gaussian functions*, **J. Comp. Phys. 26, 218 (1978)**.
   — The Boys / incomplete-Gamma evaluation in
   [`gamma.F`](../src/common/gamma.F#L142).
3. **S. F. Boys**, *Electronic wave functions. I.*, **Proc. R. Soc. Lond. A 200, 542
   (1950)**. — Origin of the Boys function $F_n(t)$ for Gaussian Coulomb integrals.
4. **G. Lippert, J. Hutter, M. Parrinello**, *A hybrid Gaussian and plane wave density
   functional scheme*, **Mol. Phys. 92, 477 (1997)**. — GPW.
5. **G. Lippert, J. Hutter, M. Parrinello**, *The Gaussian and augmented-plane-wave
   density functional method*, **Theor. Chem. Acc. 103, 124 (1999)**. — GAPW and the
   `rho_core` / soft-hard partition.
6. **J. VandeVondele, M. Krack, F. Mohamed, M. Parrinello, T. Chassaing, J. Hutter**,
   *Quickstep: Fast and accurate density functional calculations using a mixed Gaussian
   and plane waves approach*, **Comput. Phys. Commun. 167, 103 (2005)**.
7. **A. Wolf, M. Reiher, B. A. Hess**, *The generalized Douglas–Kroll transformation*,
   **J. Chem. Phys. 117, 9215 (2002)**. — Background for the optional `pVp` (DKH) block.
8. **S. Goedecker, M. Teter, J. Hutter**, *Separable dual-space Gaussian
   pseudopotentials*, **Phys. Rev. B 54, 1703 (1996)**. — The GTH pseudopotentials whose
   local/non-local parts are built by `build_core_ppl` / `build_core_ppnl`.
9. **J. L. Whitten**, *Coulombic potential energy integrals and approximations*,
   **J. Chem. Phys. 58, 4496 (1973)**; **B. I. Dunlap**, **J. Chem. Phys. 71, 3396
   (1979)**; **O. Vahtras, J. Almlöf, M. W. Feyereisen**, **Chem. Phys. Lett. 213, 514
   (1993)**. — The resolution-of-identity / Coulomb-metric density fitting (RI-V) used
   for the RIDFT electron–electron Hartree term (§5.6).

---

*Source files referenced:* [`src/core_ae.F`](../src/core_ae.F),
[`src/aobasis/ai_verfc.F`](../src/aobasis/ai_verfc.F),
[`src/common/gamma.F`](../src/common/gamma.F),
[`src/qs_core_matrices.F`](../src/qs_core_matrices.F),
[`src/qs_ks_methods.F`](../src/qs_ks_methods.F),
[`src/rt_bse_utils.F`](../src/rt_bse_utils.F),
[`src/subsys/external_potential_types.F`](../src/subsys/external_potential_types.F),
[`src/ed_analysis.F`](../src/ed_analysis.F).
