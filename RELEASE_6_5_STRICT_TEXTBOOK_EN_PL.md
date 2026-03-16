# RELEASE 6.5 STRICT TEXTBOOK EDITION (EN + PL)

**Version:** 6.5.37  
**Date:** 2026-03-16  
**Branch:** `main`  
**Predecessor:** Release 6.4 — Strict Textbook Edition

**Status discipline (no false pass):**
- This document is a **strict-only textbook projection** of the current repo state.
- It does **not** claim strict-core selector closure, global selector closure, global `QW-2191` discharge, or ToE closure.
- Update (`2026-03-15`): the canonical **diagonal/local lane** on the strict `n=12` carrier now exports a strict-derived
  **numeric value-instantiation** deciding the diagonal mode‑2 defect nonzero condition (`T166`) on `pair1`, and in fact
  decides nonzero defects on **all** Fourier-degenerate pairs `pair_m (m=1..5)` (`N482`, `N485`, `N487`).
- Update (`2026-03-15`): the repo exports an explicit strict-derived **mode-index assignment basis object**
  `ModeIndexAssignment_canonical_local_diagonal_strict_derived_v1` fixing axes on all `pair_m (m=1..5)` **up to residual
  sign** (`F453`, packaged by `N492`).
- Update (`2026-03-15`): residual `Z2` sign flips are **gauge-equivalent via conjugation** for the `QW-2190` SU(3)/SU(2)
  embedding audits (`P452`, `N493`), so the diagonal/local lane canonicalizes the `QW-2190` embedding **uniquely up to
  conjugation** in its declared scope (`N494`).
- Update (`2026-03-15`): moreover, for the same `QW-2190` embedding audits, the full continuous `O(2)` basis-rotation
  freedom is likewise **gauge-equivalent** (conjugation-only) (`N495`, audited by `P454`).
- Update (`2026-03-15`): independently, on the strict **Shannon element‑order reference** lane, the repo exports a full
  strict-core **mode-index assignment basis object** covering all Fourier-degenerate pairs `pair_m (m=1..5)` on `n=12`,
  cutting each `pair_m` `O(2)` family down to residual `Z2` using only the internal reference datum `r_ord`
  (`N480`, `N488`, `N496`, executed by `F454`, packaged by `N500`).
- Update (`2026-03-15`): in the declared `R1` target-slot semantics (`span{u_1,u_2}`), residual `Z2` sign flips of the
  representative vectors are gauge-irrelevant for the target slot as a span object (`N501`).
- Update (`2026-03-15`): residual `Z2` sign can be **frozen as a tracked gauge/convention layer** for the currently exported
  downstream objects where the sign is provably gauge-irrelevant (packaged by `N502`, including `alpha_12 mod π`), without implying any global
  sign-sensitive physical orientation datum.
- Update (`2026-03-16`): a theorem-level arithmetic closure shows the Shannon element‑order reference defect never vanishes:
  for every `n≥2` and every `k∈{1,…,n-1}`, `F_k(ord_{Z_n})≠0` (`N514`). Therefore the element‑order reference cross‑entropy objective cuts each
  Fourier‑degenerate `pair_m` `O(2)` family down to residual `Z2` on any `Z_n` (lane-scoped; **no** physical promotion of `n≠12`).
- Update (`2026-03-15`): the two independently exported mode-index assignment bases (diagonal/local vs Shannon ord-reference)
  are aligned on all `pair_m (m=1..5)` up to residual `Z2` sign (audit `P455`); this is a hygiene consistency check and
  does not promote any global discharge.
- Update (`2026-03-15`): the exported slot‑free sigma‑int theta‑pair supply (`F451/N489`) is axis‑aligned with both exported
  mode‑index assignment bases (diagonal/local `F453` and Shannon `F454`) on `pair1/pair2` up to residual `Z2` sign (audit `P456`,
  packaged by `N498`); this is an internal consistency check and does not lift residual sign to a physical convention.
- Update (`2026-03-15`): the repo exports a sigma‑int orientation‑slice restriction artifact for the declared residual local‑diagonal
  control pullback `M_control_residual_diag_declared = T_control^T D_local_residual T_control` (`P457`) and a **conditional**
  value instantiation of that 2×2 restriction using the `N477` rewrite and the current strict‑derived `(vpsi,g4,g6)` provider (`P458`);
  this remains below any host‑matching cancellation witness.
- Update (`2026-03-15`): moreover, the repo exports a **conditional** value instantiation of the full declared 4×4 residual local‑diagonal
  control pullback on the control basis `(c1,s1,c2,s2)` using the same `N477` rewrite and strict‑derived provider (`P459`); this exposes
  the `pair1`/`pair2` blocks and cross‑block couplings, without any promotion to host matching.
- Update (`2026-03-15`): the strict sigma‑int → residual target‑slot bridge/export‑map object is explicitly upgraded to
  `Upsilon_residual_datum_sigma_int_bridge_export_map_object_v2`, attaching the slot‑free theta‑pair outputs and the corresponding strict‑core
  `R1` inhabitant instance as map outputs (`F455`, packaged by `N499`); this is an explicit upgrade (v1 remains sign‑only) and does not imply selector closure.
- Update (`2026-03-15`): the repo exports a minimal strict-core downstream operator on `V_1=span{c1,s1}` derived only from the
  materialized sigma-int orientation direction `u_1`:
  `A_1(pair1) := |u_1><u_1|` (projector; residual `Z2` sign gauge invariant). This discharges the operator-stage bridge target of `P2`
  in declared scope (`F456`) but does **not** identify the operator with the extension-only `A_1_ext` of the `H/O` lane, does not imply
  selector closure, and does not discharge global `QW-2191`.
- Update (`2026-03-15`): derived transition data `alpha_12 := (theta_2 - theta_1) mod 2π` for `pair1/pair2` is explicitly exported
  from the strict sigma-int slot-free theta-pair supply (`F457`); this is lane-scoped and does not constitute a global selector transition/gluing object.
- Update (`2026-03-15`): the repo exports an explicit **lane-scoped pair1/pair2 chart-transport operator** `O_12` on the `n=12` carrier,
  derived only from the strict sigma-int `alpha_12` transition angle (`F461`), with projector-level sign-gauge-irrelevance packaged as a strict theorem (`N506`)
  and audited by a dedicated probe (`P465`). This is still below any global selector atlas/gluing structure and does not discharge global `QW-2191`.
- Update (`2026-03-15`): building on `O_12`, the repo exports a concrete **two-chart projector operator section** on `{pair1,pair2}`:
  `A_2(pair2) = O_12 A_1(pair1) O_12^T` (`F462`), packaged as a strict theorem (`N507`) and audited by `P466`. This remains lane-scoped and does not upgrade
  to a global selector atlas or global `QW-2191` discharge.
- Update (`2026-03-15`): the repo exports an explicit lane-scoped **two-chart selector atlas stub** on `{pair1,pair2}` with an overlap-domain declaration,
  transition data (`O_12`), and gluing data (`A_2 = O_12 A_1 O_12^T`) (`F463`). This is not a global selector atlas and does not discharge global `QW-2191`.
- Update (`2026-03-15`): the repo exports an explicit lane-scoped **three-chart** selector-atlas ingredient on `{pair1,pair2,pair3}` at projector level,
  including explicit cocycle/path-independence data on the glued projector section (`F464`), audited by `P467` and packaged as a strict theorem (`N508`).
  This remains lane-scoped, does not lift residual sign to a physical convention, does not export a global selector atlas, and does not discharge global `QW-2191`.
- Update (`2026-03-15`): the repo exports an explicit lane-scoped **five-chart** selector-atlas ingredient on `{pair1..pair5}` at projector level,
  including explicit **local** cocycle/path-independence data on adjacent triple overlaps (`F465`), audited by `P468` and packaged as a strict theorem (`N509`).
  This remains lane-scoped, does not export a global selector atlas, and does not discharge global `QW-2191`.
- Update (`2026-03-15`): the repo exports additional lane-scoped axis-only long-edge chart-transport operators (`O_14`, `O_15`, `O_25`) and an upgraded
  five-chart `{pair1..pair5}` projector-level selector-atlas ingredient with explicit cocycle/path-independence audit data for **all** triple overlaps (`F466`),
  audited by `P469` and packaged as a strict theorem (`N510`). This remains lane-scoped, does not export a global selector atlas, and does not discharge global `QW-2191`.
- Update (`2026-03-15`): the repo exports a lane-scoped **oriented transport (α mod 2π) lift** on `{pair1..pair5}` induced by the exported representative vectors
  `u_1..u_5` as a **tracked gauge/convention layer** (sign-tracked), with full triple cocycle/path-independence audited at vector level (`F467`, audited by `P470`,
  packaged by `N511`). This does not lift residual sign to a physical datum, does not export a global selector atlas, and does not discharge global `QW-2191`.
- Update (`2026-03-16`): a follow-up audit confirms the oriented transport cocycle holds only on the exported vector section (and hence on transported rays/projectors),
  and does **not** hold as an operator-level matrix identity `O_jk O_ij = O_ik` on the full carrier (`P471`, packaged by `N512`). This blocks any false-pass upgrade
  of the oriented lift into a full operator-level transition groupoid.
- Update (`2026-03-15`): probe-level control-lane transition data is also exported from the **declared residual local-diagonal control pullback**
  value instantiation: the `pair1 -> pair2` cross-block polar orthogonal factor `Q ∈ O(2)` and its induced rotation angle `alpha_cross` (`P460`).
  On the current instantiated matrix, `alpha_cross ≈ 0` while `alpha_12 = π/2`; these are distinct lane-scoped quantities and must not be conflated into a global gluing claim.
- Update (`2026-03-16`): the arithmetic nonvanishing condition behind the `P461` scan is now theorem-level: `N514` proves `F_k(ord_{Z_n})≠0` for all `n≥2` and all `k∈{1,…,n-1}`.
  `P461` remains a computational scope/regression scan only and does not promote any `n≠12` carrier into the physical `QW-2190` scaffold.
- Update (`2026-03-16`): the repo now exports a strict **global selector atlas** and a strict **global selector transition/gluing object** on the declared strict domain `C_v1`,
  discharging the global-scope atlas/transition target `T170` (`F469`, packaged by `N515`). This does **not** imply strict-core selector closure, global `QW-2191` discharge, or ToE closure.
- Update (`2026-03-16`): the repo now exports an explicit strict **global projective selector state object** on the declared strict domain `C_v1`
  (`F470`, packaged by `N516`). This remains projective/span (residual sign is gauge at state level), does **not** export a sign-sensitive directed orientation datum,
  does **not** discharge global `QW-2191`, and does **not** imply strict-core selector closure or ToE closure.
- Update (`2026-03-16`): `P474` audits that the exported global projective selector state is projector-level glued/transported consistently by the exported global
  selector transition operators on `{pair1..pair5}` (ray/projector level only; no sign lift). `N519` packages that residual `Z2` sign can be frozen as gauge for the
  exported global projective selector atlas/transition/state objects, without changing those objects.
- Update (`2026-03-16`): `P475` records the professorial **projective-only continuation** decision: proceed with strict closure tasks using projector/span semantics;
  defer the directed/sign-sensitive selector-state lift (`H37`/`T171`) to an explicit future branch.
- Update (`2026-03-15`): moreover, one cautious follow-up probe exports an explicit **probe-level** mode-index assignment *candidate*
  on `Z_24` induced by the same defect-angle rule (including numeric basis vectors `u_{m,±}` on each `pair_m`) (`P462`).
  This probe-level candidate does not promote `n=24` into the `QW-2190` physical scaffold.
- Update (`2026-03-15`): a minimal typed `Z_24` carrier + regular action primitive is exported as strict infrastructure for cautious scope-extension work
  (`F458`: `I_24_v1`, `Z_24_v1`, `tau_Z24_v1`), without any physical identification with the strict `n=12` scaffolds.
- Update (`2026-03-16`): the repo now exports a strict scope-extension infrastructure upgrade on typed `Z_24`: the Shannon element‑order reference datum `r_ord_z24`
  and a strict `Z_24` mode-index assignment object (`F468`), packaged by a theorem-level `O(2)->Z2` cut on all `pair_m (m=1..11)` on `Z_24` (`N513`).
  This remains non-physical and does not promote `n=24` into the `QW-2190` scaffold.
- Update (`2026-03-15`): the diagonal/local lane exports a strict-derived **full Psi-sector Hessian eigensystem value instantiation**
  `H_psi := K_total + (m0^2 I + D_local_residual)` (numeric 12×12 matrix + eigenvalues/eigenvectors) (`F459`).
  This supports the “light = linearized eigenmodes” reading in declared scope, but remains lane-scoped and does not imply host matching, selector closure, or ToE closure.
- Update (`2026-03-15`): a projection audit exports how the `H_psi` eigenmodes from `F459` decompose in the two exported
  mode-index assignment bases (`F453` diagonal/local and `F454` Shannon) (`P463`). This is a lane-scoped linear-algebra audit and does **not**
  claim that either basis diagonalizes the full `H_psi` operator.
- Update (`2026-03-15`): to avoid any false dependence on residual eigenvector sign conventions in `F459`, the repo exports
  sign‑gauge‑invariant rank‑one spectral projectors `P_j := |v_j><v_j|` for `H_psi` (`F460`) and packages the sign gauge‑irrelevance statement as a strict theorem (`N504`).
- Update (`2026-03-15`): a probe exports a **pair‑plane weight profile** for the `H_psi` eigenprojectors,
  `w_{j,label} := tr(Π_label P_j)` across labels `{e0,pair1..pair5,e6}` (`P464`), and a value‑instantiation theorem packages the conclusion
  that the current `H_psi` eigenmodes are strongly mixed across pair planes (so no implied diagonalization) (`N505`).
- The sigma-int corridor theta-supply and sigma-int → residual bridge theorems introduced in Release 6.4 remain in force
  (e.g. `F451/N489`, `P451`, `F452/N490`, `N491`) and are not re-derived here.

---

## ENGLISH VERSION

## 0) Delta Since Release 6.4

Release 6.5 adds the first strict **value-instantiated diagonal/local uniqueness cuts** on the `n=12` carrier:

1. a strict-derived per-site provider chain exists for the diagonal/local accelerator lane (`F447/N483`, audited by
   `P448`), making the previously “harness-only” diagonal computations actually computable without hidden slots,
2. the canonical diagonal/local residual mode‑2 defect on `pair1` is now decided as **nonzero** on the exported strict
   instantiation: $|F_2(d)| \approx 12.88048 \neq 0$ (`N482`), hence the diagonal/local sector cuts `O(2)` on `pair1`
   down to residual `Z2`,
3. the same instantiated diagonal/local residual profile cuts the continuous `O(2)` family on **all** Fourier-degenerate
   pairs `pair_m (m=1..5)` on `n=12` (`N485`),
4. therefore the `QW-2191` continuous `O(2)` family obstruction is **scoped-discharged** on the exported canonical
   diagonal/local lane: `O(2) -> Z2` on every `pair_m` (`N487`), while kernel-alone/global scope remains open,
5. an explicit strict-derived diagonal/local mode-index assignment basis object is exported, covering the full Fourier
   scaffold `(e0, pair1..pair5, e6)` (`F453`),
6. the diagonal/local lane exports a theorem-level “internal orientation datum exists” statement in the lane-scoped,
   axis-only sense (residual sign remains) (`N492`),
7. residual sign flips are gauge-equivalent (pure conjugations) for the `QW-2190` embedding audits (`P452/N493`), hence
   the diagonal/local `QW-2190` embedding is unique **up to conjugation** in that scope (`N494`),
8. full continuous `O(2)` rotations inside the embedded pair planes are gauge-equivalent (conjugation-only) for the
   same `QW-2190` embedding audits (`N495`, audited by `P454`),
9. independently, the strict Shannon element‑order reference lane cuts `O(2)` down to residual `Z2` on **all**
   Fourier-degenerate pairs `pair_m (m=1..5)` on `n=12` and exports a strict-core mode-index assignment basis object
   `ModeIndexAssignment_shannon_element_order_reference_strict_core_v1` (`N496`, executed by `F454`, packaged by `N500`),
10. in the declared `R1` target-slot semantics (`span{u_1,u_2}`), residual `Z2` sign flips of the representative vectors are gauge-irrelevant
    for the target slot as a span object (`N501`).
11. the diagonal/local and Shannon mode-index assignments are aligned on all pairs up to residual `Z2` sign (audit `P455`),
12. repo hygiene probes/dashboards were updated accordingly (`P432`, `P441`, `P442`), without any promotion to global
    closure.
13. sigma‑int slot‑free theta‑pair axes (pair1/pair2) align with both exported mode‑index assignment bases up to residual `Z2`
    sign (audit `P456`, packaged `N498`).
14. the declared residual local‑diagonal control pullback `M_control_residual_diag_declared` is explicitly restricted to the
    sigma‑int orientation slice basis `(u_1,u_2)` (`P457`) and value‑instantiated under the current strict‑derived provider via the
    conditional `N477` rewrite (`P458`), exporting a concrete 2×2 matrix while keeping host cancellation unclaimed.
15. the full 4×4 declared residual local‑diagonal control pullback on the control basis `(c1,s1,c2,s2)` is also value‑instantiated under
    the same conditional `N477` rewrite and strict‑derived provider (`P459`), exposing `pair1`/`pair2` blocks and cross‑block couplings
    without any host‑matching claim.
16. the strict sigma‑int → residual target‑slot bridge/export‑map object is explicitly upgraded to carry the slot‑free theta‑pair outputs and
    the corresponding strict‑core `R1` inhabitant instance as map outputs (`F455`, packaged by `N499`), without theta inputs and without implied selector closure
    (the `v1` export-map object remains sign-only).
17. the repo exports a minimal strict-core downstream operator on `V_1=span{c1,s1}` derived only from the materialized sigma‑int `u_1` direction:
    `A_1(pair1) := |u_1><u_1|` (projector; residual `Z2` sign gauge invariant). This closes the `P2` operator-stage target in declared scope (`F456`),
    but does **not** identify the operator with the extension-only `A_1_ext` of the `H/O` lane, does not imply selector closure, and does not discharge global `QW-2191`.
18. residual `Z2` sign can be frozen as a tracked gauge/convention layer for the currently exported downstream objects where sign is
    provably gauge-irrelevant (packaged by `N502`), without promoting any sign-sensitive physical orientation datum or any global discharge.
19. derived transition data `alpha_12 := (theta_2 - theta_1) mod 2π` for `pair1/pair2` is exported from the strict sigma-int slot-free theta-pair
    supply (`F457`), without implying a global selector transition/gluing object.
20. probe-level control-lane transition data is also exported from the declared residual control pullback value instantiation (`P459`): the cross-block polar orthogonal factor
    `Q ∈ O(2)` and its induced rotation angle `alpha_cross` (`P460`). On the current instantiated matrix, `alpha_cross ≈ 0` while `alpha_12 = π/2`; they are distinct lane-scoped quantities.
21. an explicit **lane-scoped** `pair1↔pair2` chart-transport operator `O_12` on the declared `n=12` Fourier carrier is exported, derived only from the strict sigma-int
    transition angle `alpha_12` (`F461`), with projector-level sign-gauge-irrelevance packaged as a strict theorem (`N506`) and audited by a dedicated probe (`P465`).
    This does not export a global selector atlas/gluing object and does not discharge global `QW-2191`.
22. building on `O_12`, the repo exports a concrete **two-chart projector operator section** on `{pair1,pair2}`:
    `A_2(pair2) = O_12 A_1(pair1) O_12^T` (`F462`), packaged as a strict theorem (`N507`) and audited by `P466`. This remains lane-scoped and does not upgrade
    to a global selector atlas / gluing structure and does not discharge global `QW-2191`.
23. the repo exports an explicit lane-scoped **two-chart selector atlas stub** on `{pair1,pair2}` with an overlap-domain declaration, transition data (`O_12`),
    and projector-level gluing data (`A_2 = O_12 A_1 O_12^T`) (`F463`). This is not a global selector atlas and does not discharge global `QW-2191`.
24. the repo exports an explicit lane-scoped **three-chart** selector-atlas ingredient on `{pair1,pair2,pair3}` at projector level, including explicit
    cocycle/path-independence data on the glued projector section (`F464`), audited by `P467` and packaged as a strict theorem (`N508`). This remains lane-scoped
    and does not export a global selector atlas nor discharge global `QW-2191`.
25. the repo exports an explicit lane-scoped **five-chart** selector-atlas ingredient on `{pair1..pair5}` at projector level, including explicit **local**
    cocycle/path-independence data on adjacent triple overlaps (`F465`), audited by `P468` and packaged as a strict theorem (`N509`). This remains lane-scoped
    and does not export a global selector atlas nor discharge global `QW-2191`.
26. the repo exports additional lane-scoped axis-only long-edge chart-transport operators (`O_14`, `O_15`, `O_25`) and upgrades the five-chart `{pair1..pair5}`
    projector-level selector-atlas ingredient to explicit cocycle/path-independence audit data for **all** triple overlaps (`F466`), audited by `P469` and packaged
    as a strict theorem (`N510`). This remains lane-scoped and does not export a global selector atlas nor discharge global `QW-2191`.
27. the repo exports a lane-scoped **oriented transport (α mod 2π) lift** at vector level on `{pair1..pair5}`, induced by the currently exported representative
    vectors `u_1..u_5`, explicitly tracked as a gauge/convention layer (sign-tracked) (`F467`). A dedicated probe audits vector transport and full triple
    cocycle/path-independence on the exported instance (`P470`), and a theorem packages the statement as a convention-layer result (`N511`). This remains lane-scoped
    and does not export a global selector atlas nor discharge global `QW-2191`.
28. a strict hygiene follow-up audits that the oriented transport cocycle holds only on the exported glued vector section (and hence on transported rays/projectors),
    and does **not** hold as an operator-level matrix identity `O_jk O_ij = O_ik` on the full carrier (`P471`, packaged by `N512`). Therefore the oriented lift must not be
    treated as an operator-level transition groupoid transporting arbitrary vectors.
29. the repo exports a strict scope-extension infrastructure upgrade on typed `Z_24`: a Shannon element‑order reference datum `r_ord_z24` and a strict `Z_24`
    mode-index assignment basis object (`F468`), packaged by a theorem-level `O(2)->Z2` cut on all `pair_m (m=1..11)` on `Z_24` (`N513`), without any physical promotion
    of `n=24` into the `QW-2190` scaffold.

## 1) One-Page Strict Status (6.5)

### 1.1 What is now structurally sharp

1. `QW-2191` remains true in its kernel-alone/global scope: translation-invariant kernel data does not canonically fix
   axes inside each degenerate Fourier pair plane without an additional strict internal source.
2. The canonical diagonal/local lane is now **numerically instantiated in strict-derived form** and decides the
   diagonal `O(2)` cuts on `n=12`:
   - `pair1` decided by `T166` (`N482`),
   - all `pair_m (m=1..5)` decided (`N485`).
3. On that lane, the continuous `O(2)` family is cut down to residual `Z2` on every degenerate pair plane (`N487`).
4. The repo exports an explicit strict-derived mode-index assignment basis object realizing that axis selection
   (up to sign) (`F453`), and packages it as a lane-scoped internal orientation datum (`N492`).
5. For the `QW-2190` SU(3)/SU(2) embedding audits, both residual sign flips and the full continuous `O(2)` basis-rotation
   freedom are gauge: invariance/closure audits are unchanged (conjugation equivalence) (`N493`, `N495`, audited by
   `P452`, `P454`), so the embedding audits are well-defined up to conjugation in the declared scope (`N494`).
6. Independently, on the strict Shannon element‑order reference lane, the cross‑entropy objective cuts `O(2)` down to
   residual `Z2` on all `pair_m (m=1..5)` on `n=12` (`N480`, `N488`, `N496`), enabling export of a full strict-core
   mode-index assignment basis object without per-site diagonal/local providers (`F454`).
7. The diagonal/local and Shannon strict mode-index assignments agree on the selected axes (up to residual sign) on all
   `pair_m (m=1..5)` on `n=12` (audit `P455`), so the axis choice is not an artifact of a single lane.

### 1.2 What is still missing (no false pass)

The repo still does **not** export:

1. strict-core selector closure / admissible `S_sel_int`,
2. a strict sign-sensitive / directed selector state datum or observable distinguishing `u` from `-u` (audits `H36`/`H37`), beyond the projective/span level,
3. an axiom-free **global** discharge of `QW-2191` beyond the declared diagonal/local lane and `n=12` scope,
4. a strict physical sign/orientation convention if any downstream claim depends on absolute sign (unless separately
   proven gauge-irrelevant on that lane),
5. ToE closure.

## 2) The diagonal/local `O(2) -> Z2` cut mechanism (n=12)

Let `D = diag(d_0,...,d_{11})` be a diagonal operator in the site basis on the `n=12` carrier.
For each degenerate Fourier pair plane `pair_m = span{c_m,s_m}` (`m=1..5`), define the diagonal defect:

$$
F_{2m}(d) := \sum_{k=0}^{11} d_k\, e^{i\frac{4\pi m k}{12}} \in \mathbb{C}.
$$

The induced diagonalization angle is:

$$
\theta_*(m) := \frac12\,\mathrm{atan2}(\mathrm{Im}\,F_{2m},\ \mathrm{Re}\,F_{2m})
\quad (\mathrm{mod}\ \pi).
$$

On the current exported strict-derived canonical diagonal/local residual profile `d_k` (computed by `P437` on strict-derived
inputs provided by `F447/N483`), the evaluated defect magnitudes are (`N485` / `P449`):

```text
|F2(d)|  ≈ 12.88048321986275
|F4(d)|  ≈ 13.37099987006163
|F6(d)|  ≈ 25.57352236877763
|F8(d)|  ≈ 13.37099987006158
|F10(d)| ≈ 12.88048321986272
```

Hence `|F_{2m}(d)| ≠ 0` for all `m=1..5`, cutting each `pair_m` continuous `O(2)` down to residual `Z2` on the diagonal/local lane.

## 3) Exported mode-index assignment object (F453)

The exported strict-derived object:

```text
ModeIndexAssignment_canonical_local_diagonal_strict_derived_v1
```

provides an explicit orthonormal basis:

- `(e0, e6)` for the nondegenerate modes, and
- for each `pair_m`, the rotated basis `(u_{m,+}, u_{m,-}) = (c_m,s_m) R(\theta_*(m))`,

with the stated hard limit: the basis is canonical only **up to residual sign** on each `pair_m` (no sign-sensitive physical convention implied).

## 4) Residual sign is gauge for `QW-2190` embedding audits (N493/N494)

For the `QW-2190` SU(3)/SU(2) embedding audits, flipping any subset of the residual signs on the exported basis vectors is:

1. projector-invariant (subspace invariance audits unchanged), and
2. generator-conjugate (Lie-closure audits unchanged),

so the diagonal/local lane canonicalizes the `QW-2190` embedding uniquely up to conjugation in the declared scope (`N494`).

## 5) Next honest strict moves (as of 2026-03-16)

Dashboards continue to mark `H37` (sign-distinction state: a sign-sensitive / directed selector state datum or observable distinguishing `u` from `-u`) as the next strict frontier for a **directed** continuation under `QW-2191` discipline (`P438`, `P441`), after the export of the global projective selector state object on `C_v1` (`F470/N516`) and the discharge of `T170` (`F469/N515`).

However, the continuation is now explicitly bifurcated (`T171`), and `P475` selects the **projective-only** continuation:

1. treat the exported global selector state as the strict physical state object at the ray/projector level for the declared closure stack, keeping residual sign as a gauge/convention layer where proven irrelevant (`N502`, `N519`),
2. proceed with strict-only ToE closure tasks that do not require a sign-sensitive orientation datum (projective-only compatible); the next concrete bottleneck is the kernel-split-robust canonical-ontology-supported direct formal `c1s1` family route, currently tracked at `P623`.

Update: a coefficient-filled **declared** control pullback `M_control = T_control^T H_PsiPsi T_control` is now exported from the coefficient-filled canonical Psi block on the declared transport support (`P476`, relying on `R11/R12`), but no host-to-canonical Psi-block matching witness is exported, so no strict existing-feedback promotion is claimed.

Update: a strict obstruction is now explicit: any direction-free `Aut(Z_12)`-invariant reference weight family (including the ord-reference weights `ord_Z12`, `r_ord`) cannot distinguish sign on the current exported `pair1` sine axis via a scalar of the form `Σ_x w(x) u_1(x)` (`F472`, packaged by `N518`; strengthening the ord-specific obstruction `F471/N517`). Therefore `H37` requires an explicit reflection-breaking/orientation source or a different observable class.

Probe hygiene: `P472` scans exported `generated/*.json` artifacts and finds no strict(-derived) **weight-like** reflection-breaking per-site array outside non-canonical marked-site `K_total` row vectors; therefore no strict sign-sensitive physical orientation datum is presently exported by repackaging existing direction-free references.

Update: as the next narrow `B3` continuation step, the repo now exports an explicit lane-scoped **chart-transport operator** `O_12` between `pair1↔pair2`
on `n=12`, derived only from the strict sigma-int `alpha_12` transition (`F461`). This is a sign-free projector-level gluing ingredient
(packaged as a strict theorem `N506` and audited by `P465`); global selector atlas/transition obligations are now discharged at object level (`F469/N515`), and a global projective selector state object is exported (`F470/N516`).

Update: the lane-scoped `O_12` transport is now used to export an explicit **two-chart glued projector operator section**
`A_2(pair2) = O_12 A_1(pair1) O_12^T` (`F462`), packaged as a strict theorem (`N507`) and audited by `P466`. This remains lane-scoped and does not imply strict-core selector closure nor global `QW-2191` discharge.

Update: the repo exports an explicit lane-scoped **two-chart selector atlas stub** on `{pair1,pair2}` with an overlap-domain declaration, transition data (`O_12`),
and gluing data (`A_2 = O_12 A_1 O_12^T`) (`F463`). This reduces the “no overlap-domain declaration” gap only in that declared scope; global atlas/transition objects on `C_v1` are exported separately (`F469/N515`) at projector-section level (no operator-level groupoid promotion).

Update: the repo exports an explicit lane-scoped **three-chart** selector-atlas ingredient on `{pair1,pair2,pair3}` at projector level, including explicit
cocycle/path-independence data on the glued projector section (`F464`), audited by `P467` and packaged as a strict theorem (`N508`). This remains lane-scoped
and does not export a global selector atlas nor discharge global `QW-2191`.

Update: continuing the same minimal strategy, the repo exports an explicit lane-scoped **five-chart** selector-atlas ingredient on `{pair1..pair5}` at projector level,
including explicit **local** cocycle/path-independence data on adjacent triple overlaps (`F465`), audited by `P468` and packaged as a strict theorem (`N509`). This remains lane-scoped
and does not export a global selector atlas nor discharge global `QW-2191`.

Update: continuing the same minimal `B3` strategy, the repo exports additional lane-scoped axis-only long-edge chart-transport operators (`O_14`, `O_15`, `O_25`) and upgrades
the five-chart `{pair1..pair5}` projector-level atlas ingredient to explicit cocycle/path-independence audit data for **all** triple overlaps (`F466`), audited by `P469` and packaged
as a strict theorem (`N510`). This remains lane-scoped and does not export a global selector atlas nor discharge global `QW-2191`.

Update: continuing the same minimal `B3` strategy, the repo exports a lane-scoped **oriented transport (α mod 2π) lift** at vector level on `{pair1..pair5}`, induced by the
currently exported representative vectors `u_1..u_5`, explicitly tracked as a gauge/convention layer (sign-tracked) (`F467`). A dedicated probe audits vector transport and full
triple cocycle/path-independence on the exported instance (`P470`), and a theorem packages the statement as a convention-layer result (`N511`). This does not lift residual sign to a
physical datum, does not export a global selector atlas, and does not discharge global `QW-2191`.

Update: a strict hygiene follow-up audits that the oriented transport cocycle holds only on the exported glued vector section (and hence on transported rays/projectors),
and does **not** hold as an operator-level matrix identity `O_jk O_ij = O_ik` on the full carrier (`P471`, packaged by `N512`). Therefore the oriented lift must not be treated as
an operator-level transition groupoid transporting arbitrary vectors.

Update: the strict Shannon element‑order reference lane now provides one explicit internal symmetry-breaking ingredient
on the full `n=12` Fourier scaffold (via `F454`), but this is still below strict-core selector closure and does not
constitute a global `QW-2191` discharge.

Update: the arithmetic nonvanishing condition behind the `P461` scan is now theorem-level: `N514` proves `F_k(ord_{Z_n})≠0` for all `n≥2` and all `k∈{1,…,n-1}`.
`P461` remains a computational scope/regression check only and does not promote any `n≠12` carrier into the physical `QW-2190` scaffold.

Update: one cautious follow-up exports a probe-level `Z_24` mode-index assignment *candidate* (numeric basis vectors)
induced by the same defect-angle rule (`P462`). This does not promote `Z_24` into the strict physical mode scaffold.

Update: the repo now exports a strict scope-extension infrastructure upgrade on typed `Z_24`: the Shannon element‑order reference datum `r_ord_z24`
and a strict `Z_24` mode-index assignment object (`F468`), packaged by a theorem-level `O(2)->Z2` cut on all `pair_m (m=1..11)` on `Z_24` (`N513`).
This remains non-physical and does not promote `n=24` into the `QW-2190` scaffold.

---

## WERSJA POLSKA

## 0) Delta względem Release 6.4

Release 6.5 dodaje pierwsze ścisłe (strict) **zainstancjonowane wartościowo** domknięcia na pasie diagonal/local dla nośnika `n=12`:

1. istnieje ścisły łańcuch dostawcy wartości per‑site dla pasa diagonal/local (`F447/N483`, audyt `P448`), co usuwa wcześniejszą
   “nieobliczalność” wynikającą z braku wartości wejściowych,
2. `T166` (decyzja czy $F_2(d)=0$) jest teraz rozstrzygnięty jako **niezerowy** na wyeksportowanej instancji strict-derived:
   $|F_2(d)| \approx 12.88048 \neq 0$ (`N482`), więc sektor diagonal/local tnie `O(2)` na `pair1` do residual `Z2`,
3. ten sam profil diagonal/local tnie ciągłą rodzinę `O(2)` na **wszystkich** parach `pair_m (m=1..5)` (`N485`),
4. w konsekwencji przeszkoda `QW-2191` (ciągła rodzina `O(2)`) jest **rozładowana w zakresie** wyeksportowanego pasa diagonal/local:
   `O(2) -> Z2` na każdej parze (`N487`), ale globalny zakres (kernel-alone) pozostaje otwarty,
5. wyeksportowano jawny obiekt bazowy przypisania indeksów modów (mode-index assignment) dla pełnego scaffoldu Fouriera (`F453`),
6. pas diagonal/local eksportuje twierdzenie “istnieje wewnętrzny datum orientacji” w sensie osiowym (axis-only; residual znak pozostaje) (`N492`),
7. residualne flippy znaku `Z2` są “gauge” dla audytów embeddingu `QW-2190` (czysta koniugacja) (`P452/N493`), więc embedding jest
   kanoniczny **z dokładnością do koniugacji** w zadeklarowanym zakresie (`N494`),
8. ponadto pełna ciągła swoboda rotacji `O(2)` w parach jest “gauge” dla tych samych audytów embeddingu `QW-2190`
   (`N495`, audyt `P454`),
9. niezależnie, pas strict Shannon “element‑order reference” tnie `O(2)` do residual `Z2` na **wszystkich** parach
   `pair_m (m=1..5)` na `n=12` i eksportuje strict-core obiekt przypisania osi
   `ModeIndexAssignment_shannon_element_order_reference_strict_core_v1` (`N496`, wykonane przez `F454`, opakowane przez `N500`),
10. w zadeklarowanej semantyce target‑slot `R1` (`span{u_1,u_2}`) residualne flippy znaku `Z2` wektorów reprezentantów są
    gauge‑irrelewant dla target‑slot jako obiektu typu “span” (`N501`).
11. diagonal/local i Shannon mode-index assignment są zgodne na wszystkich parach co do wyboru osi (z dokładnością do residualnego znaku) (audyt `P455`),
12. sondy/dashboards higieniczne zostały zaktualizowane (`P432`, `P441`, `P442`) bez promocji do globalnego domknięcia.
13. osie theta‑pair na korytarzu sigma‑int (pair1/pair2) są zgodne z oboma wyeksportowanymi mode-index assignment (diagonal/local `F453` i Shannon `F454`)
    z dokładnością do residualnego znaku `Z2` (audyt `P456`, opakowane przez `N498`).
14. zadeklarowany pullback kontrolny residualnego sektora diagonal/local
    `M_control_residual_diag_declared = T_control^T D_local_residual T_control` został jawnie ograniczony do “sigma‑int orientation slice”
    w bazie `(u_1,u_2)` (`P457`) i **warunkowo** zainstancjonowany liczbowo na aktualnym strict‑derived providerze przez przepisanie `N477` (`P458`),
    eksportując konkretną macierz 2×2 bez roszczeń o anulację host‑matching.
15. ponadto, pełna macierz 4×4 zadeklarowanego pullbacku residualnego sektora diagonal/local na bazie kontrolnej `(c1,s1,c2,s2)` została
    **warunkowo** zainstancjonowana liczbowo przez to samo przepisanie `N477` i ten sam strict‑derived provider (`P459`), ujawniając bloki
    `pair1`/`pair2` oraz sprzężenia cross‑block bez roszczeń o host‑matching.
16. ponadto, ścisły obiekt bridge/export‑map na korytarzu sigma‑int → residual target‑slot został jawnie ulepszony do
    `Upsilon_residual_datum_sigma_int_bridge_export_map_object_v2`, dołączając wyjścia slot‑free theta‑pair oraz odpowiadającego inhabitanta
    `R1` jako outputy mapy (`F455`, opakowane przez `N499`), bez wejściowych theta i bez implikacji domknięcia selektora (`v1` pozostaje sign‑only).
17. dodatkowo, repo eksportuje minimalny ścisły operator downstream na `V_1=span{c1,s1}` z już zmaterializowanego kierunku sigma‑int `u_1`:
    `A_1(pair1) := |u_1><u_1|` (projector; invariantny na residualny flip znaku `Z2`). To domyka etap operatorowy `P2` w zadeklarowanym scope (`F456`),
    ale **nie** identyfikuje tego operatora z extension-only `A_1_ext` z pasa `H/O`, **nie** implikuje selector closure i **nie** rozładowuje globalnie `QW-2191`.
18. residualny znak `Z2` może zostać jawnie **zamrożony jako warstwa gauge/konwencji** dla aktualnie wyeksportowanych obiektów downstream,
    dla których wykazano gauge‑irrelewantność znaku (opakowane przez `N502`, w tym `alpha_12 mod π`), bez promocji do sign-sensitive “fizycznej orientacji”.
19. pochodna dana przejścia `alpha_12 := (theta_2 - theta_1) mod 2π` dla `pair1/pair2` jest wyeksportowana z strict slot‑free theta‑pair na korytarzu sigma‑int (`F457`),
    bez implikowania globalnego selector transition/gluing object.
20. dodatkowo, na poziomie sondy/control‑lane wyeksportowano dane przejścia z warunkowo zainstancjonowanego zadeklarowanego pullbacku residualnego (`P459`): czynnik ortogonalny `Q ∈ O(2)`
    z dekompozycji polarnej cross‑blocku `pair1 -> pair2` oraz odpowiadający mu kąt rotacji `alpha_cross` (`P460`). Na aktualnej instancji liczbowej `alpha_cross ≈ 0`, podczas gdy `alpha_12 = π/2`;
    są to różne wielkości lane‑scoped i nie wolno ich mieszać w claim o globalnym gluing.
21. theorem-level arytmetyka domyka kryterium defektu element‑order reference: `N514` dowodzi, że dla każdego `n≥2` i każdego `k∈{1,…,n-1}` zachodzi `F_k(ord_{Z_n})≠0`,
    więc defekty `F_{2m}(ord_{Z_n})` są niezerowe na wszystkich parach Fouriera na dowolnym `Z_n` (lane-scoped; bez promocji fizycznej).
    Sonda `P461` pozostaje tylko sprawdzeniem obliczeniowym/regresyjnym (zgodność na `n ∈ {6,8,10,12,14,16,18,20,24}`) i nie promuje żadnego `n≠12` do fizycznego scaffoldu `QW-2190`.
22. ponadto, jako ostrożny follow-up, wyeksportowano jawny **probe-level** kandydat mode-index assignment na `Z_24`
    (numeryczne wektory bazowe `u_{m,±}` na każdej parze Fouriera) indukowany przez tę samą regułę kąta defektu (`P462`).
    To nadal nie jest theorem-level scope extension i nie promuje `n=24` do fizycznego scaffoldu `QW-2190`.
23. ponadto, wyeksportowano minimalną typed infrastrukturę `Z_24` (nośnik + działanie regularne) jako podparcie dla ostrożnych prac scope-extension
    (`F458`: `I_24_v1`, `Z_24_v1`, `tau_Z24_v1`), bez jakiejkolwiek identyfikacji fizycznej z nośnikiem strict `n=12`.
24. dodatkowo, pas diagonal/local eksportuje strict-derived **pełny eigensystem** liczbowej instancji Hessianu sektora `Psi`:
    `H_psi := K_total + (m0^2 I + D_local_residual)` (macierz 12×12 + wartości własne + wektory własne) (`F459`).
    To wspiera interpretację “światło = eigenmody liniaryzacji” w zadeklarowanym zakresie, ale pozostaje lane-scoped i nie implikuje host matching, selector closure ani ToE closure.
25. dodatkowo, wyeksportowano audyt projekcyjny pokazujący jak eigenmody `H_psi` z `F459` rozkładają się w dwóch wyeksportowanych bazach
    mode-index assignment (`F453` diagonal/local oraz `F454` Shannon) (`P463`). To jest lane-scoped audyt algebry liniowej i **nie**
    twierdzi, że którakolwiek z tych baz diagonalizuje pełny operator `H_psi`.
26. aby uniknąć fałszywej zależności od residualnych konwencji znaku wektorów własnych w `F459`, repo eksportuje
    sign‑gauge‑invariant rank‑one projektory spektralne `P_j := |v_j><v_j|` dla `H_psi` (`F460`) oraz pakuje twierdzenie o gauge‑irrelewantności znaku jako strict theorem (`N504`).
27. dodatkowo, sonda eksportuje **profil wag pair‑plane** dla eigenprojektorów `H_psi`:
    `w_{j,label} := tr(Π_label P_j)` po labelach `{e0,pair1..pair5,e6}` (`P464`), a value‑instantiation theorem pakuje wniosek,
    że na aktualnej instancji `H_psi` eigenmody są silnie zmieszane pomiędzy pair planes (więc nie ma prawa sugerować diagonalizacji) (`N505`).
28. ponadto, repo eksportuje jawny **lane-scoped** obiekt transportu/klejenia chartów selektora pomiędzy `pair1` i `pair2` na nośniku `n=12`:
    operator ortogonalny `O_12` zbudowany wyłącznie z `alpha_12` na korytarzu sigma‑int (`F461`). Pakuje się przy tym ścisłe twierdzenie,
    że transport na poziomie projektorów `P(u) = |u><u|` jest gauge‑irrelewant względem residualnego flipa znaku (`N506`), a sonda `P465`
    audytuje ortogonalność, transport płaszczyzn oraz próbkowany transport `u_{m,θ}`/projektorów. Nie jest to globalny selector atlas ani globalny gluing object
    i nie rozładowuje globalnie `QW-2191`.
29. ponadto, na bazie `O_12` repo eksportuje konkretną **dwu‑chartową sekcję operatorową** na `{pair1,pair2}` na poziomie projektorowym (sign‑free):
    `A_2(pair2) = O_12 A_1(pair1) O_12^T` (`F462`), opakowane jako twierdzenie strict (`N507`) i audytowane sondą (`P466`). Jest to nadal lane‑scoped i nie
    awansuje do globalnego selector atlas/gluing object ani do globalnego rozładowania `QW-2191`.
30. ponadto, repo eksportuje jawny lane-scoped **dwu‑chartowy atlas selektora** na `{pair1,pair2}` z deklaracją overlap-domain, danymi przejścia (`O_12`)
    oraz danymi klejenia sekcji operatorowej (`A_2 = O_12 A_1 O_12^T`) (`F463`). To nie jest globalny atlas selektora i nie rozładowuje globalnie `QW-2191`.
31. ponadto, repo eksportuje jawny lane-scoped **trzy‑chartowy** składnik atlasu selektora na `{pair1,pair2,pair3}` na poziomie projektorowym (sign‑free),
    wraz z jawnymi danymi cocycle/path‑independence na poziomie sklejonej sekcji operatorowej (`F464`), audytowane sondą (`P467`) i opakowane jako twierdzenie strict (`N508`).
    Jest to nadal lane‑scoped (overlap rozumiany jako overlap artefaktów), nie jest globalnym atlasem selektora i nie rozładowuje globalnie `QW-2191`.
32. ponadto, repo eksportuje jawny lane-scoped **pięcio‑chartowy** składnik atlasu selektora na `{pair1..pair5}` na poziomie projektorowym (sign‑free),
    wraz z jawnymi **lokalnymi** danymi cocycle/path‑independence dla sąsiednich trójek (1‑2‑3, 2‑3‑4, 3‑4‑5) na poziomie sklejonej sekcji projektorowej (`F465`),
    audytowane sondą (`P468`) i opakowane jako twierdzenie strict (`N509`). Jest to nadal lane‑scoped, nie jest globalnym atlasem selektora i nie rozładowuje globalnie `QW-2191`.
33. ponadto, repo eksportuje dodatkowe operatory transportu chartów w wariancie **axis-only** na brakujących “długich” krawędziach (`O_14`, `O_15`, `O_25`)
    oraz upgrade’uje pięcio‑chartowy składnik `{pair1..pair5}` do jawnych danych cocycle/path‑independence dla **wszystkich** trójek na poziomie sklejonej sekcji projektorowej (`F466`),
    audytowane sondą (`P469`) i opakowane jako twierdzenie strict (`N510`). Jest to nadal lane‑scoped, nie jest globalnym atlasem selektora i nie rozładowuje globalnie `QW-2191`.
34. ponadto, repo eksportuje jawny lane-scoped lift transportu chartów do **oriented** `α mod 2π` na poziomie wektorów reprezentantów na `{pair1..pair5}`,
    indukowany przez aktualnie wyeksportowane `u_1..u_5` i jawnie śledzony jako warstwa gauge/konwencji (sign‑tracked) (`F467`). Sonda `P470` audytuje transport
    wektorów i pełne relacje trójkowe cocycle/path‑independence na wyeksportowanej instancji, a `N511` pakuje to jako twierdzenie strict w dyscyplinie “convention layer”.
    Nie jest to promocja do fizycznego datumu znaku, nie jest to globalny atlas selektora i nie rozładowuje globalnie `QW-2191`.
35. dodatkowo, jako ścisła higiena, sonda pokazuje, że relacje cocycle/path‑independence dla transportu “oriented” zachodzą tylko na wyeksportowanej sekcji
    wektorowej (a więc na transportowanych rayach/projektorach), natomiast **nie** zachodzą jako operatorowa równość macierzowa `O_jk O_ij = O_ik` na pełnym nośniku
    (`P471`, opakowane przez `N512`). W szczególności nie wolno traktować tego liftu jako operatorowego transition groupoid transportującego dowolne wektory.
36. ponadto, repo eksportuje już ścisły upgrade infrastrukturalny scope-extension na typed `Z_24` dla pasa Shannon element‑order reference:
    datum `r_ord_z24` oraz jawny obiekt mode-index assignment na `Z_24` (`F468`), opakowane twierdzeniem theorem-level o cięciu `O(2)->Z2` na wszystkich parach
    `pair_m (m=1..11)` na `Z_24` (`N513`). To nadal jest **nie‑fizyczne** i nie promuje `n=24` do scaffoldu `QW-2190`.

## 1) Jednostronicowy status strict (6.5)

### 1.1 Co jest teraz ostre strukturalnie

1. `QW-2191` pozostaje prawdziwe w zakresie kernel-alone/global: sam kernel translacyjnie niezmienniczy nie wybiera osi w parach zdegenerowanych
   bez dodatkowego wewnętrznego źródła strict.
2. Pas diagonal/local jest teraz **zainstancjonowany liczbowo** w reżimie strict-derived i rozstrzyga cięcia `O(2)` na `n=12`:
   `pair1` (`N482`) oraz wszystkie `pair_m (m=1..5)` (`N485`).
3. Na tym pasie ciągła rodzina `O(2)` jest zredukowana do residual `Z2` na każdej parze (`N487`).
4. Repo eksportuje jawny obiekt bazy przypisania modów realizujący ten wybór osi (z dokładnością do znaku) (`F453`) i pakuje to jako datum orientacji (`N492`).
5. Dla audytów embeddingu `QW-2190` nie tylko residualny znak, ale też pełna swoboda rotacji `O(2)` jest “gauge”
   (koniugacja; `N493`, `N495`, audyty `P452`, `P454`), więc embedding jest kanoniczny do koniugacji (`N494`).
6. Niezależnie, na pasie strict Shannon “element‑order reference”, obiektyw cross‑entropy tnie `O(2)` do residual `Z2`
   na wszystkich `pair_m (m=1..5)` na `n=12` (`N480`, `N488`, `N496`), umożliwiając eksport pełnej bazy przypisania osi
   bez per‑site dostawców diagonal/local (`F454`).
7. Diagonal/local oraz Shannon przypisania osi są zgodne co do wyboru osi na wszystkich `pair_m (m=1..5)` na `n=12`
   (audyt `P455`), więc wybór osi nie jest artefaktem pojedynczej lane.

### 1.2 Czego nadal brakuje (bez false pass)

Repo nadal **nie** eksportuje:

1. strict-core selector closure / dopuszczalnego `S_sel_int`,
2. ścisłego znako‑czułego / kierunkowego datumu stanu selektora lub obserwabli rozróżniającej `u` od `-u` (audyty `H36`/`H37`), ponad poziom projektowy/span,
3. aksjomatycznie wolnego **globalnego** rozładowania `QW-2191` poza pasem diagonal/local i zakresem `n=12`,
4. ścisłej konwencji znaku/orientacji tam, gdzie downstream wymaga absolutnego znaku (o ile nie jest osobno pokazane, że znak jest gauge),
5. domknięcia ToE.

## 2) Następny uczciwy ruch (stan: 2026-03-16)

Dashboards nadal wskazują `H37` (rozróżnienie znaku: znako‑czuły / kierunkowy datum stanu selektora lub obserwabla rozróżniająca `u` od `-u`) jako następny strict frontier **dla gałęzi kierunkowej** w dyscyplinie `QW-2191` (`P438`, `P441`), po eksporcie globalnego projektowego stanu selektora na `C_v1` (`F470/N516`) i po rozładowaniu `T170` (`F469/N515`).

Jednak kontynuacja jest teraz jawnie rozdzielona (`T171`), a `P475` wybiera gałąź **projekcyjną**:

1. traktować globalny stan selektora jako obiekt promienia/projektora w zadeklarowanym stacku domknięcia, zamrażając residualny znak jako gauge/konwencję tam, gdzie pokazano gauge‑irrelewantność (`N502`, `N519`),
2. przejść do strict-only zadań domknięcia ToE, które nie wymagają znako‑czułego datumu orientacji (kompatybilne z gałęzią projekcyjną); najbliższy konkretny bottleneck to kernel-split-robust canonical-ontology-supported direct formal route rodziny `c1s1`, aktualnie śledzony na froncie `P623`.

Update: repo eksportuje teraz coefficient-filled **declared** control pullback `M_control = T_control^T H_PsiPsi T_control` z coefficient-filled canonical `H_PsiPsi` na zadeklarowanym transporcie (`P476`, oparte o `R11/R12`), ale bez witnessu dopasowania hosta do nośnika kanonicznego (`C10_B1`), więc bez promocji do strict existing-feedback.

Update: jawne jest teraz jedno ścisłe ograniczenie: żadna direction-free (tj. `Aut(Z_12)`-invariant) rodzina wag referencyjnych (w tym wagi ord-reference `ord_Z12`, `r_ord`) nie potrafi rozróżnić znaku na aktualnie wyeksportowanej osi `pair1` typu sine poprzez skalar postaci `Σ_x w(x) u_1(x)` (`F472`, opakowane przez `N518`; wzmacniając ord-specyficzną przeszkodę `F471/N517`). W konsekwencji `H37` wymaga jawnego źródła łamiącego symetrię odbicia/orientacji albo innej klasy obserwabli.

Higiena sondy: `P472` skanuje wyeksportowane artefakty `generated/*.json` i nie znajduje żadnej ścisłej (strict lub strict-derived) **weight-like** per-site tablicy łamiącej odbicie poza niekanonicznymi (marked-site) wektorami wierszy `K_total`; więc nie ma obecnie wyeksportowanego ścisłego znako‑czułego fizycznego datumu orientacji uzyskanego przez “repackaging” istniejących direction-free referencji.

Update: jako kolejny, najwęższy krok `B3` repo eksportuje teraz jawny lane‑scoped **operator transportu chartów** `O_12` pomiędzy `pair1↔pair2`
na `n=12`, wyprowadzony wyłącznie z `alpha_12` (`F461`). To jest “gluing ingredient” tylko na poziomie osi/projektorów (sign‑free),
opakowany twierdzeniem strict (`N506`) i audytowany sondą (`P465`); globalny atlas i globalny transition/gluing object na `C_v1` są już wyeksportowane (`F469/N515`), a globalny projektowy obiekt stanu selektora jest wyeksportowany (`F470/N516`).

Update: `P474` audytuje, że wyeksportowany globalny projektowy obiekt stanu selektora jest spójnie sklejony/transportowany przez globalne operatory przejścia na
poziomie projektorów (ray/projector-level; bez podnoszenia znaku do fizycznego datumu). `N519` pakuje, że residualny znak `Z2` można zamrozić jako gauge dla
wyeksportowanych globalnych obiektów atlasu/przejść/stanu (projective), bez zmiany tych obiektów.

Update: `P475` utrwala decyzję “projekcyjną”: kontynuujemy strict-only domknięcie na semantyce projektorów/span (tam gdzie znak jest gauge‑irrelewant), a gałąź
kierunkowa/znakoczuła (`H37`/`T171`) pozostaje jawnie otwarta jako osobny, przyszły krok.

Update: na bazie `O_12` repo eksportuje teraz jawnie dwu‑chartową **sekcję operatorową** na `{pair1,pair2}` na poziomie projektorowym:
`A_2(pair2) = O_12 A_1(pair1) O_12^T` (`F462`), opakowane (`N507`) i audytowane (`P466`). To nadal nie jest strict-core selector closure i nie rozładowuje globalnie `QW-2191`.

Update: repo eksportuje teraz jawny lane-scoped **atlas selektora** na `{pair1,pair2}` z deklaracją overlap-domain, danymi przejścia (`O_12`) oraz danymi klejenia
(`F463`). Jest to tylko atlas-stub w zadeklarowanym scope; globalny atlas/transition object na `C_v1` jest wyeksportowany osobno (`F469/N515`) na poziomie sekcji projektorowej (bez promocji do operator-level groupoid).

Update: jako kolejny, najwęższy krok `B3` repo eksportuje teraz jawny lane-scoped **trzy‑chartowy** składnik atlasu selektora na `{pair1,pair2,pair3}` na poziomie projektorowym,
wraz z jawnym audytem cocycle/path‑independence na poziomie sklejonej sekcji operatorowej (`F464`), audytowane sondą (`P467`) i opakowane jako twierdzenie strict (`N508`).
To nadal nie jest globalny atlas selektora i nie rozładowuje globalnie `QW-2191`.

Update: jako kolejny, najwęższy krok `B3` repo eksportuje teraz jawny lane-scoped **pięcio‑chartowy** składnik atlasu selektora na `{pair1..pair5}` na poziomie projektorowym (sign‑free),
wraz z jawnym audytem **lokalnych** relacji cocycle/path‑independence dla sąsiednich trójek (1‑2‑3, 2‑3‑4, 3‑4‑5) na poziomie sklejonej sekcji projektorowej (`F465`),
audytowane sondą (`P468`) i opakowane jako twierdzenie strict (`N509`). To nadal nie jest globalny atlas selektora i nie rozładowuje globalnie `QW-2191`.

Update: jako kolejny, najwęższy krok `B3` repo eksportuje teraz dodatkowe operatory transportu **axis-only** na brakujących “długich” krawędziach (`O_14`, `O_15`, `O_25`)
oraz upgrade’uje pięcio‑chartowy składnik `{pair1..pair5}` do jawnych danych cocycle/path‑independence dla **wszystkich** trójek na poziomie sklejonej sekcji projektorowej (`F466`),
audytowane sondą (`P469`) i opakowane jako twierdzenie strict (`N510`). To nadal nie jest globalny atlas selektora i nie rozładowuje globalnie `QW-2191`.

Update: jako kolejny, najwęższy krok `B3` repo eksportuje teraz jawny lane-scoped lift transportu chartów do **oriented** `α mod 2π` na poziomie wektorów reprezentantów
na `{pair1..pair5}`, indukowany przez aktualnie wyeksportowane `u_1..u_5` i jawnie śledzony jako warstwa gauge/konwencji (sign‑tracked) (`F467`).
Sonda `P470` audytuje transport wektorów i pełne relacje trójkowe cocycle/path‑independence na wyeksportowanej instancji, a `N511` pakuje to jako twierdzenie strict
w dyscyplinie “convention layer”. Nie jest to promocja do fizycznej orientacji znaku, nie jest to globalny atlas selektora i nie rozładowuje globalnie `QW-2191`.

Update: pas strict Shannon “element‑order reference” dostarcza teraz jawny wewnętrzny składnik symetrii‑łamacza na pełnym
scaffoldu Fouriera `n=12` (przez `F454`), ale nadal jest to poniżej strict-core selector closure i nie stanowi globalnego
rozładowania `QW-2191`.

Update: theorem-level arytmetyka domyka kryterium defektu element‑order reference: `N514` dowodzi, że `F_k(ord_{Z_n})≠0` dla wszystkich `n≥2` i `k∈{1,…,n-1}`.
Sonda `P461` pozostaje tylko probe-level/regresyjnym checkiem (na `n ∈ {6,8,10,12,14,16,18,20,24}`) i nie promuje `n≠12` do fizycznego scaffoldu `QW-2190`.

Update: dodatkowo, jako ostrożny follow-up, repo eksportuje probe-level kandydat mode-index assignment na `Z_24` (wektory bazowe liczbowe)
indukowany przez tę samą regułę kąta defektu (`P462`). To nie jest promocja `n=24` do `QW-2190`.

Update: repo eksportuje już ścisły upgrade infrastrukturalny scope-extension na typed `Z_24`: datum `r_ord_z24` oraz jawny obiekt mode-index assignment na `Z_24`
(`F468`), opakowane twierdzeniem theorem-level o cięciu `O(2)->Z2` na wszystkich parach `pair_m (m=1..11)` na `Z_24` (`N513`). To nadal jest **nie‑fizyczne**
i nie promuje `n=24` do scaffoldu `QW-2190`.
