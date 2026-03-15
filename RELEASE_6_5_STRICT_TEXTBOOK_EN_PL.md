# RELEASE 6.5 STRICT TEXTBOOK EDITION (EN + PL)

**Version:** 6.5.11  
**Date:** 2026-03-15  
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
  downstream objects where the sign is provably gauge-irrelevant (packaged by `N502`), without implying any global
  sign-sensitive physical orientation datum.
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
2. an axiom-free **global** discharge of `QW-2191` beyond the declared diagonal/local lane and `n=12` scope,
3. a strict physical sign/orientation convention if any downstream claim depends on absolute sign (unless separately
   proven gauge-irrelevant on that lane),
4. ToE closure.

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

## 5) Next honest strict moves (as of 2026-03-15)

Dashboards now recommend `B3` (topological selector bridge continuation) as the next strict move under `QW-2191` discipline (`P438`, `P441`):

1. the residual `Z2` sign freeze (as a tracked gauge/convention layer for exported downstream objects where sign is gauge‑irrelevant) is now
   packaged as a strict theorem (`N502`), and
2. continue strict-only closure under explicit `QW-2191` discipline: export only lane‑scoped transition data when needed
   (e.g. `alpha_12 := (theta_2-\u03b8_1) mod 2\u03c0` for `pair1/pair2` from strict sigma‑int theta supply, `F457`) without implying any global selector transition/gluing object.

Update: the strict Shannon element‑order reference lane now provides one explicit internal symmetry-breaking ingredient
on the full `n=12` Fourier scaffold (via `F454`), but this is still below strict-core selector closure and does not
constitute a global `QW-2191` discharge.

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
    dla których wykazano gauge‑irrelewantność znaku (opakowane przez `N502`), bez promocji do sign-sensitive “fizycznej orientacji”.
19. pochodna dana przejścia `alpha_12 := (theta_2 - theta_1) mod 2π` dla `pair1/pair2` jest wyeksportowana z strict slot‑free theta‑pair na korytarzu sigma‑int (`F457`),
    bez implikowania globalnego selector transition/gluing object.

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
2. aksjomatycznie wolnego **globalnego** rozładowania `QW-2191` poza pasem diagonal/local i zakresem `n=12`,
3. ścisłej konwencji znaku/orientacji tam, gdzie downstream wymaga absolutnego znaku (o ile nie jest osobno pokazane, że znak jest gauge),
4. domknięcia ToE.

## 2) Następny uczciwy ruch (stan: 2026-03-15)

Dashboards teraz rekomendują `B3` (kontynuacja topological selector bridge) jako następny ruch w dyscyplinie `QW-2191` (`P438`, `P441`):

1. zamrożenie residualnego znaku `Z2` (jako śledzona warstwa gauge/konwencji dla obiektów downstream, gdzie znak jest gauge‑irrelewant)
   jest już opakowane jako twierdzenie strict (`N502`), oraz
2. kontynuować strict-only closure w jawnej dyscyplinie `QW-2191`: eksportować tylko lane‑scoped dane przejścia wtedy, gdy są potrzebne
   (np. `alpha_12 := (theta_2-\u03b8_1) mod 2\u03c0` dla `pair1/pair2` z strict theta supply na korytarzu sigma‑int, `F457`) bez implikowania globalnego selector transition/gluing object.

Update: pas strict Shannon “element‑order reference” dostarcza teraz jawny wewnętrzny składnik symetrii‑łamacza na pełnym
scaffoldu Fouriera `n=12` (przez `F454`), ale nadal jest to poniżej strict-core selector closure i nie stanowi globalnego
rozładowania `QW-2191`.
