# H37: Sign-Distinction State Audit

**Date:** 2026-03-17  
**Status:** `PASS_PREMISE_BASED_DIRECTED_SIGN_DISTINCTION_EXPORTED__PLUS_W_BREAK_ROOTED_DIRECTED_CONVENTION_LIFT_EXPORTED__PLUS_GLOBAL_ORIENTED_EDGE_SIGN_LIFT_EXPORTED__PLUS_GAUGE_EQUIVALENCE_AUDITED__PLUS_GLOBAL_CHART_SIGN_FIXING_EXPORTED__NO_PHYSICAL_SIGN_DATUM`

## Goal

Check whether the current strict core contains any object that identifies `u` and `-u` inside `pair1 = (c_1,s_1)` as physically distinct states rather than coordinate-level representatives of the same undirected axis.

## Inputs

- `H31_PSI0_TO_PAIR1_REDUCTION_AUDIT`
- `H33_PAIR1_SELECTOR_TARGET_JUSTIFICATION_AUDIT`
- `H34_BASIS_COVARIANCE_TARGET_INDEPENDENCE_AUDIT`
- `H35_PAIR1_AXIS_SELECTION_AUDIT`
- `H36_DIRECTED_AXIS_ORIENTATION_AUDIT`
- `F467/N511`: lane-scoped oriented transport (α mod 2π) exists as a **tracked convention layer** (sign-tracked vector section), not a physical sign datum.
- `F470/N516`: global projective selector state object exported on `C_v1` (projector/span; residual sign is gauge at state level).
- `F473/N523`: explicit `Z_12` generator/orientation fixing datum exported as premise-based strict provenance (`T164` discharge; not Aut-invariant by `N462`).
- `F474/N524`: exports a sign-sensitive directed observable `S_dir_pair1_strict_v1` and a global directed selector state object `SelectorState_global_C_v1_directed_strict_v1` (`T171` discharge).
- `F647`: exports a strict-core constructed-source-object witness provider payload including a reflection-breaking per-site weight `w_break_by_x` (nonzero dot against the exported `u_1`), still below admissibility/closure claims.
- `N517` (+ `F471`): even ord-reference weights (`ord_Z12`, `r_ord`) cannot supply a sign-distinction scalar of the form `Σ_x w(x) u_1(x)` on the current exported `pair1` sine axis.
- `N518` (+ `F472`): more generally, any direction-free (`Aut(Z_12)`-invariant) reference weight family `w` cannot supply such a sign-distinction scalar on the current exported `pair1` sine axis (since `-1∈Aut(Z_12)` ⇒ Aut-invariant weights are even under reflection).
- `P472`: a mechanical scan of exported `generated/*.json` artifacts for **weight-like** per-site arrays that break reflection and yield a nonzero scalar `Σ_x w(x)u_1(x)` now reports at least one strictish candidate outside non-canonical marked-site `K_total` rows (probe-level hygiene; heuristic; no promotion).
- `P684`: audits a rooted sign lift on `C_v1` from `w_break` on `pair1` propagated via rooted transports `O_1m`, showing a global directed representative section exists in a tracked convention scope and descends to the strict projective state.
- `P686`: audits the same `w_break`-rooted directed representative against the **full** exported global transition object on `C_v1` (`F469`), confirming edgewise compatibility only **up to sign** and recording which global overlap edges force a sign flip under axis-only (`α mod π`) transport representatives.
- `N686`: packages the `P686` edgewise sign-flip finding as a strict boundary theorem (no physical sign datum promotion).
- `P687`: audits whether the `P686` edgewise sign flips can be eliminated by any per-chart sign relift (`u_i -> t_i u_i`, `t_i∈{±1}`) while keeping the exported transition operators fixed; reports an obstruction.
- `N687`: packages the `P687` non-solvability as a strict boundary theorem (no physical sign datum promotion).
- `F688`: exports an explicit **global oriented edge sign‑lift** object on `C_v1` in `strict_convention` scope: per-edge signs `s_ij∈{±1}` defining an oriented lift `O_ij^(oriented) := s_ij * O_ij^(axis-only)` such that the exported `w_break`‑rooted directed representative transports without sign flips on every edge.
- `P688`: audits that under the `F688` oriented sign‑lift, the exported `w_break`‑rooted directed representative is full‑edge coherent: `(s_ij * O_ij) u_i ≈ u_j` for all 10 edges (no sign flips).
- `N688`: packages `F688/P688` as a theorem-level discharge in the declared convention scope (no physical sign claim).
- `P689`: audits that the `F688` oriented edge sign‑lift pattern is gauge-equivalent (chart-level `Z2` 0‑cochain) to the oriented edge sign‑lift pattern induced by the exported premise-based directed representative `SelectorState_global_C_v1_directed_strict_v1` (no physical sign claim).
- `N689`: packages the `P689` gauge-equivalence result as a boundary-safe theorem: the oriented edge sign‑lift is a convention-layer datum and is not canonical “for free”.
- `T175`: targets an explicit chart-level sign fixing (a `Z2` 0‑cochain) from strict-core payload weights, yielding a sign-fixed directed representative as a tracked convention layer (no physical sign claim).
- `F690`: exports `SelectorState_global_C_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1` (sign-fixed directed representative on `{pair1..pair5}`) and the explicit per-chart sign-fixing data derived from exported strict-core payload weights (`F647`).
- `P690`: audits independence: applying the same chart sign-fix rule to the exported `w_break`-rooted directed representative yields the same sign-fixed directed representative.
- `N690`: packages `F690/P690` as a boundary-safe theorem-level discharge in the declared convention scope (no physical sign claim).
- `F684`: exports `SelectorState_global_C_v1_directed_rooted_transport_from_S_sel_int_w_break_strict_convention_v1` (directed vector representative section on `{pair1..pair5}` in an explicit convention scope; not a physical datum).
- Extension lane note: `AX28/AX29` export an explicit sign-fixing observable on `pair1` and an explicit global oriented vector selector state object on `C_v1` in `strict_extension_only` scope; `P473` audits that this extension-lane oriented vector state is projector-consistent with the strict global projective selector state (`F470`), so it fixes only a sign-gauge representative and does not change strict core.

## Audit

The strict core currently supports:

- a legal coordinate-level representative `u_psi0_pair1`,
- a local chart `pair1 = (c_1,s_1)`,
- a deterministic angle candidate `psi0`.

The strict core **now also supports** (premise-based via `T164`):

- an explicit directed sign-distinction observable on `pair1` (`S_dir_pair1_strict_v1`), whose value flips under `u -> -u`,
- an explicit global **directed** selector state object on `C_v1` (`SelectorState_global_C_v1_directed_strict_v1`) descending to the already exported projective state.

This resolves the original `H37` question in the declared directed scope:

- `u` and `-u` are now physically distinguishable by `S_dir`,
- the directed state object treats residual sign as part of the state (vector-level), not merely gauge.

The strict core still does **not** support:

- any *Aut(Z_12)-invariant* way to obtain such a sign distinction (ruled out as a “for-free” canonicity by `N462`),
- strict-core selector closure (`S_sel_int`) or global discharge of `QW-2191`.

Update (`2026-03-17`): the strict core now also exports an explicit **internal** reflection-breaking per-site weight payload
`w_break_by_x` (via the constructed-source-object witness provider export `F647`). Since `Σ_x w_break(x)u_1(x) ≠ 0` on the exported
`pair1` `u_1`, this supports a rooted sign choice on `pair1` without invoking the premise-based `T164` fixing datum.

`P684` audits that the rooted sign choice can be propagated to `{pair2..pair5}` via the exported rooted transports `O_1m`, producing a
global directed vector representative section that descends to the already exported strict projective state (projector/spans are unchanged under sign).

`F684` exports this directed representative as a global object on `C_v1` in an explicit **strict_convention** scope. This does **not** count
as a strict physical orientation datum and does **not** claim any `Aut(Z_12)`-invariant sign canonicity.

`P686` extends this audit from rooted edges to the full exported global overlap graph on `C_v1` (`F469`): on every exported overlap edge `pairi_to_pairj`,
the axis-only transition operator `O_ij` transports the directed state line correctly (`O_ij u_i ≈ ± u_j`), but several edges force a **sign flip**
(`O_ij u_i ≈ -u_j`) relative to the rooted sign convention. Therefore the exported `w_break`-rooted directed representative remains a section/convention choice
below any strict physical sign datum, and it highlights that axis-only global transport representatives (`α mod π`) cannot canonically support a globally
sign-consistent directed state without an additional oriented lift/convention layer.

`P687` strengthens this further: even allowing an arbitrary per-chart sign relift (`u_i -> t_i u_i`), the full-edge sign-flip pattern under the fixed exported axis-only
transition representatives is **not** eliminable. Therefore, under the fixed `α mod π` transition representatives, no globally edge-sign-consistent directed section exists at all:
achieving directed edgewise sign coherence requires an additional oriented edge-lift / convention layer.

Update (`2026-03-17`): the repo now exports exactly such an explicit oriented edge-lift as a convention layer. `F688` exports a global edge sign‑lift object on `C_v1`
selecting per-edge signs `s_ij∈{±1}` so that the lifted oriented transport `(s_ij * O_ij)` transports the exported `w_break`‑rooted directed representative without sign flips
on every overlap edge. `P688` audits full-edge coherence (all 10 edges) and `N688` packages this as a theorem-level discharge in the declared `strict_convention` scope. This does
**not** upgrade any sign choice into strict-core physics, does **not** claim `Aut(Z_12)`-invariant canonicity, does **not** promote to operator-level groupoid identities (`N512`),
and does **not** imply any kernel-alone/global `QW-2191` discharge.

Update (`2026-03-17`): the oriented edge sign‑lift pattern itself is not canonical. `P689` audits that changing the exported directed representative (from the `w_break`‑rooted
convention representative to the premise-based directed representative) changes the induced edge sign‑lift pattern only by a chart-level `Z2` gauge relift `u_i -> t_i u_i`,
and `N689` packages this as a boundary-safe theorem. This further reinforces that `T174` lives in a tracked convention layer and does not constitute a strict physical sign datum.

Update (`2026-03-17`): the repo now exports an explicit chart-level sign fixing (0‑cochain) as a tracked convention layer (`T175`). `F690` exports a sign-fixed directed representative on
`C_v1` constructed deterministically from already exported strict-core payload weights (including `w_break_by_x` and `w_ref_unnormalized_by_x`), and `P690` audits that the same rule collapses
two exported directed representatives to the same sign-fixed representative (independence), packaged by `N690`. This does **not** promote any strict physical sign datum and does **not** imply kernel-alone/global `QW-2191` discharge.

Update (2026-03-16): after exporting an explicit premise-based fixing datum (`F473/N523`) and an explicit sign-sensitive observable + directed state (`F474/N524`),
the strict core contains a directed sign lift *in the declared premise-based scope*.

Moreover, `N517` records a strict obstruction for one tempting route: the current strict `ord`-reference family (`ord_Z12`, `r_ord`)
is even under reflection and therefore cannot distinguish `u_1` from `-u_1` on the current exported `pair1` sine axis via an observable of the form `Σ_x w(x)u_1(x)`.

`N518` strengthens this to an entire **direction-free** class: any `Aut(Z_12)`-invariant reference weight family (and hence any reference weight not introducing a marked direction/generator) is even under reflection and therefore cannot distinguish sign via such a linear scalar on the current exported sine axis.

Probe-level hygiene: `P472` scans exported `generated/*.json` artifacts and reports strictish reflection-breaking length‑12 arrays with nonzero dot `Σ_x w(x)u_1(x)`,
including at least one strictish **weight-like** candidate outside non-canonical marked-site `K_total` rows on the current export set. This remains a necessary-indicator scan only:
it does not by itself supply a strict physical interpretation, does not imply `Aut(Z_12)`-invariant sign canonicity, and does not imply any `QW-2191` discharge.

## Result

`H37_B4 := strict core exports (i) an explicit sign-sensitive directed observable on pair1 and a global directed selector state object on C_v1 (premise-based via T164), and (ii) an internal w_break-rooted global directed representative section object on C_v1 in an explicit convention scope (F684), both descending to the already exported projective state`

## Hard limits

- No theorem-level PASS.
- No full-closure PASS.
- No claim of `Aut(Z_12)`-invariant canonicity (premise-based `T164` is not Aut-invariant; `w_break`-rooted lift is convention-scoped).
- No claim that `QW-2191` is discharged.
