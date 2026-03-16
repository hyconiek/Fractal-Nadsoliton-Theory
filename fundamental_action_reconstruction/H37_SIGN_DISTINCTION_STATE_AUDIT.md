# H37: Sign-Distinction State Audit

**Date:** 2026-03-16  
**Status:** `PASS_STRICT_DIRECTED_SIGN_DISTINCTION_OBSERVABLE_EXPORTED_AND_GLOBAL_DIRECTED_SELECTOR_STATE_OBJECT_EXPORTED__PREMISE_BASED_T164`

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
- `N517` (+ `F471`): even ord-reference weights (`ord_Z12`, `r_ord`) cannot supply a sign-distinction scalar of the form `Σ_x w(x) u_1(x)` on the current exported `pair1` sine axis.
- `N518` (+ `F472`): more generally, any direction-free (`Aut(Z_12)`-invariant) reference weight family `w` cannot supply such a sign-distinction scalar on the current exported `pair1` sine axis (since `-1∈Aut(Z_12)` ⇒ Aut-invariant weights are even under reflection).
- `P472`: a mechanical scan of exported `generated/*.json` artifacts for **weight-like** per-site arrays that break reflection and yield a nonzero scalar `Σ_x w(x)u_1(x)` reports none in strict(-derived) scope **outside** non-canonical marked-site `K_total` row vectors (probe-level hygiene; no promotion).
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

Update (2026-03-16): after exporting an explicit premise-based fixing datum (`F473/N523`) and an explicit sign-sensitive observable + directed state (`F474/N524`),
the strict core contains a directed sign lift *in the declared premise-based scope*.

Moreover, `N517` records a strict obstruction for one tempting route: the current strict `ord`-reference family (`ord_Z12`, `r_ord`)
is even under reflection and therefore cannot distinguish `u_1` from `-u_1` on the current exported `pair1` sine axis via an observable of the form `Σ_x w(x)u_1(x)`.

`N518` strengthens this to an entire **direction-free** class: any `Aut(Z_12)`-invariant reference weight family (and hence any reference weight not introducing a marked direction/generator) is even under reflection and therefore cannot distinguish sign via such a linear scalar on the current exported sine axis.

Probe-level hygiene: `P472` additionally reports that among exported strict(-derived) artifacts, the only **weight-like** reflection-breaking per-site arrays producing a nonzero `Σ_x w(x)u_1(x)` are non-canonical marked-site rows of an exported `K_total` matrix. In particular, it reports **zero** strict(-derived) weight-like candidates outside those `K_total` rows. This does not prove impossibility, but it supports the audit conclusion that discharging `H37` requires an explicit reflection-breaking/orientation source (e.g. a generator/orientation fixing datum) rather than a repackaging of currently exported direction-free references.

## Result

`H37_B3 := strict core exports an explicit sign-sensitive directed observable on pair1 and an explicit global directed selector state object on C_v1 (premise-based via T164), descending to the already exported projective state`

## Hard limits

- No theorem-level PASS.
- No full-closure PASS.
- No claim of `Aut(Z_12)`-invariant canonicity (the fixing datum is premise-based; `N462`).
- No claim that `QW-2191` is discharged.
