# N689 Current Strict `T174`: Oriented Edge Sign‑Lift Gauge‑Equivalence Theorem (Boundary‑Safe, No False‑PASS)

Status: `N689_CURRENT_STRICT_T174_ORIENTED_EDGE_SIGN_LIFT_GAUGE_EQUIVALENCE_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-17`

## Claim

On the current repo state, the convention-layer oriented edge sign‑lift pattern exported by `F688` (chosen for the `w_break`‑rooted directed representative)
and the oriented edge sign‑lift pattern induced by the exported premise-based directed representative `SelectorState_global_C_v1_directed_strict_v1`
are **gauge-equivalent** by a chart-level `Z2` 0‑cochain:

```text
s_ij^(B) = t_i * s_ij^(A) * t_j
```

Therefore the `T174` oriented edge sign‑lift remains a tracked **convention layer**:

- it is not canonical “for free” (`N462`),
- it does not upgrade any directed sign into strict-core physics,
- and it does not imply kernel-alone/global `QW-2191` discharge.

## Evidence

- `F688` (export of one explicit oriented edge sign‑lift pattern)
- `P689` (audit of gauge-equivalence across directed representatives)

## Hard limits

This theorem does **not** claim:

- strict-core selector closure in the directed/sign-sensitive physical sense,
- `Aut(Z_12)`-invariant sign canonicity,
- kernel-alone/global `QW-2191` discharge,
- operator-level groupoid identities (`N512` boundary),
- ToE closure.

