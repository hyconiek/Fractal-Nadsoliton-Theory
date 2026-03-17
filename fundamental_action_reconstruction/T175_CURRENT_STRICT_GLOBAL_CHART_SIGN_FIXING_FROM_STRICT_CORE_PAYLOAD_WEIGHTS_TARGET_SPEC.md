# T175 Current Strict Global Chart Sign Fixing (0-Cochain) From Strict-Core Payload Weights Target Spec

Status: `T175_CURRENT_STRICT_GLOBAL_CHART_SIGN_FIXING_FROM_STRICT_CORE_PAYLOAD_WEIGHTS_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

After `T174` exported an explicit **edgewise** oriented lift datum `s_ij in {+-1}` (a 1-cochain on overlap edges) so that one exported directed representative
transports without sign flips, the next honest continuation toward the `T173` frontier is to stabilize the **chartwise** sign gauge of directed representatives.

`T175` targets an explicit **chart-level** `Z2` sign-fixing datum (a 0-cochain):

```text
t_i in {+-1}  for each chart pair_i
```

constructed deterministically from already exported strict-core per-site payload weights, so that multiple exported directed representatives
become identical after the same sign-fixing rule is applied.

This is a convention-layer gauge fixing step. It does not claim any physical sign datum.

## Scope

`T175` is scoped to:

1. the strict domain `C_v1` charts `{pair1..pair5}`,
2. already exported strict-core per-site payload weights (e.g. `w_break_by_x`, `w_ref_unnormalized_by_x` from `F647`),
3. already exported directed representatives on `C_v1` (premise-based and convention-scoped),
4. a deterministic chartwise sign-fixing rule of the form:

```text
u_i -> sign( <w_i, u_i> ) * u_i
```

with all chosen weights `w_i` explicitly recorded.

`T175` is **below**:

- strict-core physical sign/orientation datum claims,
- `Aut(Z_12)`-invariant sign canonicity (`N462`),
- kernel-alone/global `QW-2191` discharge,
- ToE closure.

## Target objects (what counts as a discharge)

### A) Chart sign-fixing datum (strict_convention)

Export one object containing:

1. the chosen per-chart weights `w_i` (as references to already exported arrays),
2. the computed dot values `<w_i, u_i>`,
3. the resulting `t_i in {+-1}`,
4. the hard limit that this is a convention-layer gauge fixing.

### B) Sign-fixed directed selector state object (strict_convention)

Export one object of the form:

```text
SelectorState_global_C_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1
```

containing the sign-fixed directed vectors `u_1..u_5` on the `n=12` carrier.

### C) Independence audit (probe-level)

Export one probe that:

1. applies the same sign-fixing rule to at least two already exported directed representatives on `C_v1`,
2. verifies the resulting sign-fixed vectors agree on every chart (numeric tolerance),
3. verifies `<w_i, u_i_fixed> > 0` on every chart used by the rule.

### D) Theorem-level packaging (boundary-safe discharge)

Export one theorem-level package that records:

1. a sign-fixed directed representative exists as an exported object (A+B),
2. it is independent of the starting exported representative (C),
3. it is explicitly convention-scoped and does not upgrade directed sign into strict-core physics.

## Acceptance tests (no false pass)

Any `T175` discharge must satisfy:

1. **Explicit data:** the sign-fixing weights and `t_i` are exported as explicit data.
2. **Scope marking:** must be labeled as `strict_convention` (or narrower), not strict-core physics.
3. **Independence audit:** must verify at least two exported directed representatives collapse to the same sign-fixed representative.
4. **Hard limits preserved:** must keep `QW2191_kernel_alone_discharge=false`, no `Aut(Z_12)`-invariant canonicity claim, and no ToE closure claim.

## Hard limits

`T175` must not claim:

1. strict physical sign/orientation datum,
2. strict-core selector closure in the directed/sign-sensitive physical sense,
3. `Aut(Z_12)`-invariant sign canonicity,
4. kernel-alone/global `QW-2191` discharge,
5. ToE closure.

