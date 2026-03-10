# T132 Current Strict ToE Closure Provider-Object Carrier to Residual Bridge/Export-Map Object-Support Projection Candidate (Positive-Window) Spec

Status: `T132_CURRENT_STRICT_TOE_CLOSURE_PROVIDER_OBJECT_CARRIER_TO_RESIDUAL_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_PROJECTION_CANDIDATE_POSITIVE_WINDOW_SPEC_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

`T127/F280/N392` exports one explicit **provider-carrier → residual** projection
candidate, but that base form uses the unscaled orbit-depth distances:

```text
d_{i,k} := k
```

which exceeds the strict kernel local phase-barrier corridor from Release 6.2
Strict (and from `T119`), therefore it does **not** guarantee:

```text
cos(Theta(d))>0 and sin(Theta(d))>0 on all carrier paths
=> X_i^cand>0 and Y_i^cand>0
=> no atan2 degeneracy
```

`T132` strengthens the provider-carrier projection candidate **without**
claiming:

1. discharge of `N302`,
2. actual bridge/export-map object support,
3. actual theta export,
4. admissible `S_sel_int`,
5. selector closure or `QW-2191` discharge,
6. ToE closure,

by introducing one typed **positive-window corridor** restriction:

```text
d_{i,k} := k*delta_d,  with delta_d <= d_local/11
```

so that every path distance stays in `[0,d_local]`, where the strict kernel has
`cos(Theta(d))>0` and `sin(Theta(d))>0`.

This remains:

- noncyclic (no `theta` input; no populated-instance input),
- observer-free (no `K_obs` as a primary source),
- and candidate-only (no object-support discharge).

## Positive-window corridor (strict-side)

Reuse strict kernel notation:

```math
\Theta(d):=\omega d+\phi,
\qquad
D(d):=1+\beta d^{\eta}.
```

Assume the current strict working tuple has:

```math
0<\phi<\frac{\pi}{2},
\qquad
\omega>0,
\qquad
\beta>0,
\qquad
\eta>1.
```

Define:

```math
\delta^{barrier}:=\frac{\pi}{2}-|\phi|>0,
\qquad
\varepsilon^{local}:=\frac{1}{2}\delta^{barrier},
\qquad
d^{local}:=\frac{\varepsilon^{local}}{\omega}.
```

Then for all `d ∈ [0,d^{local}]`:

```math
0<\Theta(d)<\frac{\pi}{2}
\quad\Longrightarrow\quad
\cos(\Theta(d))>0,\ \sin(\Theta(d))>0.
```

Define the 12-slot step constraint:

```math
\delta_d \in (0, d^{local}/11].
```

Then for `k ∈ {0,...,11}` we have `k \delta_d <= d^{local}`.

## Candidate projection object (positive-window variant)

Export one explicit positive-window projection candidate:

```text
Pi_strict_provider_object_carrier_to_residual_bridge_export_map_object_support_projection_candidate_positive_window_v1
```

with intended type:

```text
Pi_strict_provider_object_carrier_to_residual_bridge_export_map_object_support_projection_candidate_positive_window_v1 :
  (ProviderObjectCarrier_pair^{cand,orbit}, a, b, delta_d)
    -> residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_instance
```

## Construction (delta-scaled orbit-depth carrier)

Inputs:

1. provider-object carrier candidate from `T126`,
2. contraction parameters `a,b` with `0<|a|<1`, `0<|b|<1`,
3. `delta_d ∈ (0, d^{local}/11]` (typed positive-window step).

Define:

```text
r_1 := |a|
r_2 := |b|
```

For each `k ∈ {0,1,...,11}` define:

```text
d_{i,k} := k*delta_d
w_{i,k} := r_i^k / (sum_{j=0..11} r_i^j)
```

so that:

1. `d_{i,k} ∈ [0,d^{local}]`,
2. `w_{i,k} >= 0`,
3. `sum_k w_{i,k} = 1`.

Apply the strict phasor reduction candidate `T115`:

```text
(theta_1^cand, theta_2^cand) := M_fractal_light_path_pair_map_rule_candidate_v1(E_pair)
```

and emit a residual-datum target-slot **candidate** population record:

```text
residual_datum_target_slot_candidate_population_record :=
(
  theta_1^cand, theta_2^cand,
  u_1^cand = cos(theta_1^cand)c_1 + sin(theta_1^cand)s_1,
  u_2^cand = cos(theta_2^cand)c_2 + sin(theta_2^cand)s_2,
  S_orient^cand = span{u_1^cand,u_2^cand}
).
```

Because every `d_{i,k}` lies in the positive-window corridor, every term in
the phasor sums has strictly positive `cos(Theta(d))/D(d)` and
`sin(Theta(d))/D(d)`. Therefore:

```math
X_i^{cand}>0,\qquad Y_i^{cand}>0,
```

so `atan2(Y_i^{cand},X_i^{cand})` is well-defined (no degeneracy).

## Hard limits

`T132` must not claim:

1. discharge of `N302`,
2. actual bridge/export-map object support,
3. actual export-map object export,
4. actual theta export / pair population,
5. admissible `S_sel_int`,
6. strict-core selector closure or `QW-2191` discharge,
7. ToE closure.

