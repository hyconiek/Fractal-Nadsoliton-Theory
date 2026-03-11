# T119 Current Strict Sigma-Int To Residual Datum Bridge Export-Map Object-Support Projection Candidate (Positive-Window) Spec

Status: `T119_CURRENT_STRICT_SIGMA_INT_TO_RESIDUAL_DATUM_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_PROJECTION_CANDIDATE_POSITIVE_WINDOW_SPEC_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

`T118/F271/N383` exported one explicit **candidate** object-to-map support
projection artifact from `sigma_int_input ∈ {+1,-1}` into a residual-datum target-slot
candidate population record.

The honest weakness of that construction is also explicit in `T115`:

```text
(X_i^cand, Y_i^cand) = (0,0) is a real degeneracy frontier
```

`T119` strengthens the projection **without** claiming:

1. discharge of `T2`,
2. discharge of `N302`,
3. actual bridge/export-map object support,
4. actual `theta_1`, `theta_2`,
5. selector closure or ToE closure,

by introducing one narrow, typed **positive-window corridor** that guarantees:

```text
cos(Theta(d)) > 0 and sin(Theta(d)) > 0 on every path in E_pair
=> X_i^cand > 0 and Y_i^cand > 0
=> atan2 is well-defined (no degeneracy)
```

This keeps the move:

- noncyclic (no theta input, no populated-instance input),
- observer-free (no `K_obs` primary source),
- and still candidate-only.

## Positive-window corridor (strict-side)

Reuse strict kernel notation from Release 6.2 Strict:

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

Define the phase barrier margin and its local half-margin:

```math
\delta^{barrier}:=\frac{\pi}{2}-|\phi|>0,
\qquad
\varepsilon^{local}:=\frac{1}{2}\delta^{barrier}.
```

Define the strict positive-window corridor on the `d` axis:

```math
d^{local}:=\frac{\varepsilon^{local}}{\omega}.
```

Then for all `d ∈ [0, d^{local}]`:

```math
0<\Theta(d)<\frac{\pi}{2}
\quad\Longrightarrow\quad
\cos(\Theta(d))>0,\ \sin(\Theta(d))>0.
```

## Candidate generator (delta-scaled nad12 carrier)

Inputs:

1. `sigma_int_input ∈ {+1,-1}`,
2. `eps ∈ [0,1]` (candidate parameter),
3. `delta_d ∈ (0, d^{local}/11]` (candidate step; enforces positive-window),
4. fixed masks:

   ```math
   b_{1,k}:=(-1)^k,\qquad b_{2,k}:=(-1)^{k+1},
   \qquad k\in\{0,1,\dots,11\}.
   ```

Notation: in the formulas below, `\sigma_{int}^{in}` denotes the chosen input
value `sigma_int_input`.

Define the delta-scaled distances and weights:

```math
d_{i,k}:=k\,\delta_d,
\qquad
w_{i,k}:=\frac{1+\sigma_{int}^{in}\,\varepsilon\,b_{i,k}}{12}.
```

Define:

```text
E_i(pair) := { (d_{i,k}, w_{i,k}) }_{k=0..11}.
```

Because `delta_d <= d^{local}/11`, we have `d_{i,k} <= d^{local}` for all
`k=0..11`. Therefore every term in the phasor reduction is strictly positive
in both its `cos` and `sin` parts, hence:

```math
X_i^{cand}>0,\qquad Y_i^{cand}>0,
```

and `atan2(Y_i^{cand},X_i^{cand})` is well-defined.

## Projection object (positive-window variant)

Export one explicit projection candidate:

```text
Pi_sigma_int_to_residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_v1
```

with intended type:

```text
Pi_sigma_int_to_residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_v1 :
  (sigma_int_input, eps, delta_d)
    -> residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_instance
```

Construction:

1. build `E_pair` via the delta-scaled nad12 generator above,
2. reduce via `M_fractal_light_path_pair_map_rule_candidate_v1` (`T115`),
3. emit a residual-datum target-slot **candidate** population record below `R1`.

## Hard limits

`T119` must not claim:

1. discharge of `T2`,
2. discharge of `N302`,
3. satisfaction of the `N301` bridge/export-map object target,
4. actual residual bridge/export-map object support,
5. actual `theta_1`, `theta_2`,
6. actual pair population,
7. admissible `S_sel_int`,
8. strict-core selector closure or `QW-2191` discharge,
9. ToE closure.
