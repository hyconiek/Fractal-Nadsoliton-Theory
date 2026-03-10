# T127 Current Strict ToE Closure Provider-Object Carrier to Residual Bridge/Export-Map Object-Support Projection Candidate Spec

Status: `T127_CURRENT_STRICT_TOE_CLOSURE_PROVIDER_OBJECT_CARRIER_TO_RESIDUAL_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_PROJECTION_CANDIDATE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

`T125` requires that any honest provider-object carrier layer (if it ever
exists) be **bridge-facing** under the current strict frontier:

```text
ProviderObjectCarrier_pair -> residual bridge/export-map object-support frontier
```

`N302` states that the residual bridge/export-map route is blocked below actual
object support, and `N385` states that one projection layer exists but still
remains below actual object support.

Therefore the next honest move is not an “actual bridge/export-map object
support” claim.

It is weaker:

```text
export one explicit candidate projection interface
from a provider-object carrier candidate
into the residual bridge/export-map object-support frontier,
without theta inputs and without populated-instance inputs.
```

`T127` specifies exactly such a candidate projection, using only already
exported strict-side coupling machinery (`T115`) and a finite noncyclic carrier
built from the provider contraction parameters.

## Inputs

1. Provider-object carrier candidate from `T126`:

   ```text
   ProviderObjectCarrier_pair^{cand,orbit}
   ```

   with associated contraction parameters:

   ```text
   0 < |a| < 1,  0 < |b| < 1.
   ```

2. Strict kernel operational coupling channel from the strict lane:
   `K_strict_gate` parameters (already fixed in the repo).

3. Pair-map-rule candidate form from `T115`:

   ```text
   M_fractal_light_path_pair_map_rule_candidate_v1
   ```

## Candidate projection construction

### 1. Build a finite noncyclic path-ensemble carrier `E_pair`

For each pair slot `i ∈ {1,2}` define:

```text
r_1 := |a|
r_2 := |b|
```

and for each `k ∈ {0,1,...,11}` define:

```text
d_{i,k} := k
w_{i,k} := r_i^k / (sum_{j=0..11} r_i^j)
```

so that:

1. `d_{i,k} >= 0`,
2. `w_{i,k} >= 0`,
3. `sum_k w_{i,k} = 1` (by construction),
4. **noncyclic:** no theta inputs and no populated-instance inputs are used.

This yields a concrete finite carrier:

```text
E_i(pair) := { (d_{i,k}, w_{i,k}) }_{k=0..11}.
```

### 2. Reduce `E_pair` into candidate phases (reuse T115)

Apply the exported strict-side phasor reduction candidate (`T115`) to obtain:

```text
(theta_1^cand, theta_2^cand) and S_orient^cand
```

without promoting the result into an actual theta export.

### 3. Package the result as a bridge-facing projection candidate object

Export one candidate projection object:

```text
Pi_strict_provider_object_carrier_to_residual_bridge_export_map_object_support_projection_candidate_v1
```

with intended type (record-shape reuse; no semantic identification implied):

```text
Pi_strict_provider_object_carrier_to_residual_bridge_export_map_object_support_projection_candidate_v1 :
  (ProviderObjectCarrier_pair^{cand,orbit}, a, b)
    -> residual_datum_bridge_export_map_object_support_projection_candidate_instance
```

where `residual_datum_bridge_export_map_object_support_projection_candidate_instance`
is the same **record class name** already used by the sigma-int projection
candidate lane (`T118`) for the downstream object-support frontier.

with intended meaning:

```text
one candidate-level projection from provider-object carrier data
into a residual-datum bridge/export-map object-support-facing record,
strictly below actual object support and strictly below any export-map object.
```

## Noncyclic / gauge-symmetry discipline

`T127` is admissible only if:

1. no `theta_1,theta_2` are used as inputs,
2. no populated basis-pair instance is used as input,
3. no `K_obs`-indexed selection is used as a primary source,
4. the gauge/symmetry group declared in `T126` is treated as an internal
   symmetry action (not an observer source),
5. no discharge of `N302` is implied.

## Hard limits

`T127` must not claim:

1. actual residual bridge/export-map object support,
2. actual bridge/export-map object export,
3. actual theta export,
4. admissible `S_sel_int`,
5. strict-core selector closure,
6. `QW-2191` discharge,
7. ToE closure.
