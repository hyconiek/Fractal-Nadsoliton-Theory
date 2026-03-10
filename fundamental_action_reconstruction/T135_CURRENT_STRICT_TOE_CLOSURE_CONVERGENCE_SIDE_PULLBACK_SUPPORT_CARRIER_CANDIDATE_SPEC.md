# T135 Current Strict ToE Closure Convergence-Side Pullback Support Carrier Candidate Spec

Status: `T135_CURRENT_STRICT_TOE_CLOSURE_CONVERGENCE_SIDE_PULLBACK_SUPPORT_CARRIER_CANDIDATE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

After `T134/F289/N401`, the repo exports one explicit convergence-side
**joint coherence witness candidate** comparing:

1. the provider-object carrier → residual projection lane (`N399`),
2. the sigma-int → residual projection lane (`N400`),

while preserving `N302` (no actual bridge/export-map object support).

The next honest move (still below object support and still selector-neutral)
is to export not only a *witness record*, but one explicit **pullback carrier**
candidate that packages the coherent pair of projections as a typed carrier
object usable in downstream convergence-side reasoning.

This remains strictly below:

1. discharge of `N302`,
2. any export-map object export (`N300`),
3. discharge of `T2`,
4. admissible `S_sel_int` / selector closure / `QW-2191` discharge,
5. ToE closure.

`T135` specifies that pullback carrier candidate form.

## Inputs (already exported)

1. Provider residual projection record (positive-window):
   - `Pi_strict_provider_object_carrier_to_residual_bridge_export_map_object_support_projection_candidate_positive_window_v1`
     from `T132/N399`.
2. Sigma-int residual projection record (positive-window low-delta):
   - `Pi_sigma_int_to_residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_low_delta_v1`
     from `T133/N400`.
3. Joint coherence witness predicate (record-level):
   - `T134` (used only as a predicate shape; no new claims).

## Pullback support carrier candidate (record-level)

Define a candidate pullback carrier record:

```text
PullbackSupportCarrier_pair^{cand} := {
  pair_family: [pair1,pair2],

  provider_projection: <record>,
  sigma_int_projection: <record>,

  theta_provider_cand: {pair1: θ1_prov, pair2: θ2_prov},
  theta_sigma_int_cand:{pair1: θ1_sig,  pair2: θ2_sig },
  theta_deltas:        {pair1: |θ1_prov-θ1_sig|, pair2: |θ2_prov-θ2_sig|},

  eps_theta: <tolerance>,
  coherence_verdict: (theta_deltas[pairi] <= eps_theta for i in {1,2}),

  residual_records: {
    from_provider: <residual_datum_target_slot_candidate_population_record>,
    from_sigma_int:<residual_datum_target_slot_candidate_population_record>
  },

  status: "candidate_only_pullback_support_carrier_below_object_support"
}
```

This carrier is intentionally **non-selecting**:
it does not pick one “true” residual record; it only packages two coherent
candidate residual records and the explicit coherence predicate.

## Candidate object

Export one explicit pullback carrier candidate object:

```text
Mu_strict_convergence_side_pullback_support_carrier_candidate_v1
```

with intended type:

```text
Mu_strict_convergence_side_pullback_support_carrier_candidate_v1 :
  (
    Pi_strict_provider_object_carrier_to_residual_bridge_export_map_object_support_projection_candidate_positive_window_v1,
    Pi_sigma_int_to_residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_low_delta_v1,
    eps_theta
  ) -> PullbackSupportCarrier_pair^{cand}.
```

## Noncyclic / observer-free discipline

`T135` is admissible only if:

1. the pullback carrier introduces no theta inputs and no populated-instance
   inputs (thetas appear only as candidate outputs of the two projections),
2. the pullback carrier introduces no `K_obs`-indexed selection as a primary
   source,
3. `N302` remains in force (this does not claim actual object support).

## Hard limits

`T135` must not claim:

1. discharge of `N302`,
2. actual bridge/export-map object support,
3. actual export-map object export,
4. actual theta export / pair population,
5. admissible `S_sel_int`,
6. strict-core selector closure or `QW-2191` discharge,
7. ToE closure.

