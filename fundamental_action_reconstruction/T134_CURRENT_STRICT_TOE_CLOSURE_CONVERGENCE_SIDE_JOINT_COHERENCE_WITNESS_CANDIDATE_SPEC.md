# T134 Current Strict ToE Closure Convergence-Side Joint Coherence Witness Candidate Spec

Status: `T134_CURRENT_STRICT_TOE_CLOSURE_CONVERGENCE_SIDE_JOINT_COHERENCE_WITNESS_CANDIDATE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

`T114` names a future-only convergence-side support target immediately below
`N378`. One acceptance component is:

```text
Zeta_strict_convergence_side_joint_coherence_target_v1
```

which demands an explicit coherence predicate binding the provider-object
carrier side and the internal residual-datum side **on the same pair index**
without claiming:

- actual provider-object realization,
- actual internal-orientation realization,
- actual `E_orient`,
- admissible `S_sel_int`,
- selector closure or `QW-2191` discharge,
- ToE closure.

On the current repo state, the strongest honest convergence-side move remains
below any such realizations. The admissible move is narrower:

```text
export one explicit joint coherence *witness candidate*
at the level of residual-datum projection records,
binding two noncyclic, observer-free projection lanes.
```

`T134` specifies such a witness candidate form.

## Inputs (already exported lanes)

1. Provider-object carrier projection (positive-window):
   - `N399` exports a provider-carrier → residual projection candidate
     protected against `atan2` degeneracy.
2. Sigma-int projection (positive-window low-delta):
   - `N400` exports a sigma-int → residual projection candidate protected
     against `atan2` degeneracy with a low-delta calibration.
3. Both projections land in the same residual-frontier record shape used on
   the strict residual-datum lane (record-class compatibility is an explicit
   contract in `T118/T127`).
4. `N302` remains in force:
   - no actual bridge/export-map object support is exported.

## Joint coherence witness candidate (record-level)

Define a candidate joint coherence witness record:

```text
JointCoherenceWitness_pair^{cand} := {
  pair_family: [pair1,pair2],
  provider_projection: <record>,
  sigma_int_projection: <record>,
  theta_deltas: {
    pair1: |theta_1^cand(prov) - theta_1^cand(sig)|,
    pair2: |theta_2^cand(prov) - theta_2^cand(sig)|
  },
  eps_theta: <tolerance>,
  coherence_verdict: (theta_deltas[pairi] <= eps_theta for i in {1,2}),
  status: "candidate_only_joint_coherence_witness_below_object_support"
}
```

This witness is admissible only as:

1. a **candidate** convergence-side coherence witness,
2. built from already-exported noncyclic, observer-free projection records,
3. explicitly below actual object support (`N302`),
4. explicitly neutral w.r.t. selector closure.

## Candidate object

Export one explicit coherence witness candidate object:

```text
Eta_strict_convergence_side_joint_coherence_witness_candidate_v1
```

with intended type:

```text
Eta_strict_convergence_side_joint_coherence_witness_candidate_v1 :
  (
    Pi_strict_provider_object_carrier_to_residual_bridge_export_map_object_support_projection_candidate_positive_window_v1,
    Pi_sigma_int_to_residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_low_delta_v1,
    eps_theta
  ) -> JointCoherenceWitness_pair^{cand}
```

## Hard limits

`T134` must not claim:

1. discharge of `N302`,
2. actual bridge/export-map object support,
3. actual bridge/export map export,
4. actual theta export / pair population,
5. admissible `S_sel_int`,
6. strict-core selector closure or `QW-2191` discharge,
7. ToE closure.

