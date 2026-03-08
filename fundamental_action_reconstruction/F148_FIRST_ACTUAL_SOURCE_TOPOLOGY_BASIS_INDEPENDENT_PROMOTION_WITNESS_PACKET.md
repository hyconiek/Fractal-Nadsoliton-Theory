# F148 First Actual Source Topology Basis Independent Promotion Witness Packet

Status: `F148_EXECUTED_FIRST_ACTUAL_SOURCE_TOPOLOGY_BASIS_INDEPENDENT_PROMOTION_WITNESS_PACKET_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

After `F146/P234/N254` and `F147/P235/N255`, the current repo state already
exports:

1. one actual full source-topology nontriviality witness for
   `tau_src_candidate_v1`,
2. one actual source-side selector witness for `tau_src_candidate_v1`.

But it still did not yet export:

1. an actual basis-independent selector-promotion witness,
2. a quotient-safe `QW-2191` resolution,
3. a current selector closure proof.

The next honest move is narrower:

```text
lift the current actual chart-bound selector witness
to one actual basis-free selector datum packet
using only internal class reduction on the already exported preobserver lane
strictly below quotient-safe QW-2191 resolution
and strictly below current selector closure
```

`F148` executes exactly that move.

## Fixed input

Reuse the future-only basis-independent promotion target from `F136`:

```text
Upsilon_sel_basis_target_v1 :
(Theta_src_nontriv_discharge_target_v1, Pi_sel_src_target_v1)
  -> Sigma_sel_basis_free_target_v1
```

Reuse the actual full source-topology nontriviality witness from `F146`:

```text
Theta_src_nontriv_actual_discharge_witness_v1 :
tau_src_candidate_v1 -> actual_full_source_topology_nontriviality_discharge_target_v1
```

Reuse the actual source-side selector witness from `F147`:

```text
Pi_sel_src_actual_witness_v1 :
tau_src_candidate_v1 -> Sigma_sel_src_target_v1
```

## Internal class reduction

Freeze one explicit current-repo-state class-reduction operator:

```text
Q_basis_sel_v1 :
Sigma_sel_src_target_v1 -> Sigma_sel_basis_free_target_v1
```

with the following intended action:

1. map the chart-bound selector axis realization to one basis-free axis class
   by forgetting frame labels and keeping only the normalized rank-one ray
   projector class,
2. map the chart-bound signed selector split to one basis-free signed-split
   class by keeping only the complementary projector pair together with the
   positive source-alignment witness and the vanishing minus-channel witness,
3. map the chart-bound preobserver scope realization to one basis-free
   preobserver scope tag by keeping only preobserver-only and
   observer-downstream-only scope data,
4. remove chart labels without using observer-side asymmetry as the source of
   selector asymmetry,
5. remain strictly below quotient-safe `QW-2191` resolution.

So inside the current repo state:

```text
selector_axis_basis_free_class_v1
selector_signed_split_basis_free_class_v1
preobserver_basis_free_scope_tag_v1
```

are treated as class-reduced images of the already exported chart-bound
selector witness, not as fresh observer-driven imports.

## Actual basis-independent promotion lift

Freeze one explicit support packet:

```text
W_src_basis_promotion_support_packet_v1 :=
(
  tau_src_candidate_v1,
  Theta_src_nontriv_actual_discharge_witness_v1,
  Pi_sel_src_actual_witness_v1,
  Q_basis_sel_v1,
  selector_axis_basis_free_class_v1,
  selector_signed_split_basis_free_class_v1,
  preobserver_basis_free_scope_tag_v1,
  observer_downstream_only
)
```

Interpretation on the current repo state:

1. `Theta_src_nontriv_actual_discharge_witness_v1` supplies the actual source
   nontriviality support required by `F136`,
2. `Pi_sel_src_actual_witness_v1` supplies the actual source-side selector
   datum required by `F136`,
3. `Q_basis_sel_v1` removes the dependence on the admissible chart labels by
   passing from realized selector data to basis-free classes,
4. `selector_axis_basis_free_class_v1` is represented by the normalized
   rank-one ray projector induced by `E_orient_preLM_v1`,
5. `selector_signed_split_basis_free_class_v1` is represented by the
   complementary selector projector pair induced by `B_sel_preLM_v1` together
   with the positive source response carried through `R_sel_preLM_v1` and
   `O_sel_preLM_v1`,
6. `preobserver_basis_free_scope_tag_v1` keeps the witness inside an explicit
   preobserver scope with downstream-only observer role,
7. this step does not identify `tau_src_candidate_v1` with
   `S_preLM_strict_core_source_object_v1`,
8. this step does not yet claim quotient-safe uniqueness across the
   `QW-2191` ambiguity family.

Therefore freeze one actual basis-independent promotion witness:

```text
Upsilon_sel_basis_actual_witness_v1 :
tau_src_candidate_v1 -> Sigma_sel_basis_free_target_v1
```

with current-repo-state support packet:

```text
Upsilon_sel_basis_actual_witness_v1
  := W_src_basis_promotion_support_packet_v1
```

## Meaning of the witness

This witness is intended only as:

1. one actual source-side basis-independent selector-promotion witness for
   `tau_src_candidate_v1`,
2. one actual refinement of the future-only promotion target route from
   `F136`,
3. one class-level lift from a chart-bound selector realization to a
   basis-free selector datum packet,
4. a current positive witness strictly before quotient-safe `QW-2191`
   resolution and current selector closure.

It is not yet:

1. a quotient-safe `QW-2191` witness,
2. a current selector closure proof,
3. a current global `QW-2191` discharge.

## Why this is the honest lift

`F148` is the narrowest honest lift because:

1. `F136` already froze the exact future-only basis-free codomain packet
   `Sigma_sel_basis_free_target_v1`,
2. `F146` already supplied one actual full source-topology nontriviality
   witness,
3. `F147` already supplied one actual source-side selector witness,
4. `N196/N197/N198/N199` already supplied the admissible preobserver
   orientation/bridge/reduction/output chain needed to define the selector
   classes,
5. the present step adds only the basis-class reduction from that already
   exported selector witness into the already frozen basis-free packet,
6. it does not claim quotient-safe `QW-2191` resolution or selector closure.

## Why this is still kernel-split-safe

`F148` remains kernel-split-safe because:

1. it uses only already exported strict-kernel source-limit witness data,
2. it reuses only already exported preobserver selector operators,
3. it does not identify `K_strict_gate` with `K_legacy_ont`,
4. it does not transfer any legacy physical-role semantics,
5. it does not claim a legacy-to-strict bridge.

## Observer role

Observer remains outside the witness domain:

1. the witness is exported from source-side nontriviality and source-side
   selector data,
2. observer remains downstream only,
3. no observer-side asymmetry is used as the selector source,
4. observer may appear later only as a downstream pushforward witness.

This preserves the declared order:

```text
nadsoliton -> light -> matter -> emergent observer
```

## Result

`F148` exports one actual source-side basis-independent selector-promotion
witness:

```text
Upsilon_sel_basis_actual_witness_v1 :
tau_src_candidate_v1 -> Sigma_sel_basis_free_target_v1
```

supported by:

```text
W_src_basis_promotion_support_packet_v1
```

with the declared properties:

1. actual basis-independent selector-promotion witness,
2. source-side only,
3. observer-free in the witness domain,
4. derived by internal class reduction from the actual selector witness,
5. below quotient-safe `QW-2191` resolution,
6. below current selector closure,
7. below current global `QW-2191` discharge,
8. no false pass.

## Hard limits

`F148` does not discharge:

1. quotient-safe `QW-2191` resolution,
2. current selector closure,
3. current global `QW-2191` discharge,
4. ToE closure.

## Recommended next move

The correct next move is:

1. test whether the current repo really exports this actual basis-independent
   selector-promotion witness in a guardrail-consistent way,
2. then attempt one actual quotient-safe `QW-2191` resolution witness from
   the basis-free selector packet,
3. only after that discuss any stronger selector-closure language.
