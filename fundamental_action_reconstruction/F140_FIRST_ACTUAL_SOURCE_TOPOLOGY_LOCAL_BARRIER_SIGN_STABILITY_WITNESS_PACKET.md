# F140 First Actual Source Topology Local Barrier Sign Stability Witness Packet

Status: `F140_EXECUTED_FIRST_ACTUAL_SOURCE_TOPOLOGY_LOCAL_BARRIER_SIGN_STABILITY_WITNESS_PACKET_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

After `F139/P227/N247`, the next honest move is still not:

1. a full barrier-protected sign discharge,
2. a full source-topology nontriviality discharge,
3. a basis-independent selector promotion discharge,
4. a quotient-safe `QW-2191` resolution,
5. a current selector closure proof.

It is narrower:

```text
export one actual source-side local sign-stability witness
on an explicit positive-radius neighborhood of the declared core branch
strictly below full barrier-protected sign discharge
and strictly below QW-2191 discharge
```

`F140` executes exactly that move.

## Fixed input

Reuse the candidate packet from `F127`:

```text
tau_src_candidate_v1 :=
(
  source_limit_tag_v1,
  phi_barrier_tag_v1,
  T_flow^(0)
)
```

Reuse the exported strict-kernel core datum from `F74`:

```text
phi = 0.16250
```

Reuse the actual scalar sign witness and barrier margin from `F139`:

```text
delta_src_barrier_sign_margin_v1 := (pi/2) - |phi| = 1.4082963267948965 > 0
psi_src_barrier_sign_component_witness_v1 := sign(cos(phi)) = +1
```

## Actual local stability witness

Freeze one explicit local branch radius:

```text
epsilon_src_local_barrier_radius_v1 := delta_src_barrier_sign_margin_v1 / 2
```

Numerically:

```text
epsilon_src_local_barrier_radius_v1 = 0.7041481633974482 > 0
```

Then for every real perturbation `epsilon` satisfying:

```text
|epsilon| <= epsilon_src_local_barrier_radius_v1
```

we have:

```text
|phi + epsilon|
  <= |phi| + |epsilon|
  <= |phi| + delta_src_barrier_sign_margin_v1 / 2
  < pi/2
```

and therefore:

```text
cos(phi + epsilon) > 0
sign(cos(phi + epsilon)) = +1
```

Freeze one actual local barrier-sign stability witness:

```text
chi_src_local_barrier_sign_stability_witness_v1 :
for all epsilon in R,
if |epsilon| <= epsilon_src_local_barrier_radius_v1,
then sign(cos(phi + epsilon)) = +1
```

## Meaning of the witness

This witness is intended only as:

1. one actual source-side local sign-stability witness on the declared
   source-limit core branch,
2. one actual positive-radius branch witness strictly stronger than a single
   pointwise sign component,
3. a current positive witness still strictly below full
   `barrier_protected_sign_class_v1`,
4. a current positive witness strictly before full source-topology
   nontriviality, basis-independent selector promotion, and quotient-safe
   `QW-2191` resolution.

It is not yet:

1. a full barrier-protected sign discharge for `tau_src_candidate_v1`,
2. a full source-topology nontriviality discharge,
3. a basis-independent selector witness,
4. a quotient-safe `QW-2191` witness,
5. a current selector closure proof.

## Why this is still below full barrier-protected sign discharge

`F140` remains strictly below full barrier-protected sign discharge because:

1. it proves only a local sign-stability neighborhood on the declared current
   core branch,
2. it does not yet promote that local witness to the full abstract
   `barrier_protected_sign_class_v1`,
3. it does not yet assemble sign with actual observer-free scope into a full
   source-topology nontriviality discharge,
4. it does not claim any quotient-safe selector consequence.

## Why this is still kernel-split-safe

`F140` remains kernel-split-safe because:

1. it uses only the already exported strict-kernel source-limit control datum,
2. it uses the already exported `phi_barrier_tag_v1` only on the declared
   current core branch,
3. it does not identify `K_strict_gate` with `K_legacy_ont`,
4. it does not transfer any legacy physical-role semantics,
5. it does not claim a legacy-to-strict bridge.

## Observer role

Observer remains outside the witness domain:

1. the witness is extracted before any observer-side pushforward,
2. observer remains downstream only,
3. no observer-side algebra is used as proof of local barrier sign stability.

This preserves the declared order:

```text
nadsoliton -> light -> matter -> emergent observer
```

## Result

`F140` exports one explicit positive local radius and one actual local
barrier-sign stability witness:

```text
epsilon_src_local_barrier_radius_v1 := 0.7041481633974482 > 0

chi_src_local_barrier_sign_stability_witness_v1 :
for all epsilon in R,
if |epsilon| <= epsilon_src_local_barrier_radius_v1,
then sign(cos(phi + epsilon)) = +1
```

with the declared properties:

1. actual positive local radius,
2. actual local sign-stability witness,
3. source-side only,
4. observer-free in the witness domain,
5. below full barrier-protected sign discharge,
6. below full source-topology nontriviality discharge,
7. below basis-independent selector promotion,
8. below quotient-safe `QW-2191` resolution,
9. below current selector closure,
10. no false pass.

## Hard limits

`F140` does not discharge:

1. full barrier-protected sign of `tau_src_candidate_v1`,
2. full source-topology nontriviality,
3. basis-independent selector promotion,
4. quotient-safe `QW-2191` resolution,
5. current selector closure,
6. current global `QW-2191` discharge,
7. ToE closure.

## Recommended next move

The correct next move is:

1. test whether the current repo really exports this actual local barrier-sign
   stability witness in a guardrail-consistent way,
2. keep the result strictly below full barrier-protected sign discharge,
3. only then attempt a branch-to-packet lift toward an actual
   `barrier_protected_sign_class_v1` witness.
