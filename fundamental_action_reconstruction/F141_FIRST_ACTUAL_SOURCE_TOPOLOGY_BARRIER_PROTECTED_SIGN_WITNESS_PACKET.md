# F141 First Actual Source Topology Barrier-Protected Sign Witness Packet

Status: `F141_EXECUTED_FIRST_ACTUAL_SOURCE_TOPOLOGY_BARRIER_PROTECTED_SIGN_WITNESS_PACKET_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

After `F140/P228/N248`, the next honest move is still not:

1. a full source-topology nontriviality discharge,
2. a basis-independent selector promotion discharge,
3. a quotient-safe `QW-2191` resolution,
4. a current selector closure proof.

It is narrower:

```text
lift the current branch-level barrier-sign support
to one actual source-side barrier-protected sign witness packet
for tau_src_candidate_v1
strictly below full source-topology nontriviality
and strictly below QW-2191 discharge
```

`F141` executes exactly that move.

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

Reuse the future-only target from `F131`:

```text
Psi_src_barrier_sign_target_v1 : tau_src_candidate_v1 -> barrier_protected_sign_class_v1
```

Reuse the actual sign support from `F139`:

```text
delta_src_barrier_sign_margin_v1 := (pi/2) - |phi| = 1.4082963267948965 > 0
psi_src_barrier_sign_component_witness_v1 := sign(cos(phi)) = +1
```

Reuse the actual local stability witness from `F140`:

```text
epsilon_src_local_barrier_radius_v1 := delta_src_barrier_sign_margin_v1 / 2

chi_src_local_barrier_sign_stability_witness_v1 :
for all epsilon in R,
if |epsilon| <= epsilon_src_local_barrier_radius_v1,
then sign(cos(phi + epsilon)) = +1
```

## Branch-to-packet lift

Freeze one explicit support packet:

```text
W_src_barrier_sign_support_packet_v1 :=
(
  phi_barrier_tag_v1,
  delta_src_barrier_sign_margin_v1,
  epsilon_src_local_barrier_radius_v1,
  psi_src_barrier_sign_component_witness_v1,
  chi_src_local_barrier_sign_stability_witness_v1
)
```

Interpretation on the current repo state:

1. `phi_barrier_tag_v1` fixes the relevant source-side barrier tag,
2. `delta_src_barrier_sign_margin_v1 > 0` shows noncontact with the first
   sign-flip boundary,
3. `psi_src_barrier_sign_component_witness_v1 = +1` gives the declared core
   sign,
4. `chi_src_local_barrier_sign_stability_witness_v1` promotes that sign from a
   pointwise fact to a positive-radius local branch fact.

Therefore freeze one actual branch-to-packet lift witness:

```text
Psi_src_barrier_sign_actual_witness_v1 :
tau_src_candidate_v1 -> barrier_protected_sign_class_v1
```

with current-repo-state support packet:

```text
Psi_src_barrier_sign_actual_witness_v1
  := W_src_barrier_sign_support_packet_v1
```

## Meaning of the witness

This witness is intended only as:

1. one actual source-side barrier-protected sign witness for
   `tau_src_candidate_v1`,
2. one actual lift from the previously exported local branch witness to the
   abstract sign class already targeted in `F131`,
3. a current positive witness strictly before full source-topology
   nontriviality, basis-independent selector promotion, and quotient-safe
   `QW-2191` resolution.

It is not yet:

1. a full source-topology nontriviality discharge,
2. a basis-independent selector witness,
3. a quotient-safe `QW-2191` witness,
4. a current selector closure proof.

## Why this is the honest lift

`F141` is the narrowest honest lift because:

1. `F131` already froze the exact codomain
   `barrier_protected_sign_class_v1`,
2. `F139` already supplied current scalar sign and positive barrier margin,
3. `F140` already supplied a positive-radius local sign-stability witness on
   the declared core branch,
4. the present step adds only the packet-level lift from those already
   exported support objects into the already declared sign class,
5. it does not claim full source-topology nontriviality or any selector
   consequence.

## Why this is still kernel-split-safe

`F141` remains kernel-split-safe because:

1. it uses only the already exported strict-kernel source-limit control datum,
2. it uses the already exported `phi_barrier_tag_v1` only on the declared
   current core branch,
3. it does not identify `K_strict_gate` with `K_legacy_ont`,
4. it does not transfer any legacy physical-role semantics,
5. it does not claim a legacy-to-strict bridge.

## Observer role

Observer remains outside the witness domain:

1. the witness is exported before any observer-side pushforward,
2. observer remains downstream only,
3. no observer-side algebra is used as proof of barrier-protected sign.

This preserves the declared order:

```text
nadsoliton -> light -> matter -> emergent observer
```

## Result

`F141` exports one actual source-side barrier-protected sign witness:

```text
Psi_src_barrier_sign_actual_witness_v1 :
tau_src_candidate_v1 -> barrier_protected_sign_class_v1
```

supported by:

```text
W_src_barrier_sign_support_packet_v1 :=
(
  phi_barrier_tag_v1,
  delta_src_barrier_sign_margin_v1,
  epsilon_src_local_barrier_radius_v1,
  psi_src_barrier_sign_component_witness_v1,
  chi_src_local_barrier_sign_stability_witness_v1
)
```

with the declared properties:

1. actual source-side sign witness,
2. actual branch-to-packet lift,
3. observer-free in the witness domain,
4. below full source-topology nontriviality discharge,
5. below basis-independent selector promotion,
6. below quotient-safe `QW-2191` resolution,
7. below current selector closure,
8. no false pass.

## Hard limits

`F141` does not discharge:

1. full source-topology nontriviality,
2. basis-independent selector promotion,
3. quotient-safe `QW-2191` resolution,
4. current selector closure,
5. current global `QW-2191` discharge,
6. ToE closure.

## Recommended next move

The correct next move is:

1. test whether the current repo really exports this actual
   `barrier_protected_sign_class_v1` witness in a guardrail-consistent way,
2. then attempt an actual observer-free scope witness over
   `tau_src_candidate_v1`,
3. only after both attempt any full source-topology nontriviality lift.
