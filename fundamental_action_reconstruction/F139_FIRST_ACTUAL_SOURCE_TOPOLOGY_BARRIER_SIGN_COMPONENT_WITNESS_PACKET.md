# F139 First Actual Source Topology Barrier Sign Component Witness Packet

Status: `F139_EXECUTED_FIRST_ACTUAL_SOURCE_TOPOLOGY_BARRIER_SIGN_COMPONENT_WITNESS_PACKET_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

After `F138/P226/N246`, the next honest move is still not:

1. an actual full barrier-protected sign discharge,
2. an actual full source-topology nontriviality discharge,
3. a basis-independent selector promotion discharge,
4. a quotient-safe `QW-2191` resolution,
5. a current selector closure proof.

It is narrower:

```text
export one actual source-side scalar sign component witness
for the current source-limit kernel core
together with one explicit positive barrier margin on the declared core branch
strictly below full barrier-protected sign discharge
and strictly below QW-2191 discharge
```

`F139` executes exactly that move.

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
K_strict(0) = cos(phi) = 0.9868259031903286
```

Reuse the actual nonzero-flow component support from `F138`:

```text
xi_src_nonzero_flow_component_witness_v1 := |cos(phi)| = 0.9868259031903286 > 0
```

## Actual component witness

Freeze one explicit scalar barrier margin on the current core branch:

```text
delta_src_barrier_sign_margin_v1 := (pi/2) - |phi|
```

Numerically:

```text
delta_src_barrier_sign_margin_v1 = 1.4082963267948965 > 0
```

Freeze one actual source-side scalar sign component witness:

```text
psi_src_barrier_sign_component_witness_v1 := sign(cos(phi))
```

Numerically on the current repo state:

```text
psi_src_barrier_sign_component_witness_v1 = +1
```

## Meaning of the witness

This witness is intended only as:

1. one actual scalar barrier-margin witness extracted from the current
   exported source-limit core datum,
2. one actual scalar sign component witness for the source-limit kernel core,
3. a source-side barrier-sign component witness on the declared core branch,
4. a current positive witness strictly before full barrier-protected sign
   discharge,
5. a current positive witness strictly before full source-topology
   nontriviality, basis-independent selector promotion, and quotient-safe
   `QW-2191` resolution.

It is not yet:

1. a full barrier-protected sign discharge,
2. a full source-topology nontriviality discharge,
3. a basis-independent selector witness,
4. a quotient-safe `QW-2191` witness,
5. a current selector closure proof.

## Why this is still kernel-split-safe

`F139` remains kernel-split-safe because:

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
3. no observer-side algebra is used as proof of source barrier sign.

This preserves the declared order:

```text
nadsoliton -> light -> matter -> emergent observer
```

## Result

`F139` exports one explicit scalar barrier margin and one actual source-side
scalar sign witness:

```text
delta_src_barrier_sign_margin_v1 := (pi/2) - |phi| = 1.4082963267948965 > 0
psi_src_barrier_sign_component_witness_v1 := sign(cos(phi)) = +1
```

with the declared properties:

1. actual scalar barrier margin,
2. actual scalar sign component witness,
3. source-side only,
4. observer-free in the witness domain,
5. below full barrier-protected sign discharge,
6. below full source-topology nontriviality discharge,
7. below basis-independent selector promotion,
8. below quotient-safe `QW-2191` resolution,
9. below current selector closure,
10. no false pass.

## Hard limits

`F139` does not discharge:

1. full barrier-protected sign of `tau_src_candidate_v1`,
2. full source-topology nontriviality,
3. basis-independent selector promotion,
4. quotient-safe `QW-2191` resolution,
5. current selector closure,
6. current global `QW-2191` discharge,
7. ToE closure.

## Recommended next move

The correct next move is:

1. test whether the current repo really exports this actual scalar barrier-sign
   component witness in a guardrail-consistent way,
2. keep the result strictly below full barrier-protected sign discharge,
3. avoid any promotion to full source-topology nontriviality, selector
   promotion, or `QW-2191` discharge.
