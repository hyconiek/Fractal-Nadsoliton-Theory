# F127 First Source Topology Invariant Candidate Packet

Status: `F127_EXECUTED_FIRST_SOURCE_TOPOLOGY_INVARIANT_CANDIDATE_PACKET_FUTURE_ROUTE_ONLY_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

After `T14`, the most honest next move is not a global promotion claim.

It is narrower:

```text
export one explicit future-only tau_src candidate packet at the source/kernel limit
without claiming basis-independent selector promotion
and without claiming current QW-2191 discharge
```

`F127` executes exactly that move.

## Fixed future-only packet

`F127` defines one explicit future-only candidate packet:

```text
tau_src_candidate_v1 :=
(
  source_limit_tag_v1,
  phi_barrier_tag_v1,
  T_flow^(0)
)
```

with:

```text
source_limit_tag_v1 := d -> 0
phi_barrier_tag_v1 := fixed_nonzero_phi_kernel_core_barrier
T_flow^(0) := cos(phi) * e_topo
```

## Source-side meaning

`tau_src_candidate_v1` is intended only as:

1. a candidate non-trivial source-topology packet,
2. a packet upstream of observer stages,
3. a future route component for `T14`.

It is not yet:

1. a discharged topological invariant theorem,
2. a selector datum,
3. a basis-independent promotion witness,
4. a quotient-safe `QW-2191` resolution,
5. a current selector closure proof.

## Strict-kernel role

The packet uses the strict kernel only through the already exported core-limit
control datum from `F74`:

```text
K_strict(0) = cos(phi)
```

This packet does **not**:

1. identify `K_strict` with `K_legacy_ont`,
2. claim a legacy-to-strict bridge,
3. transfer legacy physical-role semantics.

## Observer role

Observer asymmetry remains outside the packet.

Inside `F127`:

1. observer is downstream only,
2. observer is not a selector source,
3. observer may serve later only as pushforward witness if a real upstream
   selector source is exported.

## Result

`F127` exports one explicit future-only source-topology candidate packet:

```text
tau_src_candidate_v1
```

with the declared properties:

1. source-side only,
2. strict-core future-route only,
3. no external selector import,
4. observer downstream only,
5. kernel-split-safe,
6. no false pass.

## Hard limits

`F127` does not discharge:

1. existence of a non-trivial topological invariant,
2. basis-independent promotion,
3. quotient-safe selector promotion at `QW-2191`,
4. current selector closure,
5. current global `QW-2191` discharge,
6. ToE closure.

## Recommended next move

The correct next move is:

1. test whether the current repo really exports this candidate packet in a
   guardrail-consistent way,
2. keep the result below basis-independence and below quotient safety,
3. avoid any observer-based global promotion claim.
