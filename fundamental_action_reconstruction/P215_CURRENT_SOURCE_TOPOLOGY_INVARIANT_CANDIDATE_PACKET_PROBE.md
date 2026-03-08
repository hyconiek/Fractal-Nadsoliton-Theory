# P215 Current Source Topology Invariant Candidate Packet Probe

Status target:
`CURRENT_REPO_EXPORTS_ONE_FUTURE_SOURCE_TOPOLOGY_INVARIANT_CANDIDATE_PACKET_BELOW_BASIS_INDEPENDENCE_AND_QW2191_QUOTIENT_SAFETY_AFTER_P215`

## Probe question

Does the current repo already export one explicit future-only source-topology
candidate packet suitable for the `T14` route, while honestly remaining below:

1. basis-independent selector promotion,
2. quotient-safe `QW-2191` resolution,
3. current selector closure,
4. current global `QW-2191` discharge?

## Fixed inputs

1. `F74`
2. `N163`
3. `N234`
4. `F127`

## Acceptance conditions

The probe passes only if all of the following hold:

1. a future-only `tau_src_candidate_v1` packet is exported;
2. the packet is upstream of observer stages;
3. the packet uses no external selector import;
4. observer remains downstream only;
5. basis-independent selector promotion is not falsely claimed;
6. quotient-safe `QW-2191` promotion is not falsely claimed;
7. no global selector closure is claimed;
8. no global `QW-2191` discharge is claimed.

## Result shape

The only honest positive result is:

```text
current repo exports one future-only tau_src candidate packet
below basis-independence and below QW-2191 quotient safety
```

Any stronger result is a false pass.
