# F131 First Source Topology Barrier-Protected Sign Subtarget Packet

Status: `F131_EXECUTED_FIRST_SOURCE_TOPOLOGY_BARRIER_PROTECTED_SIGN_SUBTARGET_PACKET_FUTURE_ROUTE_ONLY_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

After `F130/P218/N238`, the next honest move is still not:

1. an actual barrier-protected sign discharge,
2. a full source-topology nontriviality discharge,
3. a selector-promotion discharge,
4. a `QW-2191` discharge.

It is narrower:

```text
freeze one explicit second subtarget
Psi_src_barrier_sign_target_v1 : tau_src_candidate_v1 -> barrier_protected_sign_class_v1
below actual barrier-protected sign discharge
and below full source-topology nontriviality
```

`F131` executes exactly that move.

## Fixed future-only subtarget

Reuse the explicit packet from `F127`:

```text
tau_src_candidate_v1 :=
(
  source_limit_tag_v1,
  phi_barrier_tag_v1,
  T_flow^(0)
)
```

Freeze one explicit second future-only subtarget:

```text
Psi_src_barrier_sign_target_v1 : tau_src_candidate_v1 -> barrier_protected_sign_class_v1
```

## Meaning of the subtarget

`barrier_protected_sign_class_v1` is intended only as:

1. an abstract future witness class for a source-side barrier-protected sign in
   the kernel-limit packet,
2. the second component below `Lambda_src_nontriv_target_v1`,
3. a source-side target strictly before selector promotion.

It is not yet:

1. a discharged barrier-protected sign witness,
2. a discharged nonzero-flow witness,
3. a full source-topology nontriviality witness,
4. a selector datum,
5. a basis-independent selector witness,
6. a current selector closure proof,
7. a current global `QW-2191` discharge.

## Relation to F129 and F130

`F131` sits strictly below `F129` and alongside `F130`:

```text
tau_src_candidate_v1
  -> Xi_src_nonzero_flow_target_v1
  -> future nonzero-flow discharge

tau_src_candidate_v1
  -> Psi_src_barrier_sign_target_v1
  -> future barrier-protected sign discharge

  -> Lambda_src_nontriv_target_v1
  -> future full nontriviality discharge
  -> Pi_sel_src_target_v1
  -> future selector promotion
```

So inside `F131`:

1. no actual sign discharge is claimed,
2. no actual nonzero-flow discharge is claimed,
3. no full source-topology nontriviality is claimed,
4. no selector promotion is claimed,
5. no observer-side algebra is used as proof.

## Observer role

Observer remains outside the subtarget:

1. observer is downstream only,
2. observer is not part of the sign witness domain,
3. observer may appear later only after a real upstream discharge.

This preserves the declared order:

```text
nadsoliton -> light -> matter -> emergent observer
```

## Kernel-split safety

`F131` remains kernel-split-safe:

1. `K_strict_gate` is used only through the already exported source-limit
   operational control packet,
2. no legacy physical-role transfer is claimed,
3. no silent `K_legacy_ont -> K_strict_gate` identification is introduced.

## Result

`F131` exports one explicit future-only source-topology barrier-protected sign
subtarget:

```text
Psi_src_barrier_sign_target_v1 : tau_src_candidate_v1 -> barrier_protected_sign_class_v1
```

with declared properties:

1. source-side only,
2. observer-free in the witness domain,
3. future sign target only,
4. below actual sign discharge,
5. below full source-topology nontriviality,
6. below selector promotion,
7. below quotient-safe `QW-2191` resolution,
8. no false pass.

## Hard limits

`F131` does not discharge:

1. actual barrier-protected sign of `tau_src_candidate_v1`,
2. actual nonzero-flow of `tau_src_candidate_v1`,
3. full source-topology nontriviality,
4. basis-independent selector promotion,
5. quotient-safe `QW-2191` resolution,
6. current selector closure,
7. current global `QW-2191` discharge,
8. ToE closure.

## Recommended next move

The correct next move is:

1. test whether the current repo really exports this barrier-protected sign
   subtarget in a guardrail-consistent way,
2. keep the result below actual sign discharge,
3. avoid any promotion to full nontriviality, selector closure, or
   `QW-2191` discharge.
