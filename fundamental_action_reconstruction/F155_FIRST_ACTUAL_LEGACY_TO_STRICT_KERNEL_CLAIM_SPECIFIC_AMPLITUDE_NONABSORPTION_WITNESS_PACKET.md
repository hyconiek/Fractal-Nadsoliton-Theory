# F155 First Actual Legacy-to-Strict Kernel Claim-Specific Amplitude Nonabsorption Witness Packet

Status: `F155_EXECUTED_FIRST_ACTUAL_LEGACY_TO_STRICT_KERNEL_CLAIM_SPECIFIC_AMPLITUDE_NONABSORPTION_WITNESS_PACKET_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

After `N264`, the `T16` route now contains one first actual component-level
amplitude nonabsorption witness in the legacy Weinberg amplitude scope.

The next honest move is still not to claim full
`A_abs_nonbridge_obstruction_target_v1`.

It is narrower:

```text
freeze one actual claim-specific amplitude nonabsorption witness
above the first component witness
but below full amplitude obstruction
```

`F155` executes exactly that move.

## Scope

The declared scope remains:

```text
legacy Weinberg-angle amplitude role
sin^2(theta_W)=alpha_geo/12
```

This packet does not extend beyond that claim-specific amplitude scope.

## Fixed support packet

Reuse `N264`:

```text
A_abs_nonbridge_component_witness_v1 :
  (K_legacy_ont, K_strict_gate)
    -> amplitude_nonabsorption_component_obstruction_tag_v1
```

Freeze one actual support packet:

```text
W_abs_nonbridge_claim_specific_support_v1 :=
(
  A_abs_nonbridge_component_witness_v1,
  legacy_alpha_geo_weinberg_role_marker,
  strict_side_candidate_sin2_theta_w_mz_present,
  no_explicit_role_equivalence_verdict,
  bridge_nonbridge_frontier_still_undecided
)
```

## Result

`F155` exports one actual claim-specific amplitude nonabsorption witness:

```text
A_abs_nonbridge_claim_specific_actual_witness_v1 :
  (K_legacy_ont, K_strict_gate)
    -> claim_specific_amplitude_nonabsorption_obstruction_target_v1
```

meaning only:

1. the legacy Weinberg amplitude role is now blocked at one actual
   claim-specific obstruction target,
2. the witness is stronger than the first component witness,
3. it remains below full amplitude nonabsorption obstruction.

## Hard limits

`F155` does not discharge:

1. full `A_abs_nonbridge_obstruction_target_v1`,
2. actual damping non-renormalization obstruction,
3. actual phase/frequency non-conformal obstruction,
4. actual strengthened nonbridge theorem,
5. actual legacy-to-strict bridge derivation,
6. current branch selection between bridge and nonbridge,
7. strict-core selector closure,
8. global selector closure,
9. global `QW-2191` discharge,
10. ToE closure.

## Recommended next move

The correct next move is:

1. test whether the current repo really exports this claim-specific amplitude
   witness,
2. only after that try one honest lift to full
   `A_abs_nonbridge_obstruction_target_v1`,
3. keep damping and phase layers downstream of that decision.
