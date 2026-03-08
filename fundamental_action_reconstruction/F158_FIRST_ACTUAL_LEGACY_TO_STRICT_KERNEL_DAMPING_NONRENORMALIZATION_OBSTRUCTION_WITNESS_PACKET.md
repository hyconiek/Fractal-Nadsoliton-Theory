# F158 First Actual Legacy-to-Strict Kernel Damping Nonrenormalization Obstruction Witness Packet

Status: `F158_EXECUTED_FIRST_ACTUAL_LEGACY_TO_STRICT_KERNEL_DAMPING_NONRENORMALIZATION_OBSTRUCTION_WITNESS_PACKET_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `N267`, the amplitude layer of the `T16` route is already discharged.

The next honest move is now the second kernel-comparison layer:

```text
R_damp : (1 + beta_tors * d) -> (1 + beta * d^eta)
```

`F158` executes one actual damping non-renormalization obstruction attempt,
but still below phase and below the full strengthened nonbridge theorem.

## Scope

The declared scope is:

```text
kernel-comparison damping layer only
for the current legacy-to-strict comparison question
```

This packet does not yet touch:

1. phase/frequency non-conformality,
2. the full strengthened nonbridge theorem.

## Fixed support packet

Reuse:

1. `P47/N50`:
   the current repo exports no explicit `beta_tors -> (beta, eta)` translation
   rule,
2. `N117`:
   package-level legacy-to-strict nontransfer is discharged,
3. `N267`:
   the amplitude layer is already discharged, so damping is now the next
   active obstruction layer.

Freeze one actual support packet:

```text
Lambda_damp_nonbridge_full_support_v1 :=
(
  no_exported_explicit_beta_tors_to_beta_eta_translation_rule,
  legacy_to_strict_package_nontransfer_on_current_repo_state,
  A_abs_nonbridge_actual_obstruction_witness_v1
)
```

## Result

`F158` exports one actual damping-layer obstruction witness:

```text
R_damp_nonbridge_actual_obstruction_witness_v1 :
  (K_legacy_ont, K_strict_gate)
    -> R_damp_nonbridge_obstruction_target_v1
```

meaning only:

1. at the declared damping comparison layer, the current repo exports no
   admissible explicit renormalization flow from `beta_tors` to `(beta, eta)`,
2. this absence is now packaged as one actual damping obstruction witness,
3. the result still remains below phase and below the full strengthened
   nonbridge theorem.

## Hard limits

`F158` does not discharge:

1. actual phase/frequency non-conformal obstruction,
2. actual strengthened nonbridge theorem,
3. actual legacy-to-strict bridge derivation,
4. current branch selection between bridge and nonbridge,
5. strict-core selector closure,
6. global selector closure,
7. global `QW-2191` discharge,
8. ToE closure.

## Recommended next move

The correct next move is:

1. test whether the current repo really exports this full damping-layer
   obstruction witness,
2. if it does, move next to phase rather than overclaiming strengthened
   nonbridge,
3. if it fails, package the exact missing damping ingredient.
