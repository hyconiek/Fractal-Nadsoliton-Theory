# T01 Sandbox Rho-Int Orientation Slot Attack Scope

Status: `T01_SANDBOX_RHO_INT_ORIENTATION_SLOT_ATTACK_SCOPE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Attack only:

```text
rho_int_orientation_request_slot_v0
```

in a way consistent with the current repo.

The sandbox question is:

```text
can the generic rho request slot be upgraded
to one target-slot-aligned strict-core orientation request scaffold
without pretending that the missing internal orientation datum is already derived?
```

## Repo-consistency constraints

This attack must stay consistent with all of the following:

1. `B2`
   - no strict-core internal orientation datum is currently derived,
2. `R1`
   - one packet-ready strict-core target slot exists for the residual
     orientation datum,
3. `N7`
   - the current strict-core `sigma_int` route does not derive a residual
     orientation datum,
4. `F32`
   - any future `E_orient` must satisfy a strict admissibility contract,
5. `P135`
   - the current repo still does not export an orientation-export branch
     discharge,
6. `N125/N278`
   - theory-level acceptance outside strict core may not be laundered into
     strict-core derivation.

## Narrow success condition

This attack counts as useful only if it upgrades the rho slot from:

```text
generic request placeholder
```

to:

```text
target-slot-aligned orientation request scaffold
```

while still keeping explicit:

1. no actual internal orientation datum,
2. no actual `E_orient`,
3. no strict-core bridge discharge,
4. no strict-core selector closure.

## Failure condition

The attack fails if it drifts into any of the following:

1. treating `R1` target-slot export as already populated,
2. treating `F32` contract as already discharged,
3. treating axiom/extension acceptance as strict-core derivation,
4. implicitly using source-topology witness as hidden selector source.
