# P01 Sandbox Rho-Int Orientation Slot Attack Probe

Status: `P01_SANDBOX_RHO_INT_ORIENTATION_SLOT_ATTACK_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Check whether the rho-slot attack really improves the sandbox candidate in a
repo-consistent way.

## What is checked

`P01` checks whether the refined slot:

1. is still below actual internal orientation discharge,
2. is now anchored to the real residual-orientation target-slot language from
   `R1`,
3. is now constrained by the real future `E_orient` contract from `F32`,
4. keeps explicit the strict negative boundaries from `B2`, `N7`, and
   `P135`,
5. avoids laundering extension/axiom acceptance into strict-core derivation.

## Result matrix

### Rho slot as generic placeholder

Current verdict after `F01`:

```text
NO
```

Reason:

1. the slot is no longer free-floating,
2. it is now tied to a real target-slot and contract layer.

### Rho slot as target-slot-aligned request scaffold

Current verdict after `F01`:

```text
YES
```

Reason:

1. `R1` provides a packet-ready residual-orientation target slot,
2. `F32` provides the admissibility contract for any future orientation
   export,
3. the refined slot explicitly binds itself to both.

### Rho slot as actual internal orientation datum

Current verdict after `F01`:

```text
NO
```

Reason:

1. `B2` still says no strict-core internal orientation datum is found,
2. `N7` still says the current `sigma_int` route does not derive the residual
   orientation datum,
3. `P135` still says the orientation-export branch is not discharged.

### Clause-2 impact on the `F29` admission contract

Current sandbox verdict after `F01`:

```text
PARTIAL_UPGRADE_ONLY
```

Reason:

1. before the attack, the orientation axis was only a generic request slot,
2. after the attack, that axis becomes a target-slot-aligned request
   scaffold,
3. but it still does not discharge internal orientation.

## Hard limits

`P01` does not establish:

1. actual internal orientation datum,
2. actual `E_orient`,
3. actual strict-core bridge discharge,
4. admissible `S_sel_int`,
5. actual strict-core selector closure,
6. actual ToE closure.
