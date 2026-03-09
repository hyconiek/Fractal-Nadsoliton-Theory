# T02 Sandbox Theta-Like Input Source Sketch Attack Scope

Status: `T02_SANDBOX_THETA_LIKE_INPUT_SOURCE_SKETCH_ATTACK_SCOPE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Attack the next weakest part of:

```text
rho_int_orientation_request_slot_v1
```

by asking:

```text
can the slot now acquire one repo-consistent theta-like input source sketch
without pretending that strict core already exports actual theta_1, theta_2?
```

## Repo-consistency constraints

This attack must stay consistent with:

1. `R1`
   - the residual-orientation target slot explicitly requires inputs
     `theta_1`, `theta_2`,
2. `C47`
   - a class-level candidate orientation slice already exists,
3. `C48`
   - a minimal basis-pair export skeleton is packet-ready,
4. `C35`
   - strict core still does not export actual `theta_1`, `theta_2`, and only
     an axiom-augmented source branch is currently available,
5. `B2/N7`
   - no internal orientation datum is currently derived in strict core.

## Narrow success condition

This attack counts as useful only if it upgrades the rho slot from:

```text
target-slot-aligned orientation request scaffold
```

to:

```text
target-slot-aligned orientation request scaffold
with one theta-like input source sketch
```

while still keeping explicit:

1. no actual strict-core `theta_1`, `theta_2`,
2. no actual populated basis pair `u_1`, `u_2`,
3. no actual internal orientation datum,
4. no strict-core selector closure.

## Failure condition

The attack fails if it drifts into:

1. treating class-level formulas as actual exported phases,
2. treating the axiom-augmented `theta^* = 0` branch as strict-core export,
3. treating the basis-pair export skeleton as a populated basis pair,
4. treating the source sketch as actual `E_orient`.
