# P02 Sandbox Theta-Like Input Source Sketch Attack Probe

Status: `P02_SANDBOX_THETA_LIKE_INPUT_SOURCE_SKETCH_ATTACK_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Check whether the theta-like input attack really improves the rho slot in a
repo-consistent way.

## What is checked

`P02` checks whether the refined slot:

1. is still below actual `theta_1`, `theta_2` export,
2. is now aligned to the real required-input structure from `R1`,
3. now carries the class-level orientation-slice candidate from `C47`,
4. now carries the basis-pair export skeleton from `C48`,
5. keeps explicit the `C35` split:
   strict-core theta absent,
   axiom-augmented theta branch present but non-strict.

## Result matrix

### Theta-like inputs as actual strict-core export

Current verdict after `F02`:

```text
NO
```

Reason:

1. `C35` still says strict core does not export actual `theta_1`, `theta_2`,
2. the current sketch does not populate them.

### Theta-like inputs as source sketch

Current verdict after `F02`:

```text
YES
```

Reason:

1. `R1` fixes that two such inputs are required,
2. `C47` fixes the class-level orientation-slice role,
3. `C48` fixes the minimal basis-pair export skeleton,
4. `C35` fixes the exact present absence/presence split.

### Rho slot as theta-source-sketch-aware request scaffold

Current verdict after `F02`:

```text
YES
```

Reason:

1. the slot is now no longer only target-slot-bound,
2. it is now also parameterized by the real theta-like input frontier.

### Clause-2 impact on the `F29` admission contract

Current sandbox verdict after `F02`:

```text
PARTIAL_UPGRADE_ONLY
```

Reason:

1. the orientation axis still does not discharge internal orientation,
2. but it is now written in the exact input language needed by the real
   residual-orientation target slot.

## Hard limits

`P02` does not establish:

1. actual `theta_1`, `theta_2`,
2. actual populated `u_1`, `u_2`,
3. actual internal orientation datum,
4. actual `E_orient`,
5. admissible `S_sel_int`,
6. actual strict-core selector closure,
7. actual ToE closure.
