# P05 Sandbox Populated Basis-Pair Instance Layer Artifact Schema Probe

Status: `P05_SANDBOX_POPULATED_BASIS_PAIR_INSTANCE_LAYER_ARTIFACT_SCHEMA_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Check whether the direct layer attack is honestly stronger than `F04` while
remaining below any actual populated-instance export.

## What is checked

`P05` checks whether the new layer attack:

1. writes one packet-ready artifact schema for the missing populated-instance
   layer,
2. stays strict-core-only,
3. avoids importing the axiom-augmented branch,
4. avoids claiming actual populated `u_1,u_2`,
5. avoids claiming actual strict-core `theta_1,theta_2`,
6. remains below actual strict-core instance export.

## Result matrix

### Packet-ready artifact schema for the missing populated-instance layer

Current verdict after `F05`:

```text
YES
```

Reason:

1. `F05` explicitly assembles the object, inputs, population rule, slice role,
   conditional serialization rule, current absences, and forbidden claims,
2. therefore the missing layer is now described as one packet-ready artifact
   schema rather than as a generic absence only.

### Import of axiom-augmented branch

Current verdict after `F05`:

```text
NO
```

Reason:

1. the schema uses only strict-core-side formulas and current strict-side
   absences,
2. it does not use the non-strict theta branch.

### Actual populated strict-core basis-pair instance

Current verdict after `F05`:

```text
NO
```

Reason:

1. `C49` gives only a conditional schema,
2. `C50` still denies the upstream strict-core source skeleton,
3. `F05` preserves both absences.

### Actual strict-core `theta_1`, `theta_2` export

Current verdict after `F05`:

```text
NO
```

Reason:

1. the direct layer attack does not solve the upstream provider problem,
2. the schema records that absence instead of hiding it.

### Actual strict-core basis-pair export

Current verdict after `F05`:

```text
NO
```

Reason:

1. an artifact schema is not an emitted instance,
2. `persisted_instance_absent` remains explicit in the packet.

## Hard limits

`P05` does not establish:

1. actual `theta_1`, `theta_2`,
2. actual populated `u_1`, `u_2`,
3. actual internal orientation datum,
4. actual `E_orient`,
5. admissible `S_sel_int`,
6. actual strict-core selector closure,
7. actual ToE closure.
