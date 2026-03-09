# P06 Sandbox Strict-Core Theta/Input Provider Artifact Schema Probe

Status: `P06_SANDBOX_STRICT_CORE_THETA_INPUT_PROVIDER_ARTIFACT_SCHEMA_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Check whether the upstream provider attack is honestly stronger than `F05`
while remaining below any actual strict-core theta export.

## What is checked

`P06` checks whether the new provider attack:

1. writes one packet-ready artifact schema for the missing upstream provider,
2. binds that provider schema explicitly to `F05`,
3. stays strict-core-only,
4. avoids identifying the axiom-augmented branch with strict core,
5. avoids claiming actual `theta_1`, `theta_2`,
6. remains below actual strict-core provider emission.

## Result matrix

### Packet-ready artifact schema for the missing upstream provider

Current verdict after `F06`:

```text
YES
```

Reason:

1. `F06` explicitly assembles provider object, outputs, downstream contract,
   support lane, strict absence, fallback contrast, and forbidden claims,
2. therefore the upstream provider layer is now described as one packet-ready
   artifact schema rather than as a generic residual absence only.

### Binding to downstream `F05` contract

Current verdict after `F06`:

```text
YES
```

Reason:

1. `F06` names `PopulatedBasisPairInstanceLayerArtifactSchema_v0` as the
   downstream consumer contract,
2. so the provider attack is now coupled to the already explicit missing
   instance layer.

### Identification of axiom-augmented fallback with strict core

Current verdict after `F06`:

```text
NO
```

Reason:

1. the packet records the fallback branch only as methodological contrast,
2. it explicitly refuses to treat that branch as a strict-core provider.

### Actual strict-core `theta_1`, `theta_2` export

Current verdict after `F06`:

```text
NO
```

Reason:

1. `C35/C50` still block that export,
2. `F06` records the missing provider instance instead of resolving it.

### Actual upstream provider emission

Current verdict after `F06`:

```text
NO
```

Reason:

1. an artifact schema is not an emitted provider,
2. `provider_instance_absent` remains explicit in the packet.

## Hard limits

`P06` does not establish:

1. actual `theta_1`, `theta_2`,
2. actual populated `u_1`, `u_2`,
3. actual internal orientation datum,
4. actual `E_orient`,
5. admissible `S_sel_int`,
6. actual strict-core selector closure,
7. actual ToE closure.
