# N06 Sandbox Strict-Core Theta/Input Provider Artifact Schema Status Note

Status: `N06_SANDBOX_STRICT_CORE_THETA_INPUT_PROVIDER_ARTIFACT_SCHEMA_STATUS_NOTE_NO_FALSE_PASS`
As of: `2026-03-09`

## Strongest honest conclusion

The sandbox attack succeeds in one narrow but useful sense:

```text
the missing upstream strict-core theta/input provider is now exposed as one
packet-ready artifact schema bound to the F05 downstream instance contract
```

This is stronger than `N05` because the sandbox no longer stops at the
missing populated-instance layer itself.
It now also writes the provider that would have to feed that layer as one
explicit artifact schema with:

1. provider object,
2. target outputs,
3. downstream consumer contract,
4. support lane,
5. strict absence,
6. fallback contrast,
7. forbidden claims.

## What did not happen

This attack did **not** derive:

1. actual `theta_1`, `theta_2`,
2. actual populated `u_1`, `u_2`,
3. actual basis-pair export,
4. actual internal orientation datum,
5. actual `E_orient`,
6. admissible `S_sel_int`,
7. strict-core selector closure.

## Why this attack is still useful

It sharpens the frontier into two clean residual pieces:

1. the missing provider layer is no longer semantically vague, because its
   artifact schema is explicit and tied to `F05`,
2. the remaining live gap is now even narrower:
   strict core still lacks any emitted provider instance that could satisfy
   that schema without importing the axiom-augmented branch.

## Honest next move

If the sandbox is continued, the next clean move is:

1. either write one incompatibility boundary showing that no emitted
   strict-core provider instance can be obtained on current inputs,
2. or attack a concrete candidate field inside the provider schema that could
   narrow the missing instance problem further.
