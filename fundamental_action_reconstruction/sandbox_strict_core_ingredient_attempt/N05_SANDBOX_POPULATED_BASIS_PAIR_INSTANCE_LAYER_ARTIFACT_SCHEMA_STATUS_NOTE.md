# N05 Sandbox Populated Basis-Pair Instance Layer Artifact Schema Status Note

Status: `N05_SANDBOX_POPULATED_BASIS_PAIR_INSTANCE_LAYER_ARTIFACT_SCHEMA_STATUS_NOTE_NO_FALSE_PASS`
As of: `2026-03-09`

## Strongest honest conclusion

The sandbox attack succeeds in one narrow but useful sense:

```text
the missing populated basis-pair instance layer is now exposed as one
packet-ready artifact schema rather than only as a residual absence
```

This is stronger than `N04` because the sandbox no longer stops at a
conditional theta-source rule candidate.
It now also writes the missing instance layer itself as one explicit artifact
schema with:

1. object,
2. required inputs,
3. population rule,
4. orientation-slice role,
5. conditional serialization rule,
6. current absences,
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

1. the missing populated basis-pair instance layer is no longer semantically
   vague, because its artifact schema is explicit,
2. the remaining live gap is now even narrower:
   strict core still lacks the upstream input provider that would actually
   populate that schema.

## Honest next move

If the sandbox is continued, the next clean move is:

1. either attack the upstream strict-core theta/input provider with this
   sharper instance-layer contract in hand,
2. or write one incompatibility boundary stating that the populated-instance
   layer cannot become emitted on current strict-core inputs.
