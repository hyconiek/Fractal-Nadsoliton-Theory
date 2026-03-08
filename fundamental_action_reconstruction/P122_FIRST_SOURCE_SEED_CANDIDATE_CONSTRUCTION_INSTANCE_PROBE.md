# P122 First Source-Seed Candidate Construction Instance Probe

Status: `P122_EXECUTABLE_FIRST_SOURCE_SEED_CANDIDATE_CONSTRUCTION_INSTANCE_PROBE_READY`
As of: `2026-03-07`

## Goal

Test whether the current repo now reduces the next constructive move to one
explicit first candidate construction instance:

```text
S_sel_int_candidate_seed_v0
```

and whether this instance is still kept below the admissible source-object
threshold.

## Probe rule

The probe may succeed only if all of the following are true:

1. `N132` already reduces the next constructive move to one precursor route,
2. `F36` packages exactly one candidate construction instance on that route,
3. the packaged instance is anchored only on
   `QW-2206_local_topological_protection_layer` and `sigma_int_candidate`,
4. the instance is marked `candidate_seed_instance_only`,
5. the instance is explicitly blocked from counting as admissible `S_sel_int`,
   `E_orient`, or downstream closure.

## Allowed conclusion

If the probe passes, the only allowed conclusion is:

```text
the next constructive move is reduced to one first candidate source-seed
construction instance
```

No stronger conclusion is allowed.
