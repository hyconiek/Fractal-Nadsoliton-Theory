# P124 First Source-Seed Genuinely-New Strict-Core Object Clause Probe

Status: `P124_EXECUTABLE_FIRST_SOURCE_SEED_GENUINELY_NEW_STRICT_CORE_OBJECT_CLAUSE_PROBE_READY`
As of: `2026-03-08`

## Goal

Test whether the first clause-level admissibility question is already forced to:

```text
does S_sel_int_candidate_seed_v0 count as a genuinely new strict-core source object?
```

## Probe rule

The probe may succeed only if all of the following are true:

1. `N134` already reduces the next constructive move to one first
   admissibility-upgrade target,
2. the target still reuses only current route ingredients:
   `QW-2206_local_topological_protection_layer` and `sigma_int_candidate`,
3. `F34` still requires `genuinely_new_strict_core_source_object_required`,
4. no earlier clause has priority over that clause,
5. the result is stated only as a first clause-level test target, not as a
   satisfied admissibility verdict.

## Allowed conclusion

If the probe passes, the only allowed conclusion is:

```text
the first clause-by-clause admissibility test is reduced to the
genuinely-new strict-core source-object requirement
```

No stronger conclusion is allowed.
