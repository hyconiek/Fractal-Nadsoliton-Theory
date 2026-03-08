# P123 First Source-Seed Admissibility Upgrade Target Probe

Status: `P123_EXECUTABLE_FIRST_SOURCE_SEED_ADMISSIBILITY_UPGRADE_TARGET_PROBE_READY`
As of: `2026-03-07`

## Goal

Test whether the current repo now reduces the next constructive move to one
explicit attempted admissibility-upgrade target:

```text
S_sel_int_candidate_seed_v0
against
minimal_admissible_S_sel_int_construction_contract
```

## Probe rule

The probe may succeed only if all of the following are true:

1. `N133` already reduces the next constructive move to one first candidate
   source-seed instance,
2. `F34` already exports one minimal admissibility contract for `S_sel_int`,
3. `F37` packages exactly the ordered pair:
   candidate seed instance + minimal admissibility contract,
4. the package is marked only as an attempted admissibility-upgrade target,
5. the package is explicitly blocked from counting as already admissible
   `S_sel_int`, `E_orient`, or downstream closure.

## Allowed conclusion

If the probe passes, the only allowed conclusion is:

```text
the next constructive move is reduced to one first source-seed admissibility-upgrade target
```

No stronger conclusion is allowed.
