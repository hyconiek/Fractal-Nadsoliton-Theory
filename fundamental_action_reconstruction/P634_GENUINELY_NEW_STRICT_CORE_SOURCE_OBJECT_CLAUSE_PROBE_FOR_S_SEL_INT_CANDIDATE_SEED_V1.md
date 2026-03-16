# P634 Genuinely‑New Strict‑Core Source Object Clause Probe For `S_sel_int_candidate_seed_v1`

Status: `P634_EXECUTABLE_FIRST_CLAUSE_PROBE_READY_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

Test the first admissibility clause directly, on the strict‑sigma‑int upgraded seed:

```text
does S_sel_int_candidate_seed_v1 already count as a genuinely new strict‑core
source object (for admissible S_sel_int)?
```

## Probe rule

The probe may return a positive result only if all of the following are true:

1. `S_sel_int_candidate_seed_v1` is exported as a genuinely new strict‑core
   source object rather than only as a seed candidate instance,
2. the object is not merely a naming wrapper around the current seed
   ingredients (even if one ingredient is now a strict source‑upgrade datum),
3. the result is supported inside strict core (no axiom/selector import),
4. the conclusion is stated only at the first‑clause level (no admissible
   `S_sel_int`, no selector closure, no `QW‑2191` discharge).

## Allowed conclusion

If the probe fails, the only allowed conclusion is:

```text
current repo does not yet show that S_sel_int_candidate_seed_v1 satisfies the
genuinely‑new strict‑core source‑object clause
```

