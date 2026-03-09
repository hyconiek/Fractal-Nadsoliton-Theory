# P00 Sandbox Strict-Core Ingredient Candidate Probe

Status: `P00_SANDBOX_STRICT_CORE_INGREDIENT_CANDIDATE_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Check exactly how far

```text
G_strict_core_selector_source_sandbox_candidate_v0
```

gets against the real `F29` admission contract.

## Probe matrix

### Clause 1: strict-core source export

Current sandbox verdict:

```text
PARTIAL
```

Reason:

1. the candidate is a genuinely new object symbol,
2. but it is exported only inside this removable sandbox folder,
3. therefore it is not yet an official strict-core repo export.

### Clause 2: internal orientation discharge

Current sandbox verdict:

```text
NO
```

Reason:

1. the candidate contains only
   `rho_int_orientation_request_slot_v0`,
2. that slot is a request placeholder,
3. no actual strict-core internal orientation datum is discharged.

### Clause 3: bridge discharge

Current sandbox verdict:

```text
NO
```

Reason:

1. the candidate contains only
   `beta_strict_selector_bridge_request_slot_v0`,
2. that slot names a missing bridge task,
3. no actual strict-core bridge is discharged.

### Clause 4: selector reduction discharge

Current sandbox verdict:

```text
NO
```

Reason:

1. no basis-covariant or target-independent selector reduction is exported
   from this candidate,
2. the candidate still sits strictly before reduction.

### Clause 5: downstream operator reachability

Current sandbox verdict:

```text
NO
```

Reason:

1. the candidate contains only
   `lambda_pair1_reachability_request_slot_v0`,
2. no actual strict-core route to `A_1(pair1)` is exported.

### Clause 6: no silent ontological substitution

Current sandbox verdict:

```text
YES
```

Reason:

1. the candidate does not identify `K_legacy_ont` with `K_strict_gate`,
2. the candidate does not transfer legacy physical-role semantics,
3. the candidate does not treat extension-only acceptance as strict-core
   derivation.

## Net result

The sandbox candidate currently satisfies:

1. one safety clause fully,
2. one source-export clause only partially at sandbox level,
3. four positive strict-core construction clauses not at all.

## Interpretation

This is still useful because it separates two questions that were previously
mixed together:

1. whether the project needs a new strict-core ingredient route at all,
2. which exact positive strict-core subobjects are still missing.

The current answer is:

```text
the route is non-empty,
but almost all positive strict-core work remains genuinely undone
```

## Hard limits

This probe does not claim:

1. a new official strict-core export,
2. admissible `S_sel_int`,
3. actual `E_orient`,
4. actual strict-core selector closure,
5. actual ToE closure.
