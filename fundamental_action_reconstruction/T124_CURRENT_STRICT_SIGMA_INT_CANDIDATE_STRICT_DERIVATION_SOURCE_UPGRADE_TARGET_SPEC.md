# T124 Current Strict Sigma-Int Candidate Strict-Derivation/Source-Upgrade Target Spec

Status: `T124_CURRENT_STRICT_SIGMA_INT_CANDIDATE_STRICT_DERIVATION_SOURCE_UPGRADE_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

`P4/P5` keep one upstream missing ingredient explicit on the strict-core
sigma-int route:

```text
a strict derivation or source-object upgrade for sigma_int_candidate
```

This missing ingredient is a prerequisite for any honest strict-core
equivalence/export map of the form:

```text
sigma_int_candidate -> residual orientation datum
```

because a candidate-only datum cannot be silently treated as a strict-derived
source on the selector-facing bridge lane (`B8/N7/N8`).

`T124` does **not** claim a strict derivation.

`T124` does something narrower:

- name the missing strict-derivation/source-upgrade ingredient as one explicit
  **future-only target object** with explicit acceptance tests,
- without claiming that the target is already discharged.

## Scope

`T124` is scoped only to the strict-core source status of:

```text
sigma_int_candidate
```

It does not decide:

1. theorem-level gauge-quotient safety (handled separately by `T123`),
2. discharge of `T2`,
3. export of the strict-core bridge/export map object from `N301`,
4. selector closure or `QW-2191` discharge,
5. ToE closure.

## Target object

If the answer is negative with respect to actual discharge but positive with
respect to sharp naming, export one future-only target object:

```text
Delta_sigma_int_candidate_strict_derivation_source_upgrade_target_v1
```

with the intended meaning:

```text
the current repo now names one exact future-only target object for the missing
strict derivation / source-object upgrade required before sigma_int_candidate
may be used as a strict-core bridge/export-map source datum.
```

## Acceptance tests (what would count as discharge)

An **actual** discharge of this target must at minimum provide:

1. an explicit strict-core derivation chain (or strict source upgrade) that
   upgrades sigma-int from candidate-only status to a strict-derived source
   datum on a declared domain,
2. no use of axiom-lane-only promotion (`AX3`) counted as strict-core
   derivation,
3. explicit statement of the resulting status change (candidate -> strict
   derived) as an exported object/packet/witness,
4. noncyclic input contract: no use of `theta_1,theta_2` and no use of a
   populated basis-pair instance as input (respects `N18`),
5. explicit separation from gauge-quotient safety:
   - either discharge `T123` as well, or keep gauge quotient safety explicitly
     marked as still open (no silent pass).

## Hard limits

`T124` must not claim:

1. strict derivation/source upgrade already discharged,
2. theorem-level gauge-quotient safety discharged,
3. strict-core bridge/export-map object export,
4. discharge of `T2`,
5. admissible `S_sel_int`,
6. strict-core selector closure,
7. `QW-2191` discharge,
8. ToE closure.

