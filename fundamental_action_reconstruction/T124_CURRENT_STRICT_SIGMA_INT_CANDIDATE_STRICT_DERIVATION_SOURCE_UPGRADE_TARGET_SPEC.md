# T124 Current Strict Sigma-Int Candidate Strict-Derivation/Source-Upgrade Target Spec

Status: `T124_CURRENT_STRICT_SIGMA_INT_CANDIDATE_STRICT_DERIVATION_SOURCE_UPGRADE_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

Before `F307/N418`, the strict-core sigma-int bridge lane kept one explicit
upstream missing ingredient:

```text
a strict derivation or source-object upgrade for sigma_int_candidate
```

On the current repo state (`P390/P391`), the strict sigma-int lane exports a
strict-side FR-sign source-upgrade package and the upgraded strict datum:

```text
sigma_int_strict_derived_v1 := chi_FR_strict_v1(gamma_pi1_v1) ∈ {+1,-1}
```

via `F307/N418` (explicit strict-side premise; no hybrid reuse).

Therefore `T124` is no longer a “current missing-object target spec” for the
strict sigma-int lane. It is kept as:

1. a historical target-spec / acceptance-test record for sigma-int
   strict-derivation/source-upgrade,
2. a guardrail against silently treating candidate-only sigma-int as strict,
3. a guardrail against promoting axiom-lane-only material into strict core.

This still matters as an integrity constraint because any strict-core
bridge/export-map use of a sigma-int datum must keep its strict provenance
explicit.

Historically, this missing ingredient was treated as a prerequisite for any
honest strict-core equivalence/export map of the form:

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
