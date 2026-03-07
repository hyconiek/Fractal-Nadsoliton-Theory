# P67 Strict-Side Role-Equivalence Probe For Legacy Weinberg Role

Status: `P67_EXECUTED_STRICT_SIDE_ROLE_EQUIVALENCE_PROBE_FOR_LEGACY_WEINBERG_ROLE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F8`, the next honest direct question is:

```text
does the current repo already export an explicit legacy-to-strict
role-equivalence verdict for the legacy Weinberg-angle role?
```

## Result

The route is mixed but still negative at the theorem-relevant level:

```text
CURRENT_REPO_EXPORTS_STRICT_SIDE_SIN2_THETA_W_MZ_CANDIDATE_BUT_NO_EXPLICIT_LEGACY_WEINBERG_ROLE_EQUIVALENCE_VERDICT_AFTER_P67
```

## What was checked

`P67` checks two different things and keeps them separate:

1. whether the strict side exports a real candidate object carrying the modern
   Weinberg-angle observable label,
2. whether the repo also exports an explicit semantic-transfer verdict saying
   that this object is the retained role-equivalent successor of the old
   legacy Weinberg formula.

## Why it only partially succeeds

On the current repo state:

1. `QW-2068`, `QW-2069`, `QW-2098`, and Release `4.9` together do export a
   real strict-side candidate object `sin2_theta_w_mz`,
2. but none of the current strict-side sources exports an explicit verdict
   that `sin2_theta_w_mz` is the retained role-equivalent successor of the
   legacy `alpha_geo/12` Weinberg-angle semantics,
3. therefore the retained-side role-equivalence branch is narrowed but not yet
   discharged.

## Real reduction after `P67`

The retained frontier for the legacy Weinberg role is no longer:

```text
generic role-equivalence retention blocker
```

It is now:

```text
one explicit semantic-transfer blocker attached to sin2_theta_w_mz
```

namely:

```text
explicit_legacy_to_strict_semantic_transfer_verdict_identifying_sin2_theta_w_mz_as_the_retained_strict_side_successor_of_the_legacy_weinberg_angle_role
```

## What P67 does not claim

`P67` does not claim:

- that the retained branch is already discharged,
- that the replaced branch is already discharged,
- that `sin2_theta_w_mz` automatically inherits all legacy semantics,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. formalize theorem-level that the current repo exports a strict-side
   Weinberg candidate object but no explicit legacy-to-strict role-equivalence
   verdict,
2. keep the replaced branch separate,
3. do not silently treat `sin2_theta_w_mz` as if it already carried the full
   legacy `alpha_geo/12` semantics.
