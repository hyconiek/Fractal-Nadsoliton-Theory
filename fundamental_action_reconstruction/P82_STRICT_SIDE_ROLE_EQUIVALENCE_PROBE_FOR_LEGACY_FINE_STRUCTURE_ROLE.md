# P82 Strict-Side Role-Equivalence Probe For Legacy Fine-Structure Role

Status: `P82_EXECUTED_STRICT_SIDE_ROLE_EQUIVALENCE_PROBE_FOR_LEGACY_FINE_STRUCTURE_ROLE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F15`, the next honest direct question is:

```text
does the current repo already export an explicit legacy-to-strict
role-equivalence verdict for the legacy fine-structure role?
```

## Result

The route is mixed but still negative at the theorem-relevant level:

```text
CURRENT_REPO_EXPORTS_STRICT_SIDE_ALPHA_EM_INV_MZ_CANDIDATE_BUT_NO_EXPLICIT_LEGACY_FINE_STRUCTURE_ROLE_EQUIVALENCE_VERDICT_AFTER_P82
```

## What was checked

`P82` checks two different things and keeps them separate:

1. whether the strict side exports a real candidate object carrying the modern
   fine-structure observable label,
2. whether the repo also exports an explicit semantic-transfer verdict saying
   that this object is the retained role-equivalent successor of the old
   legacy fine-structure formula.

## Why it only partially succeeds

On the current repo state:

1. `QW-2068`, `QW-2069`, `QW-2098`, and `QW-2094` together do export a
   real strict-side candidate object `alpha_em_inv_mz`,
2. but none of the current strict-side sources exports an explicit verdict
   that `alpha_em_inv_mz` is the retained role-equivalent successor of the
   legacy fine-structure semantics,
3. therefore the retained-side role-equivalence branch is narrowed but not yet
   discharged.

## Real reduction after `P82`

The retained frontier for the legacy fine-structure role is no longer:

```text
generic role-equivalence retention blocker
```

It is now:

```text
one explicit semantic-transfer blocker attached to alpha_em_inv_mz
```

namely:

```text
explicit_legacy_to_strict_semantic_transfer_verdict_identifying_alpha_em_inv_mz_as_the_retained_strict_side_successor_of_the_legacy_fine_structure_role
```

## What P82 does not claim

`P82` does not claim:

- that the retained branch is already discharged,
- that the replaced branch is already discharged,
- that `alpha_em_inv_mz` automatically inherits all legacy semantics,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. formalize theorem-level that the current repo exports a strict-side
   fine-structure candidate object but no explicit legacy-to-strict
   role-equivalence verdict,
2. keep the replaced branch separate,
3. do not silently treat `alpha_em_inv_mz` as if it already carried the full
   legacy `alpha_EM^-1` semantics.
