# F8 Legacy Weinberg Role-Equivalence Refinement Packet

Status: `F8_EXECUTED_LEGACY_WEINBERG_ROLE_EQUIVALENCE_REFINEMENT_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `P66/N69`, the retained-side frontier for the legacy Weinberg-angle role
is no longer literal retention. The only remaining retained-side question is:

```text
explicit_strict_side_role_equivalence_verdict_for_the_legacy_weinberg_angle_role
```

`F8` asks the next honest refinement question:

```text
does that remaining role-equivalence frontier reduce to
1. the presence of some strict-side Weinberg candidate object,
or instead to
2. one still-missing explicit semantic-transfer verdict from legacy to strict?
```

## Result

`F8` establishes the following refinement on the current repo state:

1. the strict side already exports a real candidate object,
   `sin2_theta_w_mz`,
2. that candidate is not yet the same thing as an explicit retained-role
   verdict for the legacy Weinberg-angle semantics,
3. therefore the remaining retained-side blocker is now narrowed to one
   explicit semantic-transfer / role-equivalence verdict linking the legacy
   role to the strict-side candidate object.

## Why this follows

The split is forced by current repo evidence:

1. `QW-2068` registers `sin2_theta_w_mz` as a target parameter,
2. `QW-2069` exports `sin2_theta_w_mz` as a strict-derived package entry,
3. `QW-2098` derives `sin2_theta_w_mz` on the strict non-anchor EW chain,
4. Release `4.9` explicitly promotes `sin2_theta_w_mz` on the strict side,
5. but the current repo still exports no explicit verdict saying that this
   strict-side object is the retained role-equivalent successor of the old
   legacy Weinberg formula.

## Real reduction after `F8`

So the retained-side frontier is no longer:

```text
generic role-equivalence retention blocker
```

It is now:

```text
one explicit semantic-transfer blocker
```

namely:

```text
explicit_legacy_to_strict_semantic_transfer_verdict_identifying_sin2_theta_w_mz_as_the_retained_strict_side_successor_of_the_legacy_weinberg_angle_role
```

## What F8 does not claim

`F8` does not claim:

- that the retained branch is already discharged,
- that the replaced branch is already discharged,
- that `sin2_theta_w_mz` automatically inherits all legacy semantics,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. probe directly whether the current strict-side materials export that one
   missing semantic-transfer verdict,
2. while keeping the replaced branch separate and without silently transferring
   legacy semantics onto `K_strict_gate`.
