# F176 First Actual Strict Source Orientation Seed One-Sided Descent Support Packet

Status: `F176_EXECUTED_FIRST_ACTUAL_STRICT_SOURCE_ORIENTATION_SEED_ONE_SIDED_DESCENT_SUPPORT_PACKET_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `N286`, the repo already exports one stronger local support packet for
the first `T26` component:

```text
strict_source_orientation_seed_target_v1
```

The strongest honest next move is still not:

1. discharge of that component,
2. actual pair-indexed population anchor,
3. actual theta values,
4. actual `E_orient`,
5. admissible `S_sel_int`,
6. strict-core selector closure,
7. ToE closure.

It is narrower:

```text
freeze one actual one-sided local descent support packet
showing that the protected positive source branch admits
one positive-radius forward half-branch
on which the strict kernel remains positive
and strictly descends away from the declared source origin
```

This remains explicitly below component discharge.

## Fixed support reused

Reuse only already exported current-state support:

1. `T26/N284`
   - the first target component remains
     `strict_source_orientation_seed_target_v1`,
2. `F175/P266/N286`
   - one stronger local branch-polarity support packet is already exported,
3. `F163/N273`
   - one actual local derivative calculation witness is exported,
4. `F141/N249`
   - one actual source-side barrier-protected sign witness is exported,
5. `F140/N248`
   - one actual positive-radius local sign-stability witness is exported,
6. `N256/N257/N258`
   - the actual source-topology selector support family remains exported,
7. `H15/N13`
   - `K_obs` remains non-identified and may not be used as primary selector
     source,
8. `N234`
   - observer remains downstream only,
9. `K1/F2`
   - kernel split remains explicit and `K_legacy_ont` is not reactivated as
     live constructive kernel.

## Local analytic facts reused

Reuse the strict local kernel formula already written in `F163`:

```text
K_strict_gate(d) = cos(0.18575 * d + 0.16250) / (1 + d^1.8)
```

and its local derivative formula:

```text
K_strict_gate'(d)
=
[
  -0.18575 * sin(0.18575 * d + 0.16250) * (1 + d^1.8)
  - 1.8 * d^0.8 * cos(0.18575 * d + 0.16250)
] / (1 + d^1.8)^2
```

Together with:

```text
K_strict_gate'(0) ≈ -0.03004 < 0
```

Every term in the displayed derivative expression is continuous on
`[0, +infty)`. Therefore there exists one positive radius:

```text
epsilon_strict_src_one_sided_descent_derivative_radius_v1 > 0
```

such that:

```text
for all d in [0, epsilon_strict_src_one_sided_descent_derivative_radius_v1],
K_strict_gate'(d) < 0
```

Separately, from the already exported barrier-protected sign route, there
exists one positive source-side sign radius:

```text
epsilon_strict_src_positive_half_branch_radius_v1 > 0
```

such that:

```text
for all d in [0, epsilon_strict_src_positive_half_branch_radius_v1],
K_strict_gate(d) > 0
```

Freeze the combined one-sided radius:

```text
epsilon_strict_src_one_sided_descent_support_radius_v1
:=
min(
  epsilon_strict_src_one_sided_descent_derivative_radius_v1,
  epsilon_strict_src_positive_half_branch_radius_v1
)
> 0
```

## Actual one-sided descent support packet

Freeze one stronger actual support packet:

```text
Omicron_strict_source_orientation_seed_one_sided_descent_support_v1 :=
(
  strict_source_orientation_seed_target_v1,
  Xi_strict_source_orientation_seed_branch_polarity_support_v1,
  epsilon_strict_src_one_sided_descent_support_radius_v1,
  forward_source_half_branch_tag_v1,
  local_positive_half_branch_support,
  local_descending_half_branch_support,
  Upsilon_sel_basis_actual_witness_v1,
  Phi_qw2191_safe_actual_witness_v1,
  T14_src_selector_declared_scope_actual_witness_v1,
  Kobs_not_primary_selector_source,
  observer_downstream_only,
  kernel_split_safe_source_side_reading
)
```

with the scoped meaning:

```text
for all d with 0 < d <= epsilon_strict_src_one_sided_descent_support_radius_v1:
  K_strict_gate(d) > 0
  and
  K_strict_gate(d) < K_strict_gate(0)
```

The second clause follows from `K_strict_gate'(d) < 0` on the same interval.

## Honest meaning of the packet

This packet establishes only:

1. the first `T26` component is now supported not only by local branch
   polarity,
2. it is also supported by one distinguished forward half-branch on which the
   strict kernel remains positive and strictly decreases away from the source
   origin,
3. this is stronger than `N286`,
4. this is still below any actual source-orientation seed discharge.

## Why this is still below component discharge

`F176` remains below discharge because it still does not export:

1. one chart-independent seed object,
2. one pair-indexed source assignment,
3. one actual theta value,
4. one actual populated basis-pair instance,
5. one actual internal orientation datum,
6. one actual `E_orient`.

So the packet is one-sided local support, not discharge.

## Hard limits

`F176` does not discharge:

1. `strict_source_orientation_seed_target_v1`,
2. pair-indexed population anchor,
3. actual theta values,
4. actual populated basis-pair instance,
5. actual `E_orient`,
6. admissible `S_sel_int`,
7. actual strict-core selector closure,
8. actual global selector closure,
9. actual global `QW-2191` discharge,
10. actual ToE closure.

## Recommended next move

The correct next move is:

1. test whether the current repo really exports this one-sided local descent
   support packet,
2. keep the result explicitly below component discharge,
3. then decide whether one chart-independent source-seed projection support
   exists at all,
   or whether the exact remaining blocker should be frozen.
