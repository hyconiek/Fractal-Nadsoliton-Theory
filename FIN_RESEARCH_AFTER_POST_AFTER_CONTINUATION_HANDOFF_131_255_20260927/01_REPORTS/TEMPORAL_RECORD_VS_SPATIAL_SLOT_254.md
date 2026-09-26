# TEMPORAL-RECORD-VS-SPATIAL-SLOT-254
## Minimum-capacity factorization produces record coordinates, not yet physical spatial sites

Date: 2026-09-26

Report 253 derives

    E ≅ Z3^m

for a minimum-capacity environment supporting m perfect resets.

This is substantial progress, but the semantics of those factors must be kept
precise.

## 1. What the factor coordinates mean operationally

The natural coordinates are

    (Y_1,...,Y_m),

where Y_t is the future trit delivered on reset number t.

Thus the factorization is canonically indexed by RESET ORDER.

It is a bank of future records.

## 2. Why this does not yet give physical space

All m coordinates coexist mathematically in the initial environment state, but
their distinguished operational role is temporal:

    slot 1 -> first future reset,
    slot 2 -> second future reset,
    ...

Nothing in the entropy theorem defines:
- adjacency between record coordinates;
- metric distance;
- simultaneous local interaction;
- spatial dimension.

Any permutation of the m record coordinates preserves:
- total entropy;
- factor count;
- reset information content.

So the product decomposition is real, while spatial ordering remains
underdetermined.

## 3. What the cyclic scheduler adds

The single-cycle permutation of reports 250-251 orders the records as

    Y_1 -> Y_2 -> ... -> Y_m -> recurrence.

This is an internally defined incidence/routing relation.

But it is selected because it schedules FUTURE resets efficiently.

Its primary meaning is therefore still memory/routing order.

Calling the same cycle physical space requires an additional role-transfer
argument.

## 4. Updated typed boundary

We can now distinguish three levels:

### Information factorization — DERIVED at minimum capacity

    E ≅ Z3^m.

### Record order — CONDITIONALLY sourced by maximal freshness

    one cycle on the m records.

### Physical spatial adjacency — STILL OPEN

    why record-neighbor relation
      =
    simultaneous physical-neighbor relation.

This is a much narrower and cleaner remaining gap than the earlier generic
"where do slots come from?" question.
