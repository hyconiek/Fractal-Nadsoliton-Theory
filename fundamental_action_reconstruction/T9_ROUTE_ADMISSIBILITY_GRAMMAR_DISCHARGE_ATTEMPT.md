# T9 Route Admissibility Grammar Discharge Attempt

Status: `T9_EXECUTED_T8_DISCHARGE_ATTEMPT_BLOCKED_BY_NO_FORMAL_ROUTE_ROLE_TYPING_RULE`
As of: `2026-03-06`

## Goal

After `T8`, the natural next move is an actual discharge attempt for:

```text
T8_Route_Admissibility_Grammar_Theorem
```

`T9` does not write another theorem spec.
It checks whether the present selector-track audit vocabulary already defines
route roles strongly enough to induce a finite admissibility grammar for actual
strict-core exports of local phase values `theta_1`, `theta_2`.

`T9` does **not** claim that `T8` is discharged.

## Target being tested

```text
For the current strict-core selector track,
every admissible actual-theta export route is typed by exactly one route role:
  raw_overlap / formula / representative / downstream_schema /
  source_skeleton / strict_to_axiom_bridge.
Hence every admissible route instantiates one of the six audited
constructor families.
```

## Strict-admissible evidence used

1. `C32`
   - raw overlap route explicitly isolated
2. `C33`
   - formula-class route explicitly isolated
3. `C34`
   - representative-class route explicitly isolated
4. `C49`
   - downstream populated-instance route explicitly isolated
5. `C50`
   - strict-core source-skeleton route explicitly isolated as absent
6. `C51`
   - strict-to-axiom bridge route explicitly isolated as fallback-only
7. `T8`
   - admissibility grammar already specified
8. `A10`
   - anti-overclaim boundary

## Discharge attempt

### Step 1. Check explicit role vocabulary

The current selector track already exposes six explicit route roles:

1. raw overlap,
2. formula,
3. representative,
4. downstream schema,
5. source skeleton,
6. strict-to-axiom bridge.

This part of `T8` is in good shape.

### Step 2. Check route-role stability across audits

Each currently audited theta-export route has already been classified by one of
those roles, and the role names stay stable across the current selector track.

This supports the weaker statement:

```text
all currently named theta-export routes carry one of six explicit route-role labels.
```

### Step 3. Test the missing typing step

The remaining question is narrower:

```text
Does the current strict-core selector track export a formal rule saying that
route admissibility is completely determined by route role, so that every
current admissible theta-export route must instantiate exactly one of the six
named roles?
```

The current repo state gives:
- a finite explicit role vocabulary,
- route-wise audits,
- stable role labels,
- anti-overclaim constraints,

but it does **not** yet give:
- a formal route-role typing rule,
- a uniqueness-of-role declaration for admissible routes,
- or a theorem-level statement that admissibility is exhausted by those six
  route-role labels.

## Strongest honest conclusion after T9

After `T9`, the strongest honest conclusion is:

- `T8` is **not discharged**,
- not because the route-role vocabulary is unclear,
- but because the repo still lacks a formal route-role typing rule proving that
  admissibility is exhausted by the six named route roles.

## Residual blocker found by the discharge attempt

```text
T9_B1 :=
no formal route-role typing rule or admissibility-by-role declaration
showing that every current strict-core theta-export route
must instantiate exactly one of the six named route roles
  {raw_overlap, formula, representative, downstream_schema,
   source_skeleton, strict_to_axiom_bridge}.
```

## Reduction of the theorem-lane frontier

Before `T9`:

- `T8_B1 := the route admissibility grammar is specified but not discharged for the current selector track`
- `T2_B1 := the bridge theorem is specified but not discharged; strict-core target slot and sign-only export-map object exist, but target-slot population (theta_1,theta_2) remains absent`
- `C32_B2 := raw overlap scalar route remains degenerate`

After `T9`:

- `T9_B1 := no formal route-role typing rule or admissibility-by-role declaration showing that every current strict-core theta-export route must instantiate exactly one of the six named route roles`
- `T2_B1 := the bridge theorem is specified but not discharged; strict-core target slot and sign-only export-map object exist, but target-slot population (theta_1,theta_2) remains absent`
- `C32_B2 := raw overlap scalar route remains degenerate`

## What T9 does not claim

`T9` does not claim:
- theorem-level PASS,
- full-closure PASS,
- that `T8` is discharged,
- that route-role labels already form a proved admissibility grammar,
- that no future theory extension can add a new admissible route role,
- that `QW-2191` is discharged.

## Product of the step

- first real discharge attempt for `T8`,
- explicit reduction of the failure to one new meta-level blocker,
- maintained no-false-pass discipline.

## Next step

Natural next move:
- write a theorem spec for the missing route-role typing rule /
  admissibility-by-role declaration,
- or explicitly construct a finite role-typing rule for the current selector
  track.
