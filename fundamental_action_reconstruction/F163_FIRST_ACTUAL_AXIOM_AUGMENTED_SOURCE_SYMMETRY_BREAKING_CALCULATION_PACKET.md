# F163 First Actual Axiom-Augmented Source Symmetry-Breaking Calculation Packet

Status: `F163_EXECUTED_FIRST_ACTUAL_AXIOM_AUGMENTED_SOURCE_SYMMETRY_BREAKING_CALCULATION_PACKET_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `N272`, the strongest honest next move is still not:

1. actual non-strict declared-scope ToE closure,
2. actual strict-core selector closure,
3. actual global `QW-2191` discharge.

It is narrower:

```text
freeze one actual source-side local derivative calculation
at the strict topological source limit
```

and keep it explicitly below:

1. selector-datum discharge,
2. basis-independent promotion,
3. quotient-safe `QW-2191` resolution,
4. actual non-strict declared-scope ToE closure.

`F163` executes exactly that move.

## The mathematical calculation

Reuse the strict source-side working kernel already admitted on the `T14`
route:

```text
K_strict_gate(d) = cos(0.18575 * d + 0.16250) / (1 + 1.0 * d^1.8)
```

At the declared source limit `d -> 0`:

### 1. Value at the origin

```text
K_strict_gate(0) = cos(0.16250) ≈ 0.9868 > 0
```

### 2. First derivative at the origin

Differentiate with respect to `d`:

```text
d/dd K_strict_gate(d)
=
[
  -0.18575 * sin(0.18575 * d + 0.16250) * (1 + d^1.8)
  - 1.8 * d^0.8 * cos(0.18575 * d + 0.16250)
] / (1 + d^1.8)^2
```

Taking the strict limit as `d -> 0`:

1. the damping derivative term vanishes,
2. the denominator tends to `1`,
3. the limit resolves to:

```text
K_strict_gate'(0) = -0.18575 * sin(0.16250) ≈ -0.03004 != 0
```

## Honest meaning of the calculation

This calculation establishes only:

1. one explicit local source-side derivative value is now written at the
   declared origin,
2. the declared source limit is not a flat critical point for this charted
   kernel expression,
3. one concrete source-side asymmetry calculation is therefore available for
   later guarded use on the non-strict lane.

This does **not** yet establish:

1. that the derivative itself is a selector datum,
2. that the derivative yields a basis-independent source selector,
3. that the derivative resolves the `QW-2191` quotient frontier,
4. that the derivative closes the non-strict declared-scope ToE lane.

In particular, the current calculation remains below the requirements still
listed in `T14`:

1. basis-independent promotion,
2. quotient-safe selector resolution,
3. theorem-level closure integration.

## Relation to existing source-side packets

`F163` must be read only as a local calculation packet extending the
source-side analysis.

It does not replace:

1. `tau_src_candidate_v1`,
2. `T_flow^(0)`,
3. `Pi_sel`,
4. `Phi_qw2191_safe_actual_witness_v1`.

It only adds one more explicit source-side chart calculation that may later be
used as a candidate supporting datum in `axiom_augmented_only` scope.

## Result

`F163` exports one actual local source-derivative calculation witness:

```text
chi_src_local_derivative_calculation_witness_v1 :
  K_strict_gate'(0) = -0.18575 * sin(0.16250) ≈ -0.03004 != 0
```

with the declared properties:

1. source-side only,
2. local-origin calculation only,
3. compatible with later axiom-augmented interpretation,
4. below selector discharge,
5. below `QW-2191` discharge,
6. below ToE closure.

## Hard limits

`F163` does not discharge:

1. actual non-strict declared-scope ToE closure,
2. actual strict-core selector closure,
3. actual global selector closure,
4. actual global `QW-2191` discharge,
5. actual legacy-to-strict bridge derivation,
6. actual strengthened nonbridge theorem.

## Recommended next move

The correct next move is:

1. test whether the current repo really exports this local derivative
   calculation packet,
2. keep the result below selector/QW-2191/ToE discharge,
3. do not relabel the derivative as an already discharged selector source.
