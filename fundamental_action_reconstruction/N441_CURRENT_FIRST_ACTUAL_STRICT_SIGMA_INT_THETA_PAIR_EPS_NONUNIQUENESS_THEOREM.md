# N441 Current First Actual Strict Sigma-Int Theta-Pair Eps Nonuniqueness Theorem

Status: `N441_DISCHARGED_CURRENT_FIRST_ACTUAL_STRICT_SIGMA_INT_THETA_PAIR_EPS_NONUNIQUENESS_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

After `F314/N425` and `F325/N436`, the repo exports:

1. one strict-input positive-window theta-candidate record (for one chosen `eps`), and
2. one dedicated strict-side **candidate** selector-ingredient packaging of that record.

Both constructions rely on the sigma-int-driven `E_pair` generator contract:

```text
eps ∈ [0,1]   (T117)
```

This theorem packages one narrow anti-false-pass statement:

```text
on the current strict tuple and strict sigma-int lane,
the computed theta-pair depends on the admissible generator amplitude eps;
therefore eps is a real selector slot and no strict-core theta export/uniqueness
may be implied from any one chosen eps.
```

## Inputs reused (strict-admissible)

1. `T117`
   - sigma-int-driven `E_pair` generator contract with `eps ∈ [0,1]`.
2. `F307/N418`
   - strict sigma-int source upgrade (`sigma_int_strict_derived_v1 = -1`).
3. `F324/N435`
   - strict-provenance nad12 sign-mask value object (`b_sigma_int_E_pair_sign_mask_strict_provenance_v1`).
4. `F328/N440`
   - strict delta_d value object (`delta_d_sigma_int_positive_window_step_strict_provenance_v1 := delta_max`, premise-based).
5. `P407`
   - eps sensitivity audit (current strict tuple; fixed strict provenance inputs; varying `eps` inside `[0,1]`).

## Theorem-level conclusion (current repo state)

From `P407`, with the strict kernel tuple and strict-provenance inputs fixed,
the computed theta-pair varies over admissible `eps` choices inside `[0,1]`.

For example, `P407` records (current strict tuple; `sigma_int=-1`, strict-provenance sign mask fixed,
and fixed corridor-saturation `delta_d = delta_max`):

1. at `eps = 0.0`:
   - `theta_1^cand ≈ 0.3470114165377394`,
   - `theta_2^cand ≈ 0.3470114165377394`,
2. at `eps = 1.0`:
   - `theta_1^cand ≈ 0.38030796605778716`,
   - `theta_2^cand ≈ 0.320088837765114`.

The observed (unwrapped) theta variation exceeds numeric tolerance `1e-12`
on the sampled admissible grid.

Therefore the strongest honest meaning is:

1. `eps` is not eliminated by the current strict sigma-int → theta candidate pipeline;
   it is a real additional selector slot (a free generator input),
2. the strict sigma-int → theta construction remains candidate-level and cannot be
   promoted to strict-core theta export nor to actual object-support discharge
   without either:
   - an explicit eps selector premise, or
   - a genuinely new strict internal eps derivation/selection source.

## Consequence (scope hygiene)

This theorem blocks the following false-pass readings:

1. “the strict sigma-int lane uniquely fixes `eps`” (it does not),
2. “the `O(2)`-cut witness of `F325/N436` is canonical without extra slots”
   (it is conditional on the exported `delta_d` choice and on the exported `eps` premise; `QW-2191` remains open).

## Hard limits

`N441` does not claim:

1. actual strict-core `theta_1`, `theta_2` export,
2. admissible `S_sel_int` or strict-core selector closure,
3. `QW-2191` discharge,
4. discharge of post-`T148` object-support targets (`N302/N395/T130`),
5. ToE closure.
