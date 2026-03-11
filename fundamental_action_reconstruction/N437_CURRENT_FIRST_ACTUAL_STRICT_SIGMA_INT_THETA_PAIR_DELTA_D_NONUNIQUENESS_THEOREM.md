# N437 Current First Actual Strict Sigma-Int Theta-Pair Delta_d Nonuniqueness Theorem

Status: `N437_DISCHARGED_CURRENT_FIRST_ACTUAL_STRICT_SIGMA_INT_THETA_PAIR_DELTA_D_NONUNIQUENESS_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After `F314/N425` and `F325/N436`, the repo exports:

1. one strict-input positive-window theta-candidate record (for one chosen `delta_d`), and
2. one dedicated strict-side **candidate** selector-ingredient packaging of that record.

Both constructions rely on the typed positive-window corridor:

```text
delta_d ∈ (0, d_local/11]    (T119)
```

This theorem packages one narrow anti-false-pass statement:

```text
on the current strict tuple and strict sigma-int lane,
the computed theta-pair depends on the admissible corridor step delta_d;
therefore delta_d remains a real selector slot and no strict-core theta export/uniqueness
may be implied from any one chosen delta_d.
```

## Inputs reused (strict-admissible)

1. `T119`
   - positive-window corridor with free step `delta_d ∈ (0, d_local/11]`.
2. `F307/N418`
   - strict sigma-int source upgrade (`sigma_int_strict_derived_v1 = -1`).
3. `F317/N428`
   - strict eps value object (`eps_sigma_int_E_pair_amplitude_strict_provenance_v1 := 1/2`, premise-based).
4. `F324/N435`
   - strict-provenance nad12 sign-mask value object (`b_sigma_int_E_pair_sign_mask_strict_provenance_v1`).
5. `P403`
   - delta_d sensitivity audit (current strict tuple; fixed strict provenance inputs; varying `delta_d` inside the corridor).

## Theorem-level conclusion (current repo state)

From `P403`, with the strict kernel tuple and strict-provenance inputs fixed,
the computed theta-pair varies over admissible `delta_d` choices inside the
positive-window corridor.

For example, `P403` records (current strict tuple; `sigma_int=-1`, `eps=1/2`,
and the strict-provenance sign mask fixed):

1. at `delta_d = delta_max/12 ≈ 0.028718469896710636`:
   - `theta_1^cand ≈ 0.19246023489971711`,
   - `theta_2^cand ≈ 0.18984234925310361`,
2. at `delta_d = delta_max ≈ 0.34462163876052765`:
   - `theta_1^cand ≈ 0.3627333053541785`,
   - `theta_2^cand ≈ 0.33287066305007096`.

The observed (unwrapped) theta variation exceeds numeric tolerance `1e-12`
on the sampled admissible grid.

Therefore the strongest honest meaning is:

1. `delta_d` is not eliminated by the positive-window corridor; it is a real
   additional selector slot (a free corridor input),
2. the strict sigma-int → theta construction remains candidate-level and cannot
   be promoted to strict-core theta export nor to actual object-support discharge
   without an explicit `delta_d` selector premise or a genuinely new strict
   internal `delta_d` selector source.

## Consequence (scope hygiene)

This theorem blocks the following false-pass readings:

1. “the positive-window corridor uniquely fixes `delta_d`” (it does not),
2. “the `O(2)`-cut witness of `F325/N436` is kernel-alone” (it is conditional on the
   recorded `delta_d` choice in the exported artifact; `QW-2191` remains open).

## Hard limits

`N437` does not claim:

1. actual strict-core `theta_1`, `theta_2` export,
2. admissible `S_sel_int` or strict-core selector closure,
3. `QW-2191` discharge,
4. discharge of post-`T148` object-support targets (`N302/N395/T130`),
5. ToE closure.

