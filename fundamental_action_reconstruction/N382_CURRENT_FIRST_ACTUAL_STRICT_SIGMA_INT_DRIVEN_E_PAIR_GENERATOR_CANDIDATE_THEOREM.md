# N382 Current First Actual Strict Sigma-Int Driven `E_pair` Generator Candidate Theorem

Status: `N382_DISCHARGED_CURRENT_FIRST_ACTUAL_STRICT_SIGMA_INT_DRIVEN_E_PAIR_GENERATOR_CANDIDATE_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

Advance Path 1 beyond the `E_pair` template-only state (`N381`) by exporting
one explicit **internal-datum-driven** (still candidate-level) generator for
`E_pair` that is:

1. noncyclic (`N18`-safe: no `theta` input and no populated-instance input),
2. observer-free (`AX9` discipline: no `K_obs` primary selection),
3. pair-indexed on `[pair1,pair2]`.

## Theorem-level conclusion

From `T117/F270`, the current repo exports one actual packaged generator
candidate:

```text
G_sigma_int_to_E_pair_generator_candidate_v1
```

with the following exact meaning:

1. it takes as input the strict-core internal datum candidate
   `sigma_int_input ∈ {+1,-1}`,
2. it outputs a finite normalized carrier field `E_pair` on `[pair1,pair2]`,
3. it is explicitly noncyclic and observer-free by contract,
4. it remains candidate-level (no strict derivation, no uniqueness).

No identification theorem between `sigma_int_candidate` (`B4`) and
`sigma_int_strict_derived_v1` (`F307/N418`) is used or implied here; both are
admissible instantiations of the abstract input `sigma_int_input ∈ {+1,-1}`.

## Operational check (scope-limited)

For one admissible parameter choice `eps = 1/2`, and for `sigma_int_input`
equal to `+1` or `-1`, the resulting `E_pair` can be plugged into the
reduction form `T115` without hitting the `(X_i^cand,Y_i^cand)=(0,0)`
degeneracy frontier on either pair slot.

One operational evaluation gives:

```text
sigma_int_input = +1, eps = 1/2:
  pair1: theta^cand ≈ 0.3554808343611710
  pair2: theta^cand ≈ 0.4625142242896770

sigma_int_input = -1, eps = 1/2:
  pair1: theta^cand ≈ 0.4625142242896770
  pair2: theta^cand ≈ 0.3554808343611710
```

This is an operational check only. It is not:

- a strict proof of global nondegeneracy,
- a strict proof of selector uniqueness,
- a discharge of `QW-2191`.

## What N382 does not prove

`N382` does not prove:

1. strict derivation or uniqueness of the generator,
2. discharge of the strict-core sigma bridge theorem `T2`,
3. resolution of the residual object-support incompatibility boundary `N302`,
4. actual strict-core `theta_1`, `theta_2`,
5. actual pair population,
6. strict-core selector closure,
7. `QW-2191` discharge,
8. ToE closure.

## Consequence (next honest step)

After `N382`, the next honest move is no longer “define an `E_pair` generator”.

It is to connect the generated carrier field `E_pair` to:

1. an explicit carrier/projection interface attacking the `N302` object-support
   layer (without cyclic reuse),
2. and then either:
   - strengthen that interface toward an actual object-support witness, or
   - prove that such an admissibility upgrade fails under strict constraints.

The current repo performs this next move only at **candidate projection**
level through `N383` and remains below actual object support.
