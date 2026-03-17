# N715 Current Strict `T176` Low-Complexity Odd-Polynomial Two-Readout Provider-Class Boundary Theorem

Status: `N715_CURRENT_STRICT_T176_LOW_COMPLEXITY_ODD_POLYNOMIAL_TWO_READOUT_PROVIDER_CLASS_BOUNDARY_THEOREM_NO_FALSE_PASS`
As of: `2026-03-18`

## Goal

Package the next honest nonlinear provider attack after `P718/N714`, while
still forbidding hidden selector tuning.

## Theorem-level conclusion

On the current repo state:

1. `P718/N714` already rules out the whole single mixed linear span
   `span{w_break, w_ref_unnormalized}`,
2. `P719` then scans the nearest untuned nonlinear extension of that same
   current two-readout carrier:
   odd polynomials in
   `(<w_break,u_root>, <w_ref,u_root>)`
   of total degree `<= 3`, with coefficient alphabet `{-1,0,1}`,
3. across all `728` nonzero candidates in that class, the repo exports:
   - `0` exact all-root directed-section candidates,
   - `576` projective-only candidates,
   - `152` candidates that fail all-root support,
   - and only the already-known negated-root patterns
     `["pair4"]` or `["pair2","pair3"]`.

Therefore the strongest honest current statement is:

```text
the minimal untuned low-complexity nonlinear extension of the current
two-readout carrier still does not discharge T176.
```

## Consequence

The provider frontier narrows again:

1. one must now go beyond both:
   - the single mixed linear span class (`P718/N714`), and
   - the low-complexity untuned odd-polynomial two-readout class (`P719/N715`),
2. so the next strict provider attempt must either:
   - raise complexity in a still-explicit, non-hidden way, or
   - leave the present two-readout carrier and introduce a genuinely new
     observable/provider class.
3. from the current Release-7 physics perspective, the preferred continuation is
   the second route: move toward a new **physically interpretable**
   observable/provider class rather than escalating coefficient tuning inside
   the current abstract two-readout algebra.

## What this theorem does not prove

`N715` does **not** prove:

1. `T176` discharge,
2. a strict-core directed/sign-sensitive physical orientation datum,
3. impossibility of all tuned higher-degree or non-polynomial candidates,
4. impossibility of all future providers using new observables,
5. kernel-alone/global `QW-2191` discharge,
6. ToE closure.
