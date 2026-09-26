# RIGOROUS-MEMORY-TAIL-BOUND-234
## Reversibility gives sector-wise exponential upper bounds on the unresolved memory tail

Date: 2026-09-26

Status:
proof-grade bound conditional on the numerically resolved N=6 coupled hidden
spectral gaps.

For one resolved Fourier sector,

    K_k(t)
      =
      <c_k,exp(-H_k t)c_k>,

with

    H_k=-QSQ >=0.

Only hidden eigenmodes with nonzero overlap with c_k contribute.

Define the COUPLED hidden gap

    gamma_k
      =
      inf{
        gamma:
        spectral weight of c_k at gamma is nonzero
      }.

Then the positive spectral representation immediately gives

    boxed:
    0 <= K_k(t)
      <=
    K_k(0) exp(-gamma_k t).

Integrating:

    boxed:
    integral_T^infinity K_k(t) dt
      <=
    [K_k(0)/gamma_k]
    exp(-gamma_k T).

This is a rigorous memory-tail inequality once gamma_k is bounded.

## 1. N=6 coupled gaps

The numerical symmetry-resolved gaps are:


    k=1:
      gamma_k≈0.809540687259

    k=2:
      gamma_k≈0.824868743147

    k=3:
      gamma_k≈0.681107123762

    k=4:
      gamma_k≈0.677188233336

    k=5:
      gamma_k≈0.683980312936

    k=6:
      gamma_k≈0.682245144156


The global hidden spectrum contains slower modes near 0.48-0.68, but many are
exactly or numerically symmetry-dark for a given memory source.

Therefore the relevant decay scale is the coupled gap, not the smallest hidden
eigenvalue of the whole Q-space.

## 2. Tail fraction relative to M0

Using

    M0=integral_0^infinity K(t)dt,

the rigorous upper bound on the remaining tail fraction at T=8 is:


    k=1:
      <= 0.834 %

    k=2:
      <= 0.717 %

    k=3:
      <= 2.936 %

    k=4:
      <= 3.003 %

    k=5:
      <= 2.786 %

    k=6:
      <= 2.849 %


Thus after eight microscopic clock units the unresolved tail contains at most
about 3.0% of M0 in every sector under the computed N=6 gap estimates.

At T=4 the bound is much looser, roughly 19-45%.

The bound is conservative; direct semigroup tests show substantially better
effective behavior.

## 3. Convolution error meaning

If the resolved observable satisfies

    |u(t)| <= U

over the relevant interval, then neglecting memory older than T changes the
convolution term by at most

    U [K(0)/gamma_k] exp(-gamma_k T).

So the memory horizon can now be selected from an explicit tolerance rather
than by visual inspection.

## 4. Relation to the one-pole model

The moment-matched one-pole decay rates are around 2.94-3.28, much faster than
the rigorous coupled gaps ~0.68-0.82.

There is no contradiction.

The one-pole rate is a weighted average emphasizing the memory spectral weight;
the rigorous gap is controlled by the slowest nonzero tail, even if that tail
has very small amplitude.

Thus:
- one-pole model predicts practical dynamics;
- coupled-gap bound certifies worst-case tail decay.

Both are useful and answer different questions.
