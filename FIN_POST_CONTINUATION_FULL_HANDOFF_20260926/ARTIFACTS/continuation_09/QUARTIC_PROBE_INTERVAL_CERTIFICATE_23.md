# QUARTIC-PROBE-INTERVAL-CERTIFICATE-23 — outward-rounded global bound

Status: **RIGOROUS GLOBAL VALUE CERTIFICATE WITHIN GENERAL-QUARTIC-THEOREM-13**.

This closes the principal methodological defect in QUARTIC-PROBE-GLOBAL-19:
the previous branch-and-bound used ordinary double arithmetic.  The present
certificate derives the reduced four-amplitude polynomial symbolically and
uses outward-rounded coefficient enclosures plus an upper-bound-only
branch-and-bound.

## 1. Exact four-sector polynomial

After the proved phase-alignment reduction, write the nonnegative real sector
amplitudes as `c3,c4,c5,c6`.  Under `c3^2+c4^2+c5^2+c6^2=1`, the general quartic
theorem gives exactly ten monomials.  `quartic_probe_symbolic_coefficients.py`
derives them from exact Z12 trigonometry:

    c3^2 c4^2 : (g lambda3 + g lambda4 - 12)/72
    c3^2 c5^2 : (g lambda3 + g lambda5 - 12)/72
    c3 c4^2 c5: (g lambda3 + 2g lambda4 + g lambda5 - 24)/72
    c3 c4 c5 c6:
        sqrt(2)(g lambda3+g lambda4+g lambda5+g lambda6-24)/36
    c3 c5^3   : (g lambda3 + 3g lambda5 - 24)/144
    c4^2 c5^2 : (g lambda4 + g lambda5 - 12)/72
    c4^2 c6^2 : (g lambda4 + g lambda6 - 12)/36
    c4 c5^2 c6:
        sqrt(2)(g lambda4+2g lambda5+g lambda6-24)/48
    c5^4       : (g lambda5 - 6)/144
    c5^2 c6^2 : (g lambda5 + g lambda6 - 12)/36.

At `g_eq=3.7183448981203875` all ten coefficients are strictly positive.
Therefore aligned phases are globally maximizing and the problem reduces to
the simplex `y_k=c_k^2 >= 0`, `sum y_k=1`.

## 2. Rigorous coefficient enclosures

The strict eigenvalues were recomputed with 80-decimal `mpmath.iv` intervals
from the strict kernel, with the exact decimal kernel parameters and exponent
`1.8=9/5`.  Each positive polynomial coefficient was replaced in the C++
certificate by the next IEEE-754 binary64 number above its interval upper
endpoint.

A feasible aligned point was independently evaluated with interval arithmetic,
giving

    C >= 0.09324871588131618645946560103157...

under coefficient-sphere normalization.

## 3. Upper-bound branch-and-bound

`probe_bnb_interval_coeffsphere.cpp` uses `FE_UPWARD`, `-frounding-math`, and
`-ffp-contract=off`.  Branch endpoints are dyadic.  For every positive monomial
it uses the minimum of two independent valid upper bounds:

1. the product of coordinate upper endpoints;
2. weighted AM-GM under the remaining simplex mass.

All square roots, products and sums used in the upper path are evaluated under
upward rounding.  No numerical maximization of a monomial is used.

The exhaustive result is

    0.093248715881316186... <= Cmax
    Cmax <= 0.093298715881316108.

6,347,882 branch boxes were processed before all remaining boxes fell below the
certified upper cutoff.  Relative to the interval-certified pure unit-k5 value,

    5.2719896341 < Cmax/C5 < 5.2748164772.

Thus the mixed probe advantage over pure k5 is now a rigorous global-value
statement in the declared normalization.

## 4. Uniform-Fisher normalization

With `w_k=sqrt(lambda_k)c_k` and `sum w_k^2=1`, every coefficient is divided by
the corresponding product of strict eigenvalue powers.  These transformed
coefficients were again interval-enclosed before branch-and-bound.

The feasible interval point gives

    C >= 0.01839077538669661224215665071915...

and the exhaustive upward-rounded branch-and-bound gives

    Cmax <= 0.018400775386696603.

Using the rigorous pure-k5 interval gives the gain bound

    5.4936489337 < Cmax/C5 < 5.4966361101.

## 5. What is and is not certified

Certified:
- the global maximum **value** lies in the stated intervals;
- a mixed aligned probe is globally more than 5.27x stronger than pure k5 under
  coefficient norm, and more than 5.49x under fixed uniform-Fisher variance.

Not yet certified:
- uniqueness of the maximizing amplitude vector;
- a tiny interval box for its coordinates;
- any experimentally privileged normalization or apparatus.

No physical observable, clock, source for `g`, QW-2191, role transfer,
`L_total`, SM/GR or ToE conclusion follows.
