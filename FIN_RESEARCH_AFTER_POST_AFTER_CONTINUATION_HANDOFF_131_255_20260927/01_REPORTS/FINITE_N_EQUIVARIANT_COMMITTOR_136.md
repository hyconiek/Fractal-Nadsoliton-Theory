# FINITE-N-EQUIVARIANT-COMMITTOR-136
## Exact small-N hitting test for the ±k6 candidate under the leave-one-out Gibbs process

Date: 2026-09-26

Status:
- exact finite-state generator and exact committor solve for N=3..6;
- basin partition is D12-equivariant and generated from deterministic descent of the same V_g;
- this is a small-N process test, not a large-N capacity theorem.

Parameter:

    g = 5.145228719489142

(the saddle-level compromise point from report 131).

## 1. Basin construction

For each finite-N count state n:
1. map to the mean-field retained coordinate
       theta0=(g/N) X^T n;
2. descend the same full-X7 potential V_g;
3. classify the endpoint as
       +k6,
       -k6,
       or one of the 12 localized minima.

To prevent solver-induced symmetry bias:
- classify one representative per D12 orbit;
- propagate the label by the exact D12 action;
- if stabilizer ambiguity would assign different endpoints, leave the state
  unclassified/interior.

The resulting partition has exact + / - symmetry.

## 2. Exact committor problem

For the exact leave-one-out Gibbs generator L_N solve

    L h = 0

away from target sets with boundary conditions

    h=1 on -k6 basin,
    h=0 on all localized basins.

Then h is the exact probability of reaching -k6 before any localized basin.

Average h over equilibrium conditioned on the +k6 basin.

## 3. Results

### N=3

    state count: 364
    P(- before localized | + basin)
      ≈ 0.0106573

    P(localized before -)
      ≈ 0.989343.

### N=4

    state count: 1365
    P(- before localized | + basin)
      ≈ 0.0313902

    P(localized before -)
      ≈ 0.968610.

### N=5

    state count: 4368
    P(- before localized | + basin)
      ≈ 9.18427e-4

    P(localized before -)
      ≈ 0.999082.

### N=6

    state count: 12376
    P(- before localized | + basin)
      ≈ 3.49427e-4

    P(localized before -)
      ≈ 0.999651.

Thus at small N the exact process overwhelmingly exits toward the localized
sector before completing a parity switch.

## 4. This does not contradict the large-N saddle ordering

At g_opt the pure-k6 order parameter is

    q=tanh J
      ≈ 0.112554.

The +k6 mean even fraction is only

    (1+q)/2
      ≈ 0.556277.

So the count displacement from the uniform saddle is

    Delta Y
      = (q/2) N
      ≈ 0.056277 N.

Consequently:

    one count of separation:
      N≈17.8;

    five counts:
      N≈88.8;

    ten counts:
      N≈177.7.

For N<=6 the two mean-field wells are not even well resolved on the discrete
count lattice.

More importantly, the switching barrier is only

    B_switch≈1.3511e-5.

The large-deviation exponent reaches order one only near

    N≈7.4e4,

and reaches five near

    N≈3.7e5.

Therefore N=3..6 is far outside the asymptotic metastable regime.

## 5. Scientific consequence

Static saddle ordering is not enough.

The exact finite-N process shows that:
- discreteness;
- basin entropy;
- transition prefactors;
- and competing channels

dominate long before the large-deviation barrier exponent becomes effective.

This validates the review instruction to use capacities/committors rather than
assign transition rates from barriers alone.

## 6. Next capacity task

Direct state enumeration cannot reach the required N~10^5 scale.

The next step must use:
- D12 symmetry;
- potential-theory variational bounds;
- local Gaussian/Eyring-Kramers-type prefactors where justified;
- or large-deviation WKB/capacity asymptotics

to bridge the gap between exact small-N committors and the large-N saddle
picture.
