# TH-005 — stopped stochastic-fold limit: derived candidate and exact blocker

## 1. Frozen double scaling

Let the supplied fold be `(theta_f,g_f,v)`. Fix a finite real lambda and set

`g_N=g_f+lambda N^(-2/3)`.

Define the soft observable

`Z_N(s)=N^(1/3) v.[theta_ret(N^(1/3)s)-theta_f]`,

stopped on exit from a fixed local fold neighbourhood. The microscopic process remains the same leave-one-out chain; only N and g are taken along an already-defined mathematical family.

## 2. Rescaled generator characteristics

A microscopic label jump i->j changes theta by

`Delta theta=(g_N/N)(X_j-X_i)`,

so `Delta Z=O(N^(-2/3))`. On s-time, the total event rate is O(N^(4/3)). Therefore:

- third and higher jump characteristics vanish (`third absolute characteristic = O(N^-2/3)` for bounded local states);
- the quadratic characteristic has an O(1) limit;
- the soft drift is multiplied by `N^(2/3)`.

Using TH-004 and the simple-fold Taylor expansion,

`N^(2/3) v.theta_dot -> -g_f[a lambda+(b/2)Z^2]`,

provided stable directions are negligible in the required averaged/perturbed-test sense. The quadratic characteristic tends to `2g_f` once the full occupation is locally driven to p_f.

The resulting candidate martingale problem is

`dZ=-g_f[a lambda+(b/2)Z^2] ds + sqrt(2g_f)dB_s`,

up to the declared stopping time.

## 3. What is already proved

The following ingredients are not conjectural:

1. simple fold and source intervals for g_f,a,b;
2. full path LDP of the microscopic chain (TH-002);
3. exact deterministic retained drift;
4. exact soft covariance identity `2g_f` at the fold;
5. leave-one-out O(1/N) drift correction, which becomes O(N^-1/3) under this scaling;
6. bounded jump sizes and vanishing third characteristic under the scaling;
7. a strictly stable six-dimensional transverse Hessian at the fold (source interval certificate plus numerical full-7D replay).

The Potts benchmark independently has an n^-2/3 spinodal window and n^(4/3) single-update critical scale. Converting its one-update clock by N updates per FIN model-time unit gives N^(1/3), matching the balance above. This is a universality diagnostic, not a transferred theorem.

## 4. Exact remaining lemma

A complete convergence theorem still requires a **FAST_MODE_AVERAGING_AND_CONTAINMENT lemma** for the actual jump chain:

For prepared local initial sequences and up to the stopped fold exit, prove that the six retained stable directions and the four occupation-fibre directions contribute `o_p(1)` to the rescaled soft martingale problem. Concretely, one needs a uniform time-integrated bound strong enough to make

`N^(2/3) integral [ |xi| ||y|| + ||y||^2 + hidden-dependent finite-N terms ] dt`

vanish after the s-time conversion, together with convergence of the soft quadratic characteristic to `2g_f`.

The mean-only path nonclosure from TH-003 prevents simply deleting the hidden occupation variables. Their deterministic fibre differences contract, and ordinary stable fluctuations are expected at O(N^-1/2), but that expectation is not substituted for the required stochastic estimate.

Katzenberger-style strong-drift methods may be useful, but the cited theorem is for SDEs and is not silently applied to this jump chain. A direct perturbed-test-function/martingale proof is the cleanest missing step.

## 5. Bounded finite-N diagnostic

For rounded fold occupations at N=120,240,480,960,1920,3840 the scaled soft variance approaches `2g_f`, the scaled absolute third jump characteristic decreases, and the drift residual behaves approximately as the predicted leave-one-out correction `(v.theta_f) N^-1/3`.

At lambda=0, `-drift_error*N^(1/3)` runs approximately

`2.663, 2.646, 2.637, 2.632, 2.630, 2.628`,

approaching the independently computed `v.theta_f≈2.627333894`. This is a diagnostic of the expansion, not a convergence proof.

## 6. Canonical dimensionless form

Let

`c=(4/b)^(1/3)`, `d=c^2/(2g_f)`, `Y=Z/c`, `tau=s/d`.

Then the candidate becomes

`dY=(Lambda-Y^2)d tau+dW`,

with

`Lambda=(-a c/2) lambda`.

Using the supplied a,b,g intervals,

- `c in [3.2274443594,3.2274456221]`,
- `d in [1.4814348295,1.4814359972]`,
- `(-a c/2) in [0.3430315477,0.3430317046]`.

Thus stopped first-exit and response laws, once the convergence lemma is paid, depend only on the dimensionless Lambda and normalized initial/boundary conditions. These are model-time combinations, not physical universal constants.
