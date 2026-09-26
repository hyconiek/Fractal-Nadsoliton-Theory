# PURE-K6-SUBCRITICAL-CASCADE-68
## The R7P-095 k3 instability is a subcritical pitchfork feeding the MP7-039 two-harmonic saddle branch

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Inputs:
- R7P-095 first pure-k6 angular instability;
- MP7-039 two-harmonic D3 transverse crossing;
- strict X7 normalization.

Status:
- exact symmetry/branch identities;
- numerical Lyapunov-Schmidt coefficients at nominal strict spectrum;
- numerical continuation connects the two previously separate repo events;
- coefficient signs are not yet interval-certified.

## 1. Pure k6 stationary branch

Write the pure alternating field as

    h_j = K (-1)^j.

With

    alpha6^2=lambda6/12,

the dual coordinate is

    s6=K/alpha6.

The stationary equation is exactly

    K = g (lambda6/12) tanh K,

hence

    boxed:
    g(K)=12 K/[lambda6 tanh K].

The nonzero branch emerges continuously from the uniform state at

    g_u = 12/lambda6
        ~ 5.12342755140.

## 2. The R7P-095 angular crossing is also a full-Hessian crossing

The k3-cosine full dual Hessian eigenvalue along the pure k6 branch is

    lambda_c(K)
      =
      1/g(K)
      -lambda3[1+tanh K]/12.

Therefore lambda_c=0 is exactly

    lambda3(1+tanh K)
      =
      lambda6 tanh K/K,

which is the R7P-095 crossing equation after identifying

    K=sqrt(lambda6/12) r.

So the certified angular event is not merely an unrelated spherical artifact:
on the stationary pure-k6 branch it is the exact zero of a full Phi Hessian
direction.

Using the nominal strict spectrum and the R7P-095 radius midpoint gives

    K_p ~ 0.182996029533,
    g_p ~ 5.18049061964.

At this point all six other full Hessian eigenvalues are positive:

    0.001481387511 (x2),
    0.004243035819,
    0.009734506073 (x2),
    0.059162676706.

Thus the pure-k6 state is a local minimum immediately below the crossing and
becomes index one immediately above it.

## 3. Cubic is symmetry-forbidden

The pure k6 state is invariant under translation by two labels.

On the k3-cosine critical amplitude x this symmetry acts as

    x -> -x.

Therefore the reduced potential is exactly even in x and every odd coefficient,
including the cubic, vanishes.

Direct cumulant evaluation gives

    D3 Phi[e,e,e] ~ 6e-17,

consistent with exact zero.

## 4. Reduced quartic is negative

Let e be the normalized k3-cosine critical direction and eliminate the six
noncritical variables by the stationary equations.

The direct fourth derivative is positive:

    D4_direct ~ 0.0486816028457.

However the cubic coupling of x^2 to the stable k6 radial mode is strong.

The Lyapunov-Schmidt formula gives

    D4_red
      =
      D4_direct
      -3 <T(e,e),H_s^{-1}T(e,e)>

      ~ -3.40061736739.

Therefore

    boxed:
    beta = D4_red/24
         ~ -0.141692390308 < 0.

So the symmetry-protected pitchfork is **subcritical**, not supercritical.

This is an important warning:
forbidding a cubic is necessary but not sufficient for stable daughter minima.

## 5. Soft eigenvalue slope

Along the pure-k6 stationary branch,

    lambda_c(g)
      =
      ell (g-g_p)+...

with

    boxed:
    ell ~ -0.291327362922.

Hence:
- for g<g_p, lambda_c>0 and the pure-k6 state is stable in this direction;
- for g>g_p, lambda_c<0 and it is unstable.

## 6. Subcritical daughter branch

With

    V_red
      =
      lambda_c x^2/2
      + beta x^4
      + ...

and beta<0, the small nonzero stationary branch exists on the side
lambda_c>0, i.e. for

    g<g_p.

Its leading amplitude is

    x^2
      ~ -lambda_c/(4 beta),

so

    x
      ~ A_theta sqrt(g_p-g),

with

    A_theta ~ 0.716947541263.

In the unscaled k3 field

    J=sqrt(lambda3/6) x,

this becomes

    boxed:
    J
      ~ 0.409916688529 sqrt(g_p-g).

Because the daughter curvature in the x direction is

    -2 lambda_c <0,

this newborn branch is index one.

It is a saddle branch, not a minimum.

## 7. Direct connection to MP7-039

MP7-039's two-harmonic branch point is

    J = 0.038434753979,
    K = 0.185272061402,
    g = 5.171841831943.

The gain distance from the pure-k6 pitchfork is

    g_p-g
      ~ 0.008648787695.

The leading subcritical formula predicts

    J_pred
      ~ 0.0381218,

already within about one percent of the actual MP7-039 J despite the event not
being asymptotically infinitesimal.

Thus the MP7-039 family is the continuation of the subcritical k3 daughter
branch born at the pure-k6 instability.

This unifies R7P-095 and MP7-039 into one bifurcation chain.

## 8. Index cascade along the k3+k6 branch

Numerical continuation of the exact two-harmonic stationary equations gives:

### Near g_p from below
The small-J branch has full Hessian index 1.

### At MP7-039, g ~ 5.17184183194
A D3-symmetric double transverse pair crosses.

The branch changes:
    index 3 for g below the crossing,
    index 1 for g above the crossing.

Report 66 shows that the emitted D3 daughter branches have index 2.

### Radial two-harmonic fold
The same base branch has a later stationary fold at approximately

    J_f ~ 0.718505246006,
    K_f ~ 0.490908854614,
    g_f2 ~ 4.621196599489.

Across that fold its index changes between 3 and 2.

No segment found in this continuation becomes a stable minimum.

## 9. Two distinct k3+k6 roots at g=5

At g=5 the two-harmonic family has two solutions.

Small-amplitude branch:
    J ~ 0.20709230,
    K ~ 0.23524693,
    full index 3.

Large-amplitude branch:
    J ~ 1.24494800,
    K ~ 0.77955551,
    full index 2.

The large branch matches the symmetry-enhanced k3+k6 stationary orbit already
visible in the g=5 numerical atlas.

## 10. Emergence implication

The high-symmetry sequence is therefore:

    pure k6 local minimum
      |
      | subcritical k3 pitchfork
      v
    k3+k6 index-1 saddle
      |
      | D3 double crossing
      v
    k3+k6 index-3 saddle
      + index-2 D3 daughter saddles
      |
      | two-harmonic fold
      v
    index-2 saddle branch.

This whole high-symmetry structure organizes the saddle network but does not
produce the globally stable localized phase.

That stable localized phase was born much earlier, at the separate
four-amplitude saddle-node near g=3.51564.

So FIN's first stable localization is a finite-amplitude fold phenomenon, not
a continuous instability of a pure Fourier mode.

## 11. Next atom

### METASTABILITY-HIERARCHY-69

Combine:
- localized fold g~3.51564;
- first global coexistence g~3.71834;
- uniform k6 spinodal g=12/lambda6;
- pure-k6 subcritical loss g~5.18049;
- subsequent D3 saddle restructuring.

Goal:
produce one rigorously scoped phase/metastability hierarchy and identify which
events change global equilibrium, local stability, or only saddle topology.
