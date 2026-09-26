# D3-NONLINEAR-DAUGHTER-BRANCHES-66
## The MP7-039 transverse crossing is a genuine cubic D3 bifurcation; its daughter branches are not new minima

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Input:
- `MP7-039_boundary_counterexample_mechanism.md`;
- `MP7-039_transverse_crossing.json`;
- the strict X7 feature normalization already used throughout MP7.

Status:
- analytic normal-form reduction plus numerical coefficient evaluation at the
  certified crossing midpoint;
- direct 7D stationary replay on both sides of the crossing;
- NOT yet an interval certificate for the cubic/quartic coefficient signs.

## 1. Certified linear crossing

MP7-039 gives the two-harmonic stationary branch

    h_j = J cos(pi j/2) + K (-1)^j

and a simple transverse crossing at

    J* = 0.038434753978898...
    K* = 0.185272061401889...
    g* = 5.171841831942818...

with two-dimensional critical space.

The base isotropy is D3.

The critical representation is the standard real 2D representation of D3.

## 2. Critical basis

Choose orthonormal critical coordinates `(x,y)` so that:

- x lies in the `(k4 cosine, k5 cosine)` block;
- y lies in the symmetry-related `(k4 sine, k5 sine)` block.

At the numerical crossing midpoint one convenient orientation is

    e_x:
      k4c = 0.390006185096
      k5c = 0.920812236879

    e_y:
      k4s = 0.390006185096
      k5s = -0.920812236879.

The remaining five Hessian directions are noncritical.

Their eigenvalues at the crossing are approximately

    -0.002827836114
     0.007649369477
     0.011861501885
     0.011861501885
     0.059903485447.

Thus one negative direction already exists before the D3 transverse pair is
considered.

## 3. Cubic coefficient

For the dual potential

    Phi(theta,g)
      = ||theta||^2/(2g)
        - log mean_j exp[(X theta)_j],

all derivatives of order >=3 are negative categorical cumulants.

At the crossing,

    D3_xxx
      = -0.0154444345753,

    D3_xyy
      = +0.0154444345753,

while

    D3_xxy ~ 0,
    D3_yyy ~ 0.

This is exactly the D3 tensor pattern

    D3 Phi
      = 6 c Re(z^3),

with

    z=x+i y

and

    boxed:
    c = -0.00257407242921.

Therefore the symmetry-allowed cubic invariant is actually present and
nonzero.

The crossing is not an accidentally cubic-free event.

## 4. Quartic coefficient after stable-mode elimination

Let H_s be the Hessian restricted to the five-dimensional noncritical
complement.

For a critical direction e define the stable cubic source

    t(e,e)
      = P_s D3Phi[e,e,.].

Lyapunov-Schmidt elimination gives

    D4_red[e,e,e,e]
      =
      D4Phi[e,e,e,e]
      -3 <t(e,e), H_s^{-1} t(e,e)>.

For e=e_x the numerical values are

    direct D4_xxxx
      = 0.0400152524039,

    reduced D4_xxxx
      = 1.42264653308.

D3 symmetry requires the reduced quartic to be radial at degree four.

Indeed,

    D4_red_xxyy
      = 0.474215511026
      = D4_red_xxxx/3

to numerical roundoff.

Thus the reduced quartic term is

    beta (x^2+y^2)^2

with

    boxed:
    beta = 0.0592769388783 > 0.

The positive quartic stabilizes large enough critical amplitude locally, but
the cubic controls the bifurcation type.

## 5. Parameter coefficient

MP7-039 certifies the total derivative of the critical 2x2 block determinant
along the base branch:

    d Delta45/dg
      ~ 0.00169588692471.

At the crossing the other eigenvalue of that block is

    lambda_other
      ~ 0.0118615018852.

Therefore the soft eigenvalue satisfies

    lambda_c(g)
      = ell (g-g*) + O((g-g*)^2),

with

    boxed:
    ell ~ 0.142974046721 > 0.

So the base branch has:
- two extra negative critical directions for g<g*;
- two positive critical directions for g>g*.

Including the pre-existing negative mode, its full local Hessian index changes

    3  ->  1

through the crossing.

## 6. D3 reduced potential

To leading nontrivial order,

    Phi_red
      =
      Phi_base
      + (ell mu/2) r^2
      + c r^3 cos(3 phi)
      + beta r^4
      + higher terms,

where

    mu=g-g*.

Because c != 0, stationarity gives three daughter rays on each side.

At leading order,

    r
      ~ |ell/(3c)| |mu|,

with

    boxed:
    |ell/(3c)| ~ 18.5146365863.

For mu>0 the daughters lie on

    cos(3phi)=+1,

and for mu<0 they lie on

    cos(3phi)=-1

in the chosen orientation.

Thus the three rays rotate by pi/3 when crossing g*.

## 7. Daughter stability

At a nonzero daughter branch, the two critical-plane Hessian eigenvalues are,
to leading order,

    lambda_radial   ~ -lambda_c,
    lambda_angular  ~  3 lambda_c.

They always have opposite signs.

Since one additional negative Hessian direction already exists in the
five-dimensional complement, every generic daughter has full Hessian index

    boxed:
    index = 2

on both sides of the crossing.

Therefore the D3 crossing does NOT create a new stable minimum.

It restructures the saddle network.

## 8. Direct 7D replay

At

    g=g*+1e-4,

three full seven-dimensional stationary roots were found at critical radii

    r ~ 0.00198635

on the three `cos(3phi)=+1` rays.

All three have full Hessian index 2.

At

    g=g*-1e-4,

three roots were found at

    r ~ 0.00174782

on the three opposite rays.

Again all have full Hessian index 2.

The leading asymptotic prediction is

    r ~ 0.00185146

at |g-g*|=1e-4.

The observed asymmetry is consistent with higher-order terms.

## 9. Orbit structure

The base state has stabilizer D3 of order six.

A generic daughter on one selected reflection ray retains only one reflection,
so its stabilizer has order two.

Therefore:
- three daughters surround each D3-symmetric base representative;
- under full D12 they form a twelve-member orbit.

This matches the expected symmetry-breaking count.

## 10. Energy splitting

Using the cubic normal form at a stationary daughter,

    Phi_daughter-Phi_base
      =
      [ell^3/(54 c^2)] mu^3
      + O(mu^4),

with

    ell^3/(54 c^2)
      ~ 8.16838770704.

Hence:
- for g>g*, daughter saddles lie above the base branch;
- for g<g*, they lie below the base branch.

This is an exchange in the saddle hierarchy, not a birth of a stable phase.

## 11. Physical/emergence interpretation boundary

The result gives a useful structural lesson for FIN:

    soft mode != automatically new stable state.

When the residual isotropy allows a nonzero cubic invariant, the soft
two-dimensional mode can generically create only saddle daughters near the
crossing.

Stable emergent daughter phases require either:
- symmetry to forbid the cubic;
- an accidental vanishing of its coefficient;
- or a different higher-order mechanism.

No physical particle, force or spacetime interpretation is asserted.

## 12. Next atom

### SYMMETRY-SELECTION-OF-STABLE-PHASES-67

Classify the lowest invariant degree for 2D critical irreps of the relevant
dihedral stabilizers.

Goal:
determine which symmetry classes generically permit cubic terms and therefore
saddle-only daughter bifurcations, and which force the first anisotropy to
quartic or higher order where stable daughter minima can occur.

This converts the MP7-039 counterexample into a general phase-selection rule.
