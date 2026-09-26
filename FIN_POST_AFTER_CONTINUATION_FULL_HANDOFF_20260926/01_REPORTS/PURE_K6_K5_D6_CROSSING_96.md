# PURE-K6–K5-D6-CROSSING-96
## Exact D6 transverse bifurcation replaces the earlier “near-k5 degeneracy” interpretation

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
- exact location of the two-dimensional k5 crossing on the pure-k6 branch;
- exact symmetry classification;
- fourth-order center-manifold coefficient from exact cumulant/Schur reduction;
- sixth-order angular coefficients from high-precision slaved center-manifold evaluation;
- daughter folds are numerical.

This report SUPERSEDES the earlier interpretation in reports 94/95 that the
d=2 branch entered an unresolved accidental near-coincidence of a fold and
transverse crossing.

## 1. Pure k6 branch

Write

    h_j = J6 (-1)^j.

The exact stationary equation is

    J6 = (g lambda6/12) tanh(J6).

The nonzero branch is born from the uniform state at

    g6 = 12/lambda6
       = 5.123427551398616.

At the later value

    g* = 12/lambda5
       = 5.220554796949195,

the nonzero pure-k6 solution has

    J6* ≈ 0.238932043126106,

or dual k6 coordinate

    theta6* ≈ 0.540822431626834.

## 2. Exact double k5 zero

Under the parity-biased pure-k6 law, both real k5 directions have zero mean and

    Var(k5c)=Var(k5s)=lambda5/12

independently of J6.

Therefore their Hessian eigenvalue is exactly

    lambda_k5(g)
      = 1/g-lambda5/12.

Hence BOTH directions vanish simultaneously and exactly at

    boxed:
    g*=12/lambda5
      =5.220554796949195.

This is not a numerical near-degeneracy.

The pure-k6 state is fixed by a D6 subgroup of D12.  The two k5 directions form
the standard real two-dimensional D6 representation.

So the correct local object is a D6-equivariant two-dimensional bifurcation.

## 3. Parent Morse index

At g=g* the full Hessian spectrum has:

    one negative direction,
    two zero k5 directions,
    four positive directions.

Thus the pure-k6 parent has index 1 immediately below the crossing.

## 4. Reduced normal form

Let

    z=r exp(i phi)

be the critical k5 coordinate.

D6 symmetry permits

    r^2,
    r^4,
    r^6,
    r^6 cos(6phi), ...

and forbids lower-order angular anisotropy.

After eliminating the five noncritical directions, the reduced potential has

    Phi_eff
      =
      Phi0
      +(lambda/2)r^2
      +a4 r^4
      +r^6[a6+c6 cos(6phi)]
      +O(r^8),

with

    lambda = 1/g-lambda5/12,

    D4 Phi_eff
      = -0.012172433674900,

so

    boxed:
    a4 = -0.000507184736454 <0.

High-precision center-manifold evaluation gives the two reflection-axis
sixth-order coefficients

    cos(6phi)=+1:
      a6+c6 ≈ 0.158452059191929,

    cos(6phi)=-1:
      a6-c6 ≈ -0.001205067803913.

Equivalently,

    a6 ≈ 0.078623495694008,
    c6 ≈ 0.079828563497921 >0.

## 5. Two inequivalent daughter families

Because a4<0, both D6 reflection-axis families are SUBCRITICAL:
small daughters exist for g<g*.

But the sixth-order terms make them very different.

### Family A: cos(6phi)=+1

The sixth-order coefficient is strongly positive:

    +0.158452...

The critical plane contributes two negative directions near birth, so together
with the parent's pre-existing negative direction the daughter has index 3.

The sixth-order truncation already predicts a nearby fold.

It predicts

    r_fold≈0.032664312,
    g_fold≈5.220525300039.

The full stationary equations give

    boxed:
    g_fold,A≈5.220523246995448,

only about

    3.15499537e-05

below the exact D6 crossing.

After this fold the branch has index 2 and continues to the strict large-g
three-label support

    S={0,2,10},

with weights

    (0.3738292547693,
     0.31308537261535,
     0.31308537261535).

### Family B: cos(6phi)=-1

Here even the sixth-order coefficient remains negative:

    -0.0012050678...

So sixth order cannot turn the subcritical branch back.

The first fold is much farther away:

    boxed:
    g_fold,B≈5.200804729029990.

Across it the daughter changes

    index 2 -> index 1,

and the return sheet tends to the large-g d=2 pair support

    S={0,2},
    p=(1/2,1/2).

## 6. Orbit counting

The D6 parent has stabilizer order 12 and full D12 orbit size 2.

Each reflection-axis daughter has stabilizer order 2, hence full orbit size 12.

The two signs

    cos(6phi)=+1
and
    cos(6phi)=-1

are two inequivalent D6 reflection classes.

Thus one two-state parent orbit connects locally to TWO distinct twelve-state
daughter orbits.

## 7. Correction to the previous d=2 picture

The d=2 pair component should now be read as

    pure k6 parent, index 1
       |
       | exact D6 k5 crossing
       | g=12/lambda5
       v
    subcritical cos(6phi)=-1 daughter, index 2
       |
       | fold g≈5.20080472903
       v
    d=2 pair sheet, index 1
       |
       -> large-g support {0,2}.

There is no need to postulate an unresolved accidental fold/crossing collision
near g5.

The apparent double softness is the exact two-dimensional D6 critical space.
