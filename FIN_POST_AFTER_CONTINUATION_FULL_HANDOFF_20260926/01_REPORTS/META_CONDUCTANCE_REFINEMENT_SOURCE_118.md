# META-CONDUCTANCE-REFINEMENT-SOURCE-118
## The metastable conductance fixes the coarse block but cannot source the refinement fiber block

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Inputs:
- reports 54-57: exact reversible maximum-entropy heat bath;
- reports 110-116: metastable conductance/resistance geometry;
- ST231: exact two-child refinement;
- ST245: complete self-adjoint exact-refinement classification.

Status:
**NO-GO FOR UNIQUE FIBER-RATE SOURCE FROM COARSE CONDUCTANCE ALONE.**

## 1. Exact refinement classification

Let L be the coarse self-adjoint Laplacian and let

    R : H_coarse -> H_refined

be the isometric uniform-fiber embedding.

ST245 proves that every self-adjoint exact refinement satisfying

    L_tilde R = R L

has the reducing-space form

    boxed:
    L_tilde
      =
      R L R^*
      direct-sum
      B

on

    Ran(R) direct-sum Ran(R)^perp,

where B is an arbitrary self-adjoint operator on the fiber complement.

Positivity only requires

    B >= 0.

Therefore exact preservation of the complete coarse semigroup determines the
coarse block exactly but leaves the complement block free.

## 2. What the new metastable conductance supplies

The reversible heat-bath reduction supplies a concrete coarse Dirichlet form

    E_meta(f)
      =
      (1/2) sum_ij C_ij (f_i-f_j)^2,

and hence a concrete coarse Laplacian L_meta.

This is genuine new information.

But L_meta acts only on coarse observables.

Under the uniform-fiber lift, every coarse observable is constant inside a
fiber. Therefore all vertical/fiber edges contribute zero to its Dirichlet
energy.

So no measurement of the coarse conductance form can determine the action of B
on Ran(R)^perp.

## 3. Identifiability theorem

Take any two positive complement generators

    B1 != B2.

Define

    L1 = R L_meta R^* direct-sum B1,
    L2 = R L_meta R^* direct-sum B2.

Then for every coarse vector f and every t>=0,

    exp(-t L1) R f
      =
    exp(-t L2) R f
      =
    R exp(-t L_meta) f.

Thus all:
- coarse diffusion distances;
- coarse heat kernels;
- coarse conductances after projection;
- coarse metastable transition statistics

are identical.

Therefore the fiber generator is statistically invisible to the entire coarse
metastable process.

## 4. Consequence

The new conductance bridge of report 115 is real, but it does NOT by itself
remove the ST231/ST245 refinement nonuniqueness.

The obstruction is now precisely localized:

    coarse metastable dynamics
        determines
    R L R^*

but not

    B on Ran(R)^perp.

Any unique fiber-rate law must therefore add a principle that acts on genuinely
fiber-resolving observables.
