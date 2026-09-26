# COMMUNICATION-BARRIER-REFINEMENT-NOGO-114
## The metastable barrier ultrametric cannot itself be the additive FIN refinement length

Date: 2026-09-26

Status:
- exact algebraic no-go;
- compares the new communication-height composition law with the already
  accepted arbitrary-split static refinement law.

## 1. Two different composition algebras

The large-g metastable communication height composes along a path by minimax:

    B(path)=max_e B_e,

and between endpoints by minimizing this maximum over paths.

By contrast, the existing FIN edge-refinement theorem uses a quadratic edge
coefficient c with series composition

    c_eff = c1 c2/(c1+c2),

or equivalently resistance

    R=1/c

with

    R_eff=R1+R2.

So the two natural candidate "length" variables obey:

    barrier:
      B_series=max(B1,B2);

    refinement resistance:
      R_series=R1+R2.

## 2. No scalar reparameterization repairs this

Suppose a positive scalar map f converted barrier height into an additive
refinement length:

    f(max(x,y))=f(x)+f(y)

for all x,y>0.

Set y=x. Then

    f(x)=2f(x),

hence

    f(x)=0

for every x.

Therefore there is no nontrivial positive scalar transformation of the barrier
ultrametric that turns minimax composition into arbitrary-split additive
length.

This is an exact no-go.

## 3. Consequence

The new communication ultrametric of reports 110-113 must NOT be identified
directly with the earlier FIN metric length.

In particular:

    barrier hierarchy != additive spatial metric.

The communication hierarchy remains a valid dynamical/energetic geometry, but
it belongs to a different composition category.

## 4. Why this is scientifically useful

The no-go prevents an attractive but incorrect shortcut:

    "FIN produced a barrier ultrametric, therefore FIN has produced the same
    geometry selected by its refinement theorem."

That inference is false.

A different object must be sought if the metastable dynamics is to connect to
the arbitrary-split refinement law.
