# Z3-MARKOV-COUPLING-SIMPLEX-236
## Exact one-unit dynamics does not determine the joint dynamics of two effective Z3 units

Date: 2026-09-26

Input:
the established effective three-state generator

    Q3
      =
      k *
      [ -2  1  1
         1 -2  1
         1  1 -2 ].

Status:
exact classification inside the maximally symmetric translation-invariant
two-unit coupling grammar.

Let the two effective states be

    x,y in Z3.

Require:
1. each marginal generator is exactly Q3;
2. unit exchange symmetry x<->y;
3. Z3 reflection symmetry;
4. translation invariance in state space;
5. uniform detailed balance.

All eight nonzero increment vectors on Z3^2 fall into three symmetry classes.

### A: one-coordinate jumps

    (+1,0),(-1,0),(0,+1),(0,-1)

each at rate A.

### B: same-direction simultaneous jumps

    (+1,+1),(-1,-1)

each at rate B.

### C: opposite-direction simultaneous jumps

    (+1,-1),(-1,+1)

each at rate C.

Marginal consistency gives exactly one equation:

    boxed:
    A+B+C=k,

with

    A,B,C >=0.

Therefore the admissible joint laws form a TWO-DIMENSIONAL simplex.

## Collective eigenmodes

For the relative character

    chi_rel(x,y)=exp[2 pi i(x-y)/3],

the eigenvalue is

    boxed:
    lambda_rel=-6A-3C.

For the sum/common character

    chi_sum(x,y)=exp[2 pi i(x+y)/3],

    boxed:
    lambda_sum=-6A-3B.

Yet either single-unit character always has

    lambda_single=-3k.

So all admissible couplings have exactly the same one-unit process while their
collective physics differs.

## Extremes

### Independent continuous-time units

    A=k,
    B=C=0.

Then

    lambda_rel=lambda_sum=-6k.

### Perfect same-direction synchrony

    B=k,
    A=C=0.

Then

    lambda_rel=0.

The difference x-y is exactly conserved.

### Perfect opposite-direction synchrony

    C=k,
    A=B=0.

Then

    lambda_sum=0.

The sum x+y is exactly conserved.

## No-go

Even after:
- deriving Q3 from the microscopic FIN process;
- fixing the clock;
- imposing maximal state and unit symmetries;
- imposing detailed balance;

the two-unit collective law is still not unique.

A new relational condition must select a point (A,B,C) in this simplex.

This is the effective-Z3 analogue of the earlier quadratic eta/kappa
nonuniqueness, now expressed directly at the accepted coarse level.
