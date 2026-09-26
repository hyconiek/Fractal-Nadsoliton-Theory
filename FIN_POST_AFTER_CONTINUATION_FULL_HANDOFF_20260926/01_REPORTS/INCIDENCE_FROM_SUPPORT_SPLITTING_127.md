# INCIDENCE-FROM-SUPPORT-SPLITTING-127
## Ordinary support-addition daughters do not supply a stable binary refinement fiber

Date: 2026-09-26

Status:
- exact Morse-index obstruction using the certified Z2 support crossings;
- identifies one exceptional stable binary split: the uniform -> ±k6 pitchfork.

The proposed refinement interpretation was:

    parent state
      -> residual Z2
      -> two symmetry-related daughters

and then treat the two daughters as a binary fiber over the parent.

The already certified support-addition events are:

    S2/index1
      -> S3/index2,

and

    S4/index3
      -> S5/index4.

In both cases the broken-symmetry daughters gain one negative Hessian direction.

Therefore they are MORE unstable than the parent.

They cannot serve as a pair of metastable child states separated by the parent
as a transition saddle.

So the generic rule

    Z2 support breaking = physical refinement child pair

is false.

## The exceptional useful event

The first uniform-state k6 bifurcation is qualitatively different.

At

    g6=12/lambda6,

the uniform state undergoes a supercritical one-dimensional Z2 pitchfork.

For g>g6:

    uniform parent:
      index 1;

    +k6 daughter:
      index 0;

    -k6 daughter:
      index 0.

Thus here the parent really is the transition saddle between two locally stable
children.

This is the first existing FIN event with the correct Morse topology for a
binary metastable fiber.

## Consequence

FIN currently has:
- many symmetry-breaking support daughters;
- but only a subset can be interpreted as stable refinement children.

Morse topology supplies a necessary incidence filter:

    stable binary fiber candidate requires
      daughter index < parent index,

not merely two symmetry-related solutions.
