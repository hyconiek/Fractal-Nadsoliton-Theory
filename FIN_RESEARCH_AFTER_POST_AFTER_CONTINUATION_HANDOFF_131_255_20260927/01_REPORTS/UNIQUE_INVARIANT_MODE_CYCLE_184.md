# UNIQUE-INVARIANT-MODE-CYCLE-184
## A single invariant mode is equivalent to one transformation cycle

Date: 2026-09-26

Status:
exact spectral theorem for a finite permutation.

Let P be the permutation matrix of the information-preserving transformation T.

Define

    L_T
      =
      (I-P)^*(I-P)
      =
      2I-P-P^*.

This is the undirected transformation-graph Laplacian.

## 1. Kernel counts transformation orbits

A vector f lies in ker(I-P) iff

    f(Tx)=f(x)

for every x.

So f is constant on each permutation orbit.

Therefore

    boxed:
    dim ker(I-P)
      =
    number of cycles of T.

The same holds for L_T.

Hence:

    boxed:
    dim ker L_T=1
      iff
    T has one orbit
      iff
    the transformation graph is one cycle C_n.

## 2. Why this is useful

Connectedness can therefore be phrased without inserting a graph target:

    unique invariant mode.

This is structurally similar to a connected Markov/Laplacian system having one
constant zero mode.

So the conditional two-port architecture can be stated as:

    exact information preservation
      -> permutation P;

    unique invariant mode
      -> one cycle;

    L_T=(I-P)^*(I-P)
      -> cycle Laplacian.

## 3. Remaining source question

FIN already has unique zero modes in its accepted connected Laplacian objects.

But transferring “one zero mode” from the internal strict operator to a
multicell transformation P requires a typed composition theorem.

The spectral similarity alone does not authorize the role transfer.
