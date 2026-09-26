# RELATIONAL-CONSERVATION-LAW-162
## Conserving total edge number can enforce sparsity, but does not select locality

Date: 2026-09-26

Status:
exact fixed-charge ensemble theorem.

Let

    P=n(n-1)/2

be the number of possible unordered unit pairs.

Introduce binary relations b_e and impose the exact conserved charge

    boxed:
    M=sum_e b_e.

Condition on a fixed M.

With no further pair structure, maximum entropy gives the uniform ensemble over
all graphs with exactly M edges: G(n,M).

## 1. Sparsity from an extensive conserved charge

Every pair is equivalent, so

    P(edge e active)
      =
      M/P.

Hence

    E[degree]
      =
      (n-1) M/P
      =
      2M/n.

If the conserved charge scales as

    M=(d/2)n,

then

    boxed:
    E[degree]=d

independently of n.

This avoids the log(n) chemical-potential tuning of report 158.

So a conserved extensive relation charge is a genuine mechanism for sparsity.

## 2. But the value of the charge is not selected

The theory must still explain why the physical sector has

    M proportional to n

and which density

    d=2M/n

is realized.

The conservation law preserves M.

It does not choose M.

Thus the missing parameter has moved from:
- an edge chemical potential

to:
- a superselection / initial-charge density.

## 3. No locality in the MaxEnt fixed-M ensemble

At fixed M, all unordered pairs remain exchangeable.

The mean adjacency is

    E[W]
      =
      p(J-I),

with

    p=M/P.

The mean graph Laplacian is therefore

    E[L]
      =
      p(nI-J),

whose nonzero spectrum is completely degenerate:

    boxed:
    n p
    repeated n-1 times.

For M=dn/2,

    n p
      =
      nd/(n-1)
      ->
      d.

So the annealed large-n theory is sparse in degree but still has no preferred
neighbor shells or directions.

## 4. Result

Global relation-number conservation solves:

    dense vs sparse

but not:

    which pairs are neighbors?

It generates a sparse exchangeable network ensemble, not physical local
geometry.
