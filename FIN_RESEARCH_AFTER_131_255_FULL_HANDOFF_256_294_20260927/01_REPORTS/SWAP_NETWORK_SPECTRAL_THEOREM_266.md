# SWAP-NETWORK-SPECTRAL-THEOREM-266
## The full conservative color-exchange process has the same spectral gap as the underlying random walk

Date: 2026-09-27

Status:
mathematical upgrade using the established Aldous spectral-gap theorem for interchange processes, plus an elementary quotient/eigenfunction argument.

Consider any finite connected weighted graph G.

On each edge {i,j}, swap the endpoint records at the edge rate.

For fully labelled particles this is the interchange process.

A theorem of Caputo, Liggett and Richthammer proves:

    boxed:
    gap(interchange process)
      =
    gap(single-particle random walk)

on the same weighted graph.

## From labelled particles to three FIN colors

The Z3 SWAP model is a quotient of the labelled interchange process:
particle identities are forgotten while only three color classes are retained.

Taking a Markov quotient cannot introduce an eigenvalue slower than the full process, so

    gap(color process)
      >=
    gap(interchange)
      =
    gap(random walk).

Conversely, suppose at least one color a is present but does not occupy every site.

Let phi be a random-walk eigenfunction at the random-walk spectral gap.

Define the color-density observable

    F_a(eta)
      =
      sum_i
      phi(i) 1_{eta_i=a}.

A direct generator calculation gives the SAME eigenvalue.

Hence

    gap(color process)
      <=
    gap(random walk).

Combining:

    boxed:
    gap(Z3 SWAP color process)
      =
    gap(random walk)

in every nontrivial color-count sector.

## Cycle consequence

For the n-cycle with edge swap rate rho/2:

    boxed:
    gap
      =
      rho[
        1-cos(2 pi/n)
      ]

for ALL n and every nontrivial fixed-color-count sector.

Therefore report 241's n=3,6,9 exact diagonalizations were not merely small-size evidence.

They are instances of an all-size theorem.

At large n:

    gap
      ~
      2 pi^2 rho/n^2.

## General graph consequence

Once an incidence graph is sourced, the slowest conservative SWAP scale is no longer an independent many-body calculation.

It is determined by the graph random-walk Laplacian.

This sharply separates two questions:

1. FIN local transport law:
   SWAP and rho;

2. FIN geometry/incidence law:
   the weighted graph whose random-walk gap controls collective relaxation.

The second remains the main locality blocker.
