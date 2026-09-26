# GENERATOR-SET-DIMENSION-NONUNIQUENESS-246
## A transformation group does not determine spatial dimension until the elementary generator set is sourced

Date: 2026-09-26

Status:
exact finite-group counterexample using the already established Z12 carrier.

The internal label set is

    Z12
      ≅
    Z3 × Z4.

One might be tempted to say that the CRT product automatically supplies two
dimensions.

It does not.

## 1. One elementary generator

Take the bijection

    j -> j+1.

Its Cayley graph with ±1 edges is

    C12.

It has degree two and one elementary transformation direction.

The same is true for the generator

    j -> j+5,

because gcd(5,12)=1.

In CRT coordinates the step +5 changes BOTH factors:

    5 -> (2,1),

yet the graph is still only one 12-cycle.

So:

    "acts on two CRT coordinates"

does NOT mean

    "two spatial dimensions".

## 2. Two elementary generators

Instead choose the commuting steps

    +3 -> (0,-1)
    +4 -> (+1,0).

Using both ±3 and ±4 edges gives exactly

    C3 square C4,

with degree four and two elementary product directions.

## 3. Same vertices, same information preservation, different geometry

Both constructions:
- use the same 12 vertices;
- use bijective information-preserving transformations;
- are transitive/connected;
- require no external coordinates.

Yet one is a cycle and the other a two-factor torus graph.

Therefore:

    boxed:
    group + orbit structure
      is insufficient to determine dimension.

One must additionally source WHICH transformations are elementary/local.

This is the transformation-language version of the old incidence problem.
