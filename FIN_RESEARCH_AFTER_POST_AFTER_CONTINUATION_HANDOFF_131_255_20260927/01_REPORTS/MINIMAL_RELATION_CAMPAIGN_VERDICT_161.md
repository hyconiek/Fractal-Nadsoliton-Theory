# MINIMAL-RELATION-CAMPAIGN-VERDICT-161
## Sparsity and persistence cannot both be obtained from the minimal exchangeable bond model without new structure

Date: 2026-09-26

Reports 158-160 test the smallest possible relational extension:

    identical FIN units
      +
    pair-specific binary bonds.

The result is restrictive.

## Independent Gibbs bonds

With a fixed pair energy:
- edge probability is O(1);
- mean degree is O(n);
- the graph is dense.

Finite mean degree requires

    epsilon_n ~ log n,

which is system-size-dependent tuning.

## MaxEnt bond dynamics

With the same microscopic refresh clock:
- bond memory decays as exp(-t);
- neighborhoods rewire on an O(1) time;
- long-time interactions are annealed/mean-field.

Making bonds persistent by hand adds a new kinetic ratio.

## Fixed valence

A degree constraint can enforce sparsity without log(n) tuning, but:
- introduces an unsourced integer d;
- does not produce spatial locality or dimension.

## Updated missing law

The next fundamental law must therefore do two things at once:

    1. select sparse pair relations;
    2. make those relations dynamically persistent on the relevant coarse
       time scale.

And it must do so without:
- coordinates inserted by hand;
- n-dependent retuning;
- a separately fitted relation clock;
- a separately chosen valence for every model.

## Next atom

`RELATIONAL-CONSERVATION-LAW-162`

Test whether an edge-count/valence constraint can arise as a conserved charge
or exact bookkeeping consequence of a reversible microscopic extension.

This is more promising than adding another arbitrary relational potential:
a conservation law could simultaneously control sparsity and persistence.
