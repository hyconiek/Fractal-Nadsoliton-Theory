# BARRIER-COMMUNICATION-ULTRAMETRIC-111
## The leading communication height induces a 3×4 hierarchical geometry on the twelve minima

Date: 2026-09-26

Status:
- exact minimax statement on the asymptotic index-one saddle graph;
- interpretation as physical space is explicitly NOT made.

Define the leading communication cost between minima i,j as the smallest
possible maximal O(g) edge coefficient along a path in the pair-saddle graph.

Because d=3 edges connect each residue class modulo 3 and d=4 is the next
barrier class that joins those components, the minimax coefficient is

    c(i,j)=0                         if i=j,

    c(i,j)=alpha_3                  if i!=j and i=j mod 3,

    c(i,j)=alpha_4                  otherwise,

with

    alpha_3=0.275091465633409,
    alpha_4=0.281135945088269.

This function obeys the strong triangle inequality

    c(i,k)
      <= max[c(i,j),c(j,k)].

So the leading barrier hierarchy defines a two-level ultrametric on the twelve
localized minima.

The hierarchy is:

    level 0:
      12 individual minima;

    level 1:
      3 clusters of 4 states,
      linked internally by d=3 saddles;

    level 2:
      one 12-state component,
      linked by d=4 saddles.

This is not the original cyclic label metric.

For example labels separated by 6 are expensive as a DIRECT pair saddle, but
their minimal communication path uses two d=3 edges and therefore has the
cheapest communication height.

So the dynamical/energetic geometry emerging from the FIN landscape differs
from naive label distance.

This is potentially relevant to the broader emergence programme:

    underlying labels
      !=
    effective communication geometry.

But no identification with physical spatial distance is currently licensed.
