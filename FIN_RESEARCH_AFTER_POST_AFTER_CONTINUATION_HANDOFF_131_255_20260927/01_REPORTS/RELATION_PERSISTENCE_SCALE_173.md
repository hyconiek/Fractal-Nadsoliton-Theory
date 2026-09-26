# RELATION-PERSISTENCE-SCALE-173
## Persistent geometry requires a relational barrier tied to the slow coarse scale

Date: 2026-09-26

Status:
scale-separation requirement.

The first effective three-sector FIN dynamics becomes slow with increasing
internal copy number N.

Reports 141-144 find a coarse transition rate that decreases rapidly with N,
while the coupled projection-memory sector stays O(1).

Now consider a dynamical edge state.

If an edge activation/deactivation event changes the relational action by only

    Delta E_rel = O(1)

and uses the common microscopic clock, detailed-balance transition rates are
also O(exp[-O(1)]).

Therefore the relation lifetime remains O(1) in the same dimensionless clock.

But the effective unit-transition time grows rapidly with N.

Hence eventually

    tau_relation
      <<
    tau_unit_transition.

The graph becomes annealed many times before a unit changes coarse state.

## Requirement for persistent locality

A neighborhood that acts as quasi-static geometry for the slow unit dynamics
needs at least

    tau_relation
      >=
    tau_unit_transition

over the intended observational regime.

This requires one of:

1. an exact conservation of pair identity;
2. a relational barrier growing with N;
3. a collective relation metastability whose barrier is derived from the same
   microscopic law.

Choosing a separate small relation refresh rate is not acceptable; it inserts
a new kinetic ratio.

## Important distinction

There are now TWO size variables that must never be conflated:

    N
      internal copy number of one FIN unit;

    n
      number of simultaneous effective units.

Sparsity is mainly an n-scaling problem.

Persistence relative to the already derived three-state dynamics is mainly an
N-scaling problem.

A successful multicell law has to solve both at once.
