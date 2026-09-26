# TYPED-CONDUCTANCE-INCIDENCE-BOUNDARY-156
## Metastable conductance cannot be reused as inter-unit incidence without a new typed map

Date: 2026-09-26

Status:
typed no-go / category distinction.

FIN now contains a dynamically generated conductance matrix among metastable
states of ONE unit.

Those nodes are mutually exclusive alternatives of the same microscopic
system.

A multicell incidence graph would instead connect SIMULTANEOUS subsystems.

These are different object types:

    metastable-state node
      = alternative macrostate of one unit;

    physical-unit node
      = simultaneously existing subsystem.

The conductance

    C_ab

between alternative states a,b is a transition statistic.

It does not answer whether two simultaneous units x,y are incident.

A map

    (state-transition conductance)
      ->
    (inter-unit adjacency)

would require an additional typing rule identifying:
- which alternative-state relation corresponds to a bond between subsystems;
- how bonds compose when both subsystems change state;
- whether bond variables have independent storage/memory.

No such map is derived in the current repo.

Therefore report 115's conductance/refinement algebra bridge remains valid but
cannot be promoted to an incidence-source theorem.
