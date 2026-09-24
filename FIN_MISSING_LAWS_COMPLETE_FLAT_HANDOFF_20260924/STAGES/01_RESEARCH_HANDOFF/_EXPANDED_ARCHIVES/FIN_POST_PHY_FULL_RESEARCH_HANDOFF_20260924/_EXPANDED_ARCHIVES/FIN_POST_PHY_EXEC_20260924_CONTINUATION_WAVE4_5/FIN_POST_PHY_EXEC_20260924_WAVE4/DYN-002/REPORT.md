# DYN-002 — Coarse kinetic coefficients

Status: **STOP_NO_FIXED_MICROSCOPIC_UPDATE_LAW__KINETIC_COEFFICIENTS_UNSOURCED**.

The task is conditional on one fixed reversible local jump law. No such new law is supplied by the current programme, and DYN-001 shows that the state-dependent rewiring version additionally needs an explicit collision/update regularizer before a uniform coarse generator can be justified.

The obstruction is exact already at the simplest level: if `Q` is any reversible local Markov generator with stationary law pi, then `c Q` for arbitrary `c>0` has the same pi, detailed balance, support and equilibrium geometry while rescaling every kinetic coefficient and continuous time by c. Static FIN and PHA-001 therefore cannot determine mobility or seconds. Different accepted kinetic categories (gradient, inertial, complex) remain inequivalent as the predecessor audit established.

Because the required microscopic update primitive is absent, fitting a mobility to recover a desired coarse law would violate the stop rule. Outcome: **STOP / NEW TEMPORAL ATOM REQUIRED**. Once a law plus collision handling is independently frozen, generator expansion can be run as a new conditional campaign.
