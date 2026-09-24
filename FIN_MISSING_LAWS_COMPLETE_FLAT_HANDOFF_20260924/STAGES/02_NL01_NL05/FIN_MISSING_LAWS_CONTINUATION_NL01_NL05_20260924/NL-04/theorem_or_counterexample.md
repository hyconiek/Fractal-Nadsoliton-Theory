# Nonuniqueness theorem

Fix a refinement-compatible spatial mass M_s and stiffness K_s and a d>=2
internal deformation vector r. For every positive definite dxd matrix G,

H_G = 1/2 dot(r)^T (M_s tensor G) dot(r)
      + 1/2 r^T (K_s tensor H) r

defines a local, time-reversal invariant, energy-conserving second-order model.
Spatial subdivision acts only on M_s,K_s and therefore leaves G free. Direct
composition is block diagonal and also leaves G free. Two non-proportional G
matrices generically change ratios of normal-mode frequencies, so they are not
related by a single clock rescaling.

Therefore spatial refinement + composition + conservation do not source a
unique reversible kinetic metric when the internal fiber has dimension >1.
