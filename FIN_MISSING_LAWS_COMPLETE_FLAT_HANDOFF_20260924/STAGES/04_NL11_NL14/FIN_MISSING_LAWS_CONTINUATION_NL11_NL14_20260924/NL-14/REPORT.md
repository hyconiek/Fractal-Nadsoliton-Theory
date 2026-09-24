# NL-14 — reversible moving-wall gate

Status: **STOP_MOVING_WALL_NOT_ACTIVATED_KINETIC_SOURCE_STILL_MISSING**

NL-13 is enough to justify a **static** radial-wall research branch under a
declared stiffness.  It is not enough to launch moving-wall numerics.

The missing pieces are still structural:

- NL-11 leaves four D12-compatible internal constitutive weights;
- NL-12 does not provide the Markov/refinement theorem needed to turn Fisher
  uniqueness into an internal FIN law;
- the earlier P3077/P3078 audits still do not export an intrinsic momentum
  variable, nondegenerate symplectic form, kinetic normalization and clock.

If we simply impose `M s_tt = G s_xx - grad Phi`, choose M and G to be
proportional and then Lorentz-boost NL-13, the resulting moving wall would be a
property of the **imported wave equation**, not a derivation from FIN.

Therefore the correct falsification-first action is STOP.
