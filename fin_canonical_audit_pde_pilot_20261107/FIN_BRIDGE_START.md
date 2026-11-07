# Starting the FIN-to-field bridge: an explicit identifiability test

This is a first mathematical obstruction and matching specification, not
a derivation of the RP action. It follows the user's priority change.
The strict-core spectral generator and the supplied RP model must remain
distinct. No lower informational substrate is introduced.

The current target potential is

    V(phi,psi)=(phi^2-1)^2/2 +(a+6phi^2)psi^2/2
                                  +g4 psi^4/4+g6 psi^6/6.

Its variational field equations with the canonical wave kinetic term
are precisely the supplied RP equations. But neither the kinetic choice
nor the nonlinear coefficients are selected merely by reusing a spectral
operator.

## A concrete non-identifiability witness

Fix any g4,g6 and change them to g4+delta4,g6+delta6. The potential changes
by delta4 psi^4/4+delta6 psi^6/6. Both its first and second field derivatives
vanish at every point with psi=0. Thus the static kink, its energy, and
both entire linearized kink operators remain unchanged for this family,
while the nonlinear equations change. This is stronger than merely
matching a few eigenvalues: the full quadratic fluctuation operators
coincide. The symbolic bridge checker verifies these identities.

Therefore matching the strict FIN quadratic spectral information to the
RP kink operators, even if achieved, cannot by itself identify g4,g6.
This is a limitation of that proposed spectral matching route, not a
no-go theorem for every possible nonlinear FIN construction. Additional
nonlinear information already derivable from the strict core could in
principle break the degeneracy; it must be explicitly supplied.

## Concrete next matching obligations

A candidate bridge must export a state/observable map and admissible
history before the RP equation can be read as its effective evolution.
For a canonical local potential, required derivatives include

    V_phiphi(1,0)=4, V_psipsi(1,0)=a+6,
    V_phiphipsi psi=12,
    V_psi psi psi psi(1,0)=6g4,
    V_psi psi psi psi psi psi(1,0)=120g6.

The mixed derivative line means two phi and two psi derivatives. These
are targets to derive with an error estimate, not conditions to impose
by fitting the desired equation. One must also derive the wave kinetic
term, propagation coordinate/locality and preparation map. The fixed
12-cycle by itself does not supply a continuum refinement family or a
physical clock; those cannot be silently added as strict-core consequences.

Next search the actual strict nonlinear compression/dynamics definitions
for a noncircular provider of the fourth- and sixth-order response data.
If no provider exists, record the precise additional premise required for
that candidate and keep it explicitly non-strict. Do not make another
empty general source inventory or rename the supplied PDE as its own origin.

No measurement map, SI calibration, quantum field identification or
empirical FIN confirmation follows from the numerical pilot.
