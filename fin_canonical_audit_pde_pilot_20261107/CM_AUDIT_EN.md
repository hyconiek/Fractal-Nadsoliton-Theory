# CA: final targeted same-author audit of CM

Date: 2026-11-07. Parent e2c8a788252d3cedc21fd8f1939d59e00cb8350e.
Disposition: no substantive error found at the stated dependency boundary.
AW remains a conditional theorem requiring independent review. Routine
same-author audit chaining ends here under the user's revised priority.

## Exact algebra reconstructed

The companion canonical.py independently constructs the 19-dimensional
even-pair quartic space, its oscillator derivative D and angle projector
Pi. It solves the three blocks of (D+Pi)G4=Q4b-Pi Q4b rather than copying
an existing list of G4 coefficients. The checker confirms:

- D G3=E3b, with OLD=Phi3(NEW) cancelling the bound cubic term.
- Pi^2=Pi, Pi D=D Pi=0, and the full determinant
  2^20 w^3 3^3 (w-3)^2, nonzero for w in [5/4,7/4].
- D G4=Q4b-Z4, Pi G4=0 and all three displayed coherent coefficients.
- Momentum reversal and psi-sign symmetries of both generators.
- The exact exceptional fixed point at w=3/2 and preservation of the
  zero-shape subspace by Phi4, without declaring the physical forced shape
  identically zero.
- The DC and second-harmonic identities in the physical shape reconstruction.

The convention OLD=Phi3 composed Phi4(NEW) implies that numerical physical
pullback first applies Phi3 inverse and then Phi4 inverse. The auxiliary
flow refinement at the largest pilot amplitude changes coordinates by
less than 1.6e-14. This is a finite-dimensional implementation check, not
a guarantee for arbitrary amplitudes.

## Energy-space remainder and unchanged continuum

The energy space has H1 fields and L2 momenta. The nonlinear position
energies of degrees three, four and six are smooth there by the 1D H1
algebra and bounded kink coefficients. Each generator depends only on
four bounded modal coordinates. Its vector field is finite rank and
smooth, with quadratic or cubic growth. A sufficiently small ball admits
the forward and inverse auxiliary flows for unit auxiliary time, uniformly
on the compact parameter band. No unbounded continuum derivative is
introduced by the Poisson brackets with these generators.

Taylor's integral remainder for the composed energy therefore gives
R5=O(||state||^5) and dR5=O(||state||^4) in the dual energy norm.
This is the bound MA uses after contraction with finite modal action
gradients. It is not a claim that an arbitrary continuum Hamiltonian
vector field is bounded in the same norm.

The finite-dimensional Hamiltonian flows preserve the modal symplectic
form and extend by identity on the symplectically orthogonal continuum.
The whole map is canonical. Degree counting retains E3m and
Q4m=E4m+{E3m,G3}; the G4 interaction with E3m starts at degree five.
In particular no radiative mixed term is removed. A finite energy jet
alone would not prove damping; the later response/work estimates remain
necessary inherited arguments.

## Preparation and verdict

Zero momenta are preserved by both auxiliary flows because the generators
are odd under simultaneous momentum reversal. Psi-sign symmetry makes
Z0 odd and C0 even in epsilon. This gives Z0=epsilon+O(epsilon^3),
C0=-d A epsilon^2+O(epsilon^4), with exactly zero new momenta and continuum.
At w=3/2 the cubic generator fixes the raw preparation and the quartic
flow preserves zero shape, so C0=R0=0 exactly. The physical shape is still
forced during PDE evolution. Canonical momenta must not be replaced by
canonical position derivatives.

The reconstruction found no substantive CM error. Fifty checks passed
(33 input hashes and 17 algebra/flow controls). This closes the requested
same-author CM audit, not independent review of AW or the full historical
chain. Further audits are triggered by a concrete failure or new result,
not scheduled as a continuing default research program.
