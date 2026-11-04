# FIN: resonances on the critical epsilon^-2 time scale

Date: 2026-11-04. Source commit: `1dcba969e12fabbfd43695f800687f5c93076b75`.

**Outcome: partial analytical progress; actual critical-time damping remains open.**

The current AGENTS routing was checked against newer execution records.
The previously recommended RP004 full-input self-audit had already been
performed. CR selected ST's `NEXT_CRITICAL_GATE.md`, so neither the old
audit nor the Gamma calculation was repeated. The subsequent continuation
found the existing CD/VD packages and corrected the routing in
`../FIN_RESEARCH_STATUS_20261104.md`. The initial lookup was not exhaustive.

## Results

1. The existing causal memory functional, evaluated on a constant-amplitude
   carrier, has an explicit secular primitive plus a uniformly bounded
   remainder. Its cubic coefficient has real part `-Gamma/omega` in the
   full-cosine-amplitude convention. This is a frozen-response theorem,
   not damping of the actual nonlinear PDE trajectory. Initial transients
   and the nondecaying shape response are retained.
2. The shape and direct quartic coherent frequency shift agrees exactly
   with CM's canonical normal form. The continuum correction is separate.
   The factor nine between sources F and 3F and the factor four between
   full and half amplitudes are explicitly checked.
3. Exact actual canonical modal equations retain both bound modes,
   continuum coupling and the true initial slips. Canonical momenta are
   not identified with derivatives of transformed positions. The coherent
   phase correction can be order one on critical time.
4. The inherited C1 bound for the fifth-order energy remainder pays an
   O(epsilon^2) contribution to modal integral equations through time
   S epsilon^-2. This is smaller than the O(epsilon) modal amplitude.
5. A cumulative-defect stability estimate specifies the missing theorem:
   an o(epsilon) oscillatory integral bound would imply critical modal
   approximation. Existing absolute bounds give only O(1) and O(epsilon)
   for the cubic and mixed quartic energy contributions, respectively.

## Remaining obligation

`OPEN-CR` in section 5 of the proof is an actual, uniform cumulative
mixed-interaction defect estimate after subtracting the declared resonant
candidate. It must include the shape channel, all mixed terms and the
actual initial state. The candidate may require additional secular terms;
its completeness has not been proved.

A full local field/physical-velocity approximation also needs a separate
whole-continuum o(epsilon) bound. Modal comparison alone cannot supply it.
One route and one refinement were attempted: retain the complex signed
primitive instead of integrating the absolute forcing. The attempt stops
at the unpaid inequality without assuming a slowly varying PDE envelope.

## Verification and scope

The standard-library Python checker tests exact polynomial and harmonic
algebra, signs, normalization factors, exponent boundaries and negative
controls. It separately verifies source SHA-256 and pinned Git blob IDs.
Results are recorded in `verification.json`. These finite controls do not
replace analytical proofs or independent review.

No actual epsilon^-2 damping law, all-time scattering or physical light
production is established. The action, geometry, coefficients and
preparation remain supplied inputs to a conditional classical model.
F01, strict/legacy, selector and physical-interpretation boundaries remain
unchanged.

## Publication update

The original research package was prepared without publication under the
then-current execution scope. The user subsequently explicitly instructed
commit/push for completed work and English for all new repository content.
This English report replaces the unpublished Polish draft; scientific
claims and input hashes are unchanged. That later instruction supersedes
the old no-publication restriction for this session.
