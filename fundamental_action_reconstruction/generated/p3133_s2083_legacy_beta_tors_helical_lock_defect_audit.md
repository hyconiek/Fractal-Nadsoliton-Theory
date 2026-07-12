# P3133/S2083 Legacy beta_tors helical-lock defect audit

Status: `P3133_LEGACY_BETA_TORS_HELICAL_LOCK_DEFECT_GAP_CERTIFICATE`

## Question
Does the legacy beta_tors plus sinusoid kernel already export the P3132-required strict helical-lock defect D_HL?

## Fresh comparison

### alpha_geo amplitude
- legacy role: multiplicative numerator normalization
- strict successor: implicit/absorbed operational normalization in K_strict_gate
- D_HL-like nontranslation trace: `False`
- nonzero inversion-odd value: `False`
- support-representative coupling: `False`
- gap: amplitude can scale a kernel but cannot by itself choose a helical support origin or chirality polarity

### cos(omega*d+phi) sinusoid
- legacy role: oscillatory phase/resonance carrier
- strict successor: omega,phi refrozen by the QW2038-QW2050 spectral/intersection chain
- D_HL-like nontranslation trace: `True`
- nonzero inversion-odd value: `False`
- support-representative coupling: `False`
- gap: the sinusoid preserves a phase/winding trace, but cos is even and translation-orbit blind unless a separate phase-origin localizer is exported

### beta_tors linear denominator
- legacy role: 1/(1+beta_tors*d) torsion/topological path-summation damping
- strict successor: 1/(1+beta*d**eta) nonlinear damping/compression with beta=1, eta=1.8 in the working strict gate
- D_HL-like nontranslation trace: `False`
- nonzero inversion-odd value: `False`
- support-representative coupling: `False`
- gap: beta_tors is a scalar damping/torsion marker; current artifacts do not make it an inversion-odd, support-local helical-lock defect

### combined alpha_geo*cos/(1+beta_tors*d)
- legacy role: effective ontological bridge kernel
- strict successor: completed/enriched strict gate only under explicit completion-map evidence
- D_HL-like nontranslation trace: `True`
- nonzero inversion-odd value: `False`
- support-representative coupling: `False`
- gap: the combined formula can motivate D_HL, but it still lacks an exported odd sign value and an absolute support representative after the diagonal Z12 quotient

## Certificate
- legacy atoms tested: `4`
- accepted D_HL sources: `0`
- sinusoid retains phase trace: `True`
- beta_tors retains scalar torsion/damping trace: `True`
- bridge completion exported: `False`

## Decision
The user's suspicion is partly right: something important became under-visible in the move from the legacy formula to the strict operational gate. The legacy numerator keeps an explicit sinusoidal phase/resonance carrier and beta_tors keeps an explicit torsion-damping label, whereas the strict refreeze records omega, phi, beta, eta as operational gate parameters. However, this is not yet the missing P3132 D_HL source. The legacy formula provides motivation and a gap diagnosis, but current artifacts do not export beta_tors or the cosine phase as a nonzero inversion-odd, support-local helical-lock defect coupled to an absolute support representative.

## Recommendation
Build exactly one explicit candidate formula for D_HL from the legacy phase/torsion split, for example an oriented phase-gradient/torsion residual that is odd under inversion and is evaluated on a support-local representative. Then audit only its provenance, diagonal-translation behavior, sign polarity, and support coupling before any Zeta_OS/Gamma_SO retest. Do not start physical role transfer or claim bridge completion from beta_tors alone.

## Source availability note
The external DIAGRAMS_KERNEL_TRANSFORMATION.md path was not available in this container; the audit therefore uses repo-restored legacy-kernel packets and QW1729/QW2041 provenance records.
