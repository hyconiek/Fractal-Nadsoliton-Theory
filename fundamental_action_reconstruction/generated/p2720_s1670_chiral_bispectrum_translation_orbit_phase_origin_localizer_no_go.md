# P2720/S1670 chiral-bispectrum translation-orbit phase-origin localizer no-go

Status: `P2720_TRANSLATION_ORBIT_PHASE_ORIGIN_LOCALIZER_NO_GO`

## Orbit result
- accepted_localizer_count=0
- marker_constant_on_each_translation_orbit=True
The P2718 bispectrum sign separates orientation but is constant across the full Z12 translation orbit for each orientation.  Any localizer using only this orbit-quotiented marker data cannot select one source/phase origin without importing an external source label or gauge convention.

## Decision
P2720 attacks the P2719 localizer obligation directly.  The exact P2718 marker keeps its orientation sign, but for each orientation it is constant on the full 12-element translation orbit.  Thus it cannot select a non-premise source/phase origin without importing an external source label or translation-gauge convention.

## Next honest step
Do not repeat the translation-orbit localizer attempt for Im(B_{1,5}).  A next admissible move must supply a genuinely new origin-sensitive strict invariant/source law that is not translation-gauge/source-label data, or an independent exported torsor-coupling theorem; otherwise preserve the P2697-P2720 no-new-live-frontier certificate.
