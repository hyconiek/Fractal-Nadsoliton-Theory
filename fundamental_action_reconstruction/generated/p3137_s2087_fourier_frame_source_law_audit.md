# P3137/S2087 Fourier-frame source law F_DHL audit

Status: `P3137_FOURIER_FRAME_SOURCE_LAW_F_DHL_BOUNDED_NO_GO`

## Constructed object
- `primitive character/frame plus phase-zero reference for J_DHL`
- `primitive characters U(12)={1,5,7,11}, primitive pairs {1,11}/{5,7}, and 12 phase-zero cells`

## Repo backscan
- P2992: no nonpremise frequency/source localizer for Z12 Fourier characters.
- P2994: exact receiver matrix exists, but no atom-specific source-coupling theorem is exported.
- P3039: chi_i sine projector is real and inversion-odd, but phase origin remains chart-relative.
- P2869/P2870: Aut-character idempotents represent endpoints but import projector/polarity and do not intrinsically select the needed character.

## Finite certificate
- profiles tested: `120`
- primitive character orbit size: `4`
- primitive pair orbit size: `2`
- profiles with primitive active pair: `48`
- profiles with nonprimitive active pair: `72`
- phase-zero cells: `12`
- source candidates tested: `9`
- candidates passing all gates: `0`
- accepted F_DHL sources: `0`

## Decision
P3137 constructs the missing F_DHL frame-source candidate space and tests it against repo-backed Fourier/character/localizer results. The finite receiver side is real: every P3134/P3136 profile has an active ±k character pair, and 48 rows live on primitive pairs. But the primitive characters form one Aut(Z12) orbit and the primitive pairs form one orbit; no invariant frame data chooses a unique primitive character/pair. Phase-zero has 12 translated cells. Across nine source candidates, zero simultaneously select a primitive character/frame, select phase-zero, and remain import-free. Thus F_DHL is not exported on current artifacts.

## Recommendation
Do not continue Fourier-frame variants unless a genuinely new strict frame-breaking source is supplied. The next admissible move is either one new non-Fourier joint source object for (r,lambda), or a no-new-live-frontier reconciliation for the D_HL lane after P3133-P3137. If continuing constructively, the object must not be a Fourier receiver, character projector, phase gauge, lexicographic label rule, or prior chi_i/localizer replay; it must compute a nonconventional support origin and polarity in one law.
