# P3136/S2086 Fourier-phase J_DHL candidate audit

Status: `P3136_FOURIER_PHASE_J_DHL_CANDIDATE_CHART_RELATIVE_POSITIVE_STRICT_SOURCE_NO_GO`

## Constructed object
`J reads C_k(D)=sum_x D(x)*exp(-2*pi*i*k*x/12); phase(C_k) reconstructs r modulo gcd(k,12), and sign convention reconstructs lambda when a Fourier frame is fixed.`

## Repo backscan
- P2992/P2994 Fourier-character source-atom lanes already warn that exact character receivers do not provide strict provenance or atom-specific source coupling.
- P3039 already tested chi_i=sin(2*pi*i/12): the torsor is real and inversion-odd, but translations move the phase origin and current receivers choose chart-relative representatives.
- P2869 warns that Aut-character Fourier idempotents can represent endpoints only by importing projector/polarity coefficients.

## Finite certificate
- candidate profiles tested: `120`
- nonzero Fourier coefficients: `120`
- unit-mode profiles: `48`
- nonunit-mode profiles: `72`
- chart-relative full pair recoveries: `0`
- aliased nonunit recoveries: `120`
- translation phase-rotation witnesses: `120`
- accepted import-free J_DHL sources: `0`
- mode alias counts: `{'1': [2], '2': [4], '3': [6], '4': [4], '5': [2]}`

## Decision
P3136 constructs an actual formula-level J_DHL candidate rather than another matrix label. The Fourier phase extractor is mathematically strong: all 120 P3134 profiles have nonzero Fourier coefficient at their generating mode. It also exposes a sharper obstruction than expected: even the primitive unit modes k=1 and k=5 have a half-period ambiguity (r,lambda) ~ (r+6,-lambda), while nonunit modes have larger aliases. Thus no row recovers an import-free full pair. The extraction assumes a labelled Fourier frame, a chosen mode/character, and a polarity convention; translation rotates the coefficient phase in all 120 rows. P3136 therefore proves a positive receiver theorem and a stricter joint-source no-go boundary at the same time.

## Recommendation
Do not repeat Fourier receiver extraction. The next admissible proof-grade object is one strict Fourier-frame/source law F_DHL that selects a primitive mode/character and phase-zero reference without importing chart labels. Its audit should be a finite frame-source obstruction/witness test against prior Fourier-character no-go lanes P2992/P2994, the chi_i localizer boundary P3039, and the Aut-character idempotent warning P2869. Without such F_DHL, preserve the P3135-P3136 conditional-positive/no-strict-source certificate.
