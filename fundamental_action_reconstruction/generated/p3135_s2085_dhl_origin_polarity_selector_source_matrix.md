# P3135/S2085 D_HL joint origin/polarity selector-source matrix

Status: `P3135_DHL_JOINT_ORIGIN_POLARITY_SELECTOR_SOURCE_MATRIX_BOUNDED_NO_GO`

## Repo backscan
- P3130/P3131/P3132 already block translation-origin, origin-twist, and relative-helix section routes (do not replay those as new D_HL sources).
- P2718-P2721 already establish nonzero chiral sign but missing origin/localizer and polarity law (candidate lambda support without r remains insufficient).
- P3059/P3064 already isolate polarity-source atoms and paired sign obstructions (lambda selection requires exported source atoms, not sign-even scoring).

## Finite orbit certificate
- pair space: `Z12 origins r crossed with lambda in {-1,+1}`
- pair count: `24`
- orbit count: `1`
- largest orbit: `24`
- invariant unique pair selectors: `0`

## Source matrix certificate
- candidate sources tested: `10`
- sources selecting r: `3`
- sources selecting lambda: `5`
- import-free sources: `4`
- sources passing all three gates: `0`
- accepted joint r/lambda sources: `0`

## Decision
P3135 builds the joint (r,lambda) selector-source matrix requested by P3134 and backscans the repo before testing. The finite action on Z12 origins crossed with the D_HL polarity has one 24-element orbit under translations, Aut(Z12) units, and lambda pairing, so invariant data alone cannot uniquely select a pair. Ten repo-supported source candidates are tested. Some conditionally select r, some conditionally select lambda, and some are import-free as scoped no-go objects, but zero candidates satisfy all three gates simultaneously. Thus the constructed D_HL remains a real local object with a precise missing source: a single strict law must choose both support origin and polarity together.

## Recommendation
Do not replay separate origin-only or polarity-only lanes. The next admissible proof-grade object is exactly one joint source law J_DHL with formula-level content, e.g. a nadsoliton-internal functional J_DHL(support, field) -> (r,lambda), whose first audit must test whether it breaks the single 24-element orbit without importing chart order, selector premises, apparatus, observed light, Lagrangian normalization, bridge completion, or role transfer. If no new formula for J_DHL is supplied, preserve the P3134-P3135 bounded no-go certificate.
