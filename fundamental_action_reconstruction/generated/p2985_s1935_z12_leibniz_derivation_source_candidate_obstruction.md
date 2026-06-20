# P2985/S1935 Z12 Leibniz derivation source-candidate obstruction

Status: `P2985_Z12_LEIBNIZ_DERIVATION_SOURCE_CANDIDATE_OBSTRUCTION_BOUNDED_NO_GO`

## Derivation certificate
- modulus: `12`
- additive candidates: `12`
- products tested per candidate: `144`
- accepted derivations: `['D_0']`
- nonzero accepted derivations: `[]`
- acceptance matrix rows/accepted: `256/1`

## Lay summary
P2985 introduces a new finite typed object outside the closed nilradical and CRT lanes: the internal Leibniz derivation algebra of Z/12Z.
Bounded no-go: exhaustive testing of all 12 additive endomorphisms on all 144 products proves that only the zero derivation satisfies Leibniz, so no nonzero strict flow/source is exported.

## Recommendation
Do not replay nilradical, CRT, ratio-package, incidence, selector, or bridge lanes through the zero derivation.  The Z12 Leibniz derivation object is now bounded no-go as a source candidate; the next proof-grade move must introduce a genuinely new strict typed object/theorem/provider outside nilradical/CRT/zero-derivation lanes, or preserve the P2929-P2985 no-strict-export certificate.
