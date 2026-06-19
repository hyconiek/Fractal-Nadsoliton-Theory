# P2894/S1844 translation-breaking source alphabet capacity no-go

Status: `P2894_TRANSLATION_BREAKING_SOURCE_ALPHABET_CAPACITY_NO_GO_NO_CLOSURE`

## Finite capacity audit
- low-capacity alphabet sizes tested: `[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]`
- low-capacity Z12-set types tested: `138`
- low-capacity types with equivariant maps to free 12-orbit: `0`
- minimum source orbit size required: `12`

## Named source classes
- `scalar/trivial` [1] -> maps `0`
- `binary sign` [2] -> maps `0`
- `ternary phase` [3] -> maps `0`
- `quarter phase` [4] -> maps `0`
- `half-cycle phase` [6] -> maps `0`
- `full Z12 phase torsor` [12] -> maps `12`

## Boundary
P2894 tests the finite source-alphabet capacity required after P2893.  Every Z12-set source alphabet of total size below 12, including scalar/sign/low-period phase sources, has zero equivariant maps to a free 12-origin carrier orbit.  A free 12-torsor is necessary, but it is not a strict source law without an exported basepoint/polarity and coupling theorem.

## Recommendation
Do not continue with scalar, binary sign, 3/4/6-period phase, low-capacity alphabet, canonical representative, or quotient-section conventions as translation-breaking sourcehood.  A next proof-grade move must either construct an explicit strict free-12-torsor source with a nonimported basepoint/polarity law and coupling theorem to the 9/5 variational density, or pivot to a genuinely different typed object outside the quotient-section/source-alphabet/orbit/Fourier/inventory family; otherwise preserve no-new-live-frontier.
