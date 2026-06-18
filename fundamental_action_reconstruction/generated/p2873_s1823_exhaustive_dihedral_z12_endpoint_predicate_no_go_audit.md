# P2873/S1823 exhaustive dihedral Z12 endpoint-predicate no-go audit

Status: `P2873_EXHAUSTIVE_DIHEDRAL_Z12_ENDPOINT_PREDICATE_NO_GO_AUDIT_NO_CLOSURE`

## Exhaustive endpoint-predicate audit
- candidate class: `all Boolean endpoint predicates on Z12 filtered by full D12 symmetry, with a reflection-only comparison boundary`
- predicate count: `4096`
- full dihedral invariant predicates: `[[], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]]`
- reflection-only atoms: `[[0], [1, 11], [2, 10], [3, 9], [4, 8], [5, 7], [6]]`

## Boundary
P2873 exhausts all Boolean endpoint predicates on Z12.  Full D12 invariance leaves only empty/all predicates; reflection-only invariance keeps the unordered pair {1,11} but cannot choose predecessor 11.  Singleton 11 therefore remains an orientation/endpoint-label import, not a strict source law.

## Recommendation
Do not replay full-dihedral endpoint predicates, reflection-only adjacent pairs, cyclic predecessor, or endpoint-label imports as sourcehood.  A next proof-grade move must provide a new strict non-premise chiral/orientation source law selecting singleton 11 and a unit-bearing 9/5 coupling theorem, or pivot to a genuinely different typed object; otherwise preserve no-new-live-frontier.
