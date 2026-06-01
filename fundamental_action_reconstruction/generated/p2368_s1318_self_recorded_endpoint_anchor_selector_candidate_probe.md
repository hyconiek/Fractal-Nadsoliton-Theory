# P2368 S1318: self-recorded endpoint anchor selector candidate

Status: `OPEN_PROGRESS_SELF_RECORDED_ENDPOINT_ANCHOR_CANDIDATE_NO_QW2191_DISCHARGE`

## Result

P2368/S1318 promotes the next conditional selector candidate above the phase-origin boundary: a self-recorded ordered d5 endpoint anchor.

## Computed certificate

- Positive five-part ledgers of total 8 enumerated: `35`.
- Unique lexicographic `(ripple, arrow)` winner: `True` with `[[2, 2, 2, 1, 1]]`.
- Source/orientation rows recovered: `24`; all sources `True`; all orientations `True`.
- D12 equivariance checked cases: `576`; mismatches `0`.
- Value multiset source/orientation blind: `True`.
- Reversed ledger moves source to the opposite endpoint and flips orientation: `True`.

## Hard limits

- This is conditional on d5 support, arrow action, endpoint-valued ledger, and endpoint-source convention.
- It is D12-equivariant self-record extraction, not D12-invariant absolute-origin selection.
- No `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.
