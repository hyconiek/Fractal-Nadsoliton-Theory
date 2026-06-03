# P2480/S1430 strict pointwise interval-Decimal P2459 segment-start cell dyadic-refinement replay audit

Status: `PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_SEGMENT_START_CELL_DYADIC_REFINEMENT_REPLAY_AUDIT_NO_FULL_COMPLEMENT_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM`

## Segment-start weakest-cell dyadic refinement

Parent cell: segment `2`, uncovered index `0`.
Dyadic subcells replayed: `128`.
Parent Decimal separation inherited from P2479: `7.32627829355898709930677580395402985536736133528410588493626676803142121650677344853550E-7`.
Minimum subcell Decimal separation: `9.79535117533551311215134182427543374492277688667778450229304324543911164760697296203721E-7`.
Maximum subcell Decimal separation: `0.000001016933134107930802946894248688882194434615869381378311750866075781669598796178650793367`.
Minimum consecutive positive delta: `2.94456392681960153204638259310811620888142803461710584581507681993440388338218834395E-10`.
Minimum is leftmost subcell: `True`.
All subcells exclude zero: `True`.

## Coverage budget

P2480 subcell replay fraction against inherited P2459 finite universe: `128/99846`.
P2480 parent-cell fraction against inherited P2459 finite universe: `1/99846`.
Inherited P2479 prefix replay count: `1116`.
Inherited P2479+P2478 union count: `1132`.
Full complement replay exported by P2480: `False`.

## Plain-language progress note

This packet opens the weakest P2479 segment-start cell into 128 smaller dyadic subcells. Every subcell still excludes zero with positive Decimal separation, and separations increase from left to right, but the weakest subcell is still the leftmost one. This is finite cell-internal refinement, not a continuum or full-complement proof.

## Hard limits / negative controls

This is a finite dyadic refinement replay of one weakest parent cell only.  It is not a full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, global continuum root-exclusion theorem, selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure.
