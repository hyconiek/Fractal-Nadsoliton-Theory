# P3012/S1962 strict-kernel phase-gradient selector-source obstruction

Status: `P3012_STRICT_KERNEL_PHASE_GRADIENT_SELECTOR_SOURCE_OBSTRUCTION_BOUNDED_NO_GO`

## Selector certificate
- label-based pair pick: `1`
- label-based unit pick: `1`
- unit score order descending: `[1, 5, 7, 11]`
- K(1)-K(11): `0.477799987576`
- Aut equivariance rows/failures: `16/12`
- pair-action leaves-pair count: `4`
- Aut-invariant directed unit count: `0`
- acceptance matrix rows/accepted: `128/1`

## Decision
Bounded no-go: the label score can prefer +1 over -1, but the preference is not Aut(Z12)-equivariant and lacks a strict chart/metric source theorem, so it is not a nonpremise selector source.

## Recommendation
Do not replay kernel label-gradient selector attempts. A next proof-grade move must either supply a genuine strict chart/metric source theorem for a directed selector and rerun equivariance, or introduce a different new typed object outside cube-map, exhausted moment-provenance, selector replay, bridge, role-transfer, and L_total lanes.
