# RAPORT QW-1849: PTA STRICT PATH SELECTION GATE

- Data UTC: 2026-03-02T23:06:09.476492+00:00
- Strict score: 0.250
- Hard gate: **PARTIAL**
- Readiness: **REPARAM_PROTOCOL_REQUIRED_BEFORE_STRICT_CONFIRMATORY**
- Best path: **PATH_B_VERSIONED_CRITERION_REPARAM_WITH_EXTERNAL_CONFIRMATORY**

## Diagnostics
- split-level prob positive: 0.988
- pair-level prob positive: 0.820
- compression gap: 0.167
- p-value H0(prob<=0.90): 0.993306

## Decision Rationale
- Current PTA V1 probability criterion (0.9) is not supported under pair-level analysis; split-level positivity appears inflated by averaging compression.

## Alternative A (keep V1)
- target_n: 76
- additional_needed: 62

## Alternative B (reparam)
- requires new prereg version: True
- requires external confirmatory dataset: True

## Artifacts
- JSON: `report_qw1849_pta_strict_path_selection_gate.json`
