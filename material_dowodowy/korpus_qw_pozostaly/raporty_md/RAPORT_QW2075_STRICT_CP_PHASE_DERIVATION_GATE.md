# RAPORT QW-2075: STRICT CP PHASE DERIVATION GATE

- Data UTC: 2026-03-04T03:43:33.009095+00:00
- Verdict: **STRICT_CP_PHASE_DERIVATION_PARTIAL_PMNS_ONLY**
- pass_count: 7/8
- updates promoted: 1

## CKM Phase
- delta_primary [0,pi]: 2.600513 rad
- delta_secondary [0,pi]: 0.541079 rad
- reference: 1.200000 rad
- rel_err primary/secondary: 116.709% / 54.910%

## PMNS Phase
- delta_primary [0,pi]: 0.015591 rad
- sin(delta): 0.015591

## Gate Rule
- CKM update is promoted only if branch-ambiguous phase is compatible with registry tolerance.
- PMNS update is promoted if deterministic and numerically stable (registry has no fixed central value here).

## Artifact
- JSON: `report_qw2075_strict_cp_phase_derivation_gate.json`
