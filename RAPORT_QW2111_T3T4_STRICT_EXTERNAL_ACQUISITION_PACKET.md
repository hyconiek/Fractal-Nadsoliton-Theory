# RAPORT QW-2111: T3/T4 STRICT EXTERNAL ACQUISITION PACKET

- Date UTC: 2026-03-04T13:28:39.381496+00:00
- Verdict: **T3T4_STRICT_EXTERNAL_ACQUISITION_PACKET_READY**

## Status Snapshot
- QW-2105: `T3T4_STRICT_INPUT_GAPS_CLOSED_READY_FOR_STRICT_RERUN`
- QW-2106: `STRICT_EXTERNAL_INPUT_INTAKE_GATE_PASS` (18/18)
- QW-2109: `STRICT_EXTERNAL_EVIDENCE_MANIFEST_GATE_PASS` (29/29)

## H(z) Acquisition
- Current n_nodes: `5` (threshold `5`)
- Current z_span: `1.15` (threshold `0.8`)
- Current E_span: `4.276715605330501` (threshold `1.0`)
- Current cond([E,1]): `5.7956053033256225` (threshold `< 8.0`)
- Suggested added z pairs (top 10): `[[0.1, 0.9], [0.1, 0.92], [0.12, 0.92], [0.1, 0.94], [0.12, 0.94], [0.1, 0.96], [0.12, 0.96], [0.14, 0.96], [0.1, 0.98], [0.12, 0.98]]`

## G Bridge Acquisition
- Target mu_ref_gev: `1.0`
- Target g_dimensionless: `6.708830750342087e-39`
- Accepted range: `{'min': 6.373389212824983e-39, 'max': 7.044272287859191e-39}`
- Hard origin requirement: `external_dimensionless_observable`

## Artifact
- JSON: `report_qw2111_t3t4_strict_external_acquisition_packet.json`
