# RUNBOOK QW-2050: Spectral Micro Bridge Independent Confirmatory

- Bundle generated UTC: 2026-03-04T01:00:04.514738+00:00
- Source verdict: SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE_PASS
- Source readiness: TOE_INTERNAL_BRIDGE_STRICTLY_CLOSED_PENDING_EXTERNAL_MULTITEAM_AUDIT

## Fixed Kernel
- omega/phi/beta/eta: 0.185750 / 0.162500 / 1.000000 / 1.800000

## Integrity Step
1. Verify all SHA256 entries from `manifest_qw2050.json`.
2. Reject bundle if any hash mismatch is found.

## External Data Sources (Not Frozen In Git)
- Large public raw archives are intentionally not included in this bundle.
- Source list and download commands: `DATA_SOURCES_EXTERNAL_DOWNLOADS.md`
- Reproducibility is based on fixed scripts/manifests plus external source provenance.

## Independent Execution
1. Run `python3 QW_2048_SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION.py`.
2. Run `python3 QW_2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE.py`.
3. Confirm `SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE_PASS`.

## Decision Rule
- PASS only if all QW-2049 flags are True and external primary/stress criteria pass.
- Any change in kernel vector invalidates this bundle and requires new freeze.

## Notes
- No sector retune is allowed.
- Use only files listed in manifest.
