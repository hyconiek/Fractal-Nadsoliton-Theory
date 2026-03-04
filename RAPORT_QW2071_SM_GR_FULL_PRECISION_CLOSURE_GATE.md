# RAPORT QW-2071: SM+GR FULL PRECISION CLOSURE GATE

- Data UTC: 2026-03-04T10:56:25.714680+00:00
- Verdict: **SM_GR_FULL_PRECISION_CLOSURE_PARTIAL_STRONG_INTERNAL**
- Gate pass count: 3/6

## Gate Flags
- strict_internal_strengthened_pass: True
- full_derivation_package_pass: False
- radiative_program_pass: False
- no_missing_parameters: True
- no_strict_unresolved_parameters: False
- all_radiative_channels_implemented: True

## Coverage Summary
- strict-derived parameters: 28
- model-formula-only parameters: 1
- missing parameters: 0
- strict-unresolved parameters: 7
- implemented radiative channels: 7
- missing radiative channels: 0
- implemented but non-closing radiative channels: 0

## Required Next Steps
- Close all strict-unresolved parameters (non-missing but non-closing statuses still open).
- Upgrade model-formula-only entries to strict derivation status.
- Run truly independent multiteam external confirmation package.

## Artifact
- JSON: `report_qw2071_sm_gr_full_precision_closure_gate.json`
