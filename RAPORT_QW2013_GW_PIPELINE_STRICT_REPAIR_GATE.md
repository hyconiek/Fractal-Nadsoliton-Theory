# RAPORT QW-2013: GW PIPELINE STRICT REPAIR GATE

- Data UTC: 2026-03-03T13:24:51.536263+00:00
- Verdict: **GW_PIPELINE_STRICT_REPAIR_PASS**
- checks pass: 7/7

## Legacy -> Repair Mapping
- [1] pass=True | legacy: projection contradiction + unstable null tail | control: bounded-coupling deep null audit
- [2] pass=True | legacy: lack of scale/local invariance | control: local deterministic perturbation pass-rates
- [3] pass=True | legacy: inter-observatory mismatch | control: explicit control-gap threshold in lockable GW gate
- [4] pass=True | legacy: weak identifiability / unstable significance | control: bootstrap pass-rate + AUC threshold
- [5] pass=True | legacy: possible hard-coded observation drift | control: exclude legacy v61 branch from lockable source pipeline
- [6] pass=True | legacy: phase/lag inconsistency | control: strict lag sanity + no-lag-dependent decision rule
- [7] pass=True | legacy: mixed GW acceptance criteria | control: single explicit deterministic GW flag set

## Required Next Step
- RUN_TRUE_EXTERNAL_BLIND_CONFIRMATORY_WITH_FROZEN_LOCKABLE_PACKAGE

## Artifacts
- JSON: `report_qw2013_gw_pipeline_strict_repair_gate.json`
