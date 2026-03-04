# RUNBOOK QW-2033: Independent Confirmatory Replication

## Scope
- Reproduce combined-branch closure and external preconfirmatory status from frozen artifacts.
- Do not retune sector-specific parameters between mass/flavor/GW.

## Required Checks
1. Verify SHA256 hashes against `manifest_qw2033.json`.
2. Re-run in this order:
   - `python3 QW_2030_FINAL_STAGE_C_GATE_COMBINED_BRANCH.py`
   - `python3 QW_2031_V2_ETA_TRIAD_BLIND_EXTERNAL_VALIDATION.py`
   - `python3 QW_2032_COMBINED_BRANCH_CONFIRMATORY_GATE.py`
3. Confirm final verdict in `report_qw2032_combined_branch_confirmatory_gate.json`.

## Expected Final Status
- verdict: `COMBINED_BRANCH_CONFIRMATORY_GATE_PASS_STRONG`
- readiness: `STAGE_C_PLUS_EXTERNAL_PRECONFIRMATORY_CLOSED`

## External Data Sources (Not Frozen In Git)
- Large raw archives are external by policy (no binary payload push in the bundle).
- Canonical source list: `DATA_SOURCES_EXTERNAL_DOWNLOADS.md`

### Required public sources
- NANOGrav 15yr timing archive:
  - dataset: `NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz`
  - url: `https://zenodo.org/records/16051178/files/NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz?download=1`
  - local example: `external_data/NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz`
  - note: Keep outside git; pass via --nanograv-archive in autocollector scripts.
- GWOSC GWTC event catalog API:
  - dataset: `GWTC event metadata (JSON)`
  - url: `https://www.gw-openscience.org/eventapi/json/GWTC/`
  - local example: `external_data/gwosc_gwtc_eventapi.json`
  - note: Used as external intervention-event source for beta-channel builds.

## Frozen File List
- [OK] `QW_2015_TRUE_EXTERNAL_BETA_CHANNEL_V2_READINESS_GATE.py`
- [OK] `QW_2017_V2_BETA_OBSERVABLE_BLIND_EXTERNAL_INTERVENTION.py`
- [OK] `QW_2021_V2_ETA_OPERATOR_BETA_CONSTRAINT_SCAN.py`
- [OK] `QW_2027_GW_CONTROL_GAP_STRUCTURAL_TERM_SCAN.py`
- [OK] `QW_2029_CKM_BLOCKER_SHARED_FLAVOR_BASIS_SCAN.py`
- [OK] `QW_2030_FINAL_STAGE_C_GATE_COMBINED_BRANCH.py`
- [OK] `QW_2031_V2_ETA_TRIAD_BLIND_EXTERNAL_VALIDATION.py`
- [OK] `QW_2032_COMBINED_BRANCH_CONFIRMATORY_GATE.py`
- [OK] `report_qw2015_true_external_beta_channel_v2_readiness_gate.json`
- [OK] `report_qw2017_v2_beta_observable_blind_external_intervention.json`
- [OK] `report_qw2021_v2_eta_operator_beta_constraint_scan.json`
- [OK] `report_qw2027_gw_control_gap_structural_term_scan.json`
- [OK] `report_qw2029_ckm_blocker_shared_flavor_basis_scan.json`
- [OK] `report_qw2030_final_stage_c_gate_combined_branch.json`
- [OK] `report_qw2031_v2_eta_triad_blind_external_validation.json`
- [OK] `report_qw2032_combined_branch_confirmatory_gate.json`
- [OK] `external_confirmatory_v2/beta_channel_true_external_v2/beta_channel_pairs.csv`
- [OK] `external_confirmatory_v2/beta_channel_true_external_v2/intervention_events.csv`
- [OK] `external_confirmatory_v2/beta_channel_true_external_v2/manifest_beta_channel.json`
- [OK] `external_confirmatory_v2/beta_channel_true_external_v2/protocol_freeze.json`
- [OK] `gw1831_window_features.csv`
