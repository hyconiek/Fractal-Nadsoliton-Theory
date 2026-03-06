# RAPORT QW-1827: GW REDESIGN GATE

- Data UTC: 2026-03-02T22:06:52.421142+00:00
- Global score: 0.176
- Hard gate: **FAIL**
- Readiness: **GW_BRANCH_REDESIGN_REQUIRED**
- Recommendation: **DO_NOT_RUN_JOINT_GW_PTA_CAMPAIGN_YET**

## Checks
- Null rejection (phase + shift): FAIL | score=0.000 | priority=critical
- Lag coherence at +10ms: FAIL | score=0.082 | priority=critical
- Scale/window stability: FAIL | score=0.374 | priority=major
- Detector consistency (H1-L1 vs H1-V1): FAIL | score=0.000 | priority=critical
- Shared-control identifiability (1826): FAIL | score=0.424 | priority=major

## Priority Order
- Null rejection (phase + shift)
- Detector consistency (H1-L1 vs H1-V1)
- Lag coherence at +10ms
- Scale/window stability
- Shared-control identifiability (1826)

## Action Items
- Rebuild GW statistic around event-windowed, detector-coherent features (chirp-conditioned) instead of global cross-Hurst only.
- Introduce explicit shared-vs-unshared contrast objective and require positive near-target prevalence advantage across SNR sweep.
- Enforce multi-detector consistency gate (H1-L1, H1-V1, L1-V1) before any projection claim.
- Require lag-coherence threshold around physical delays as hard precondition for structure claims.

## Artifacts
- JSON: `report_qw1827_gw_redesign_gate.json`
