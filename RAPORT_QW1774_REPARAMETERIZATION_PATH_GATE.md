# RAPORT QW-1774: REPARAMETERIZATION PATH GATE

- Data UTC: 2026-03-02T18:34:53.437818+00:00
- Global score: 0.804
- Hard gate: **FAIL**
- Readiness: **REPARAMETERIZATION_PATH_STRONG_PARTIAL**

## Core Flags
- mechanism_branch_core: True
- reparameterization_core: True

## Checks
- Leakage-controlled mechanism branch (1768): PASS | score=1.000 | note=LEAKAGE_CONTROLLED_NONENVELOPE_SUPPORTED
- Pre-reparameterization bridge baseline (1772): FAIL | score=0.413 | note=KERNEL_BRIDGE_INTEGRATION_OPEN
- Omega-suppressed projection quality (1773): PASS | score=1.000 | note=OMEGA_SUPPRESSED_LEGACY_PROJECTION_SUPPORTED

## Artifacts
- JSON: `report_qw1774_reparameterization_path_gate.json`
