# RELEASE 4.9: Spectral Closure + First-Principles Internal Strengthening

**Version:** 4.9  
**Date:** 2026-03-04  
**Branch:** `main`

## Executive Summary

Release 4.9 started as spectral micro-bridge closure (QW-2048..QW-2051).  
It has now been extended in the same branch by first-principles internal strengthening (QW-2063..QW-2067):

- deterministic no-scan triad reconstruction (mass + flavor + GW) passes physical thresholds,
- micro-derived renormalization constants gate passes,
- internal strict first-principles closure gate passes,
- dispersion warning is tightened by compatibility-filtered micro aggregation.

Current top-line status:

- internal strict first-principles closure path: **strengthened-pass** (`QW-2067`),
- explicit SM+GR package audit and closure gating: **partial** (`QW-2069`, `QW-2071`),
- external independent multiteam confirmation: **still required**.

## Main Additions in This Path

1. `QW_2048_SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION.py`
- verdict: `SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION_PASS`
- pass_count: `8/8`

2. `QW_2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE.py`
- verdict: `SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE_PASS`
- pass_count: `7/7`

3. `QW_2050_SPECTRAL_MICRO_BRIDGE_FREEZE_BUNDLE.py`
- verdict: `SPECTRAL_MICRO_BRIDGE_FREEZE_BUNDLE_READY`
- pass_count: `4/4`

4. `QW_2051_INDEPENDENT_REHEARSAL_GATE.py`
- verdict: `INDEPENDENT_REHEARSAL_GATE_PASS`
- pass_count: `7/7`

5. `QW_2063_DERIVATIONAL_RECONSTRUCTION_SHARED_FLAVOR_BASIS.py`
- verdict: `DERIVATIONAL_RECONSTRUCTION_TRIAD_PASS_PHYSICAL_PROVISIONAL_FIRST_PRINCIPLES`
- pass_count: `11/12`

6. `QW_2064_MICRO_DERIVED_RENORMALIZATION_CONSTANTS_GATE.py`
- verdict: `MICRO_DERIVED_RENORMALIZATION_CONSTANTS_GATE_PASS_WITH_WIDE_CI_WARNING`
- pass_count: `8/8`

7. `QW_2065_STRICT_FIRST_PRINCIPLES_INTERNAL_CLOSURE_GATE.py`
- verdict: `STRICT_FIRST_PRINCIPLES_INTERNAL_CLOSURE_PASS`
- pass_count: `12/12`

8. `QW_2066_COMPATIBILITY_FILTERED_MICRO_CONSTANTS_TIGHTENING.py`
- verdict: `COMPATIBILITY_FILTERED_MICRO_CONSTANTS_TIGHTENING_PASS`
- pass_count: `6/6`

9. `QW_2067_STRICT_FIRST_PRINCIPLES_INTERNAL_CLOSURE_STRENGTHENED_GATE.py`
- verdict: `STRICT_FIRST_PRINCIPLES_INTERNAL_CLOSURE_STRENGTHENED_PASS`
- pass_count: `3/3`

10. `QW_2068_SM_GR_PARAMETER_REGISTRY.py`
- verdict: registry artifact created
- parameters in scope: `32`

11. `QW_2069_FULL_SM_GR_DERIVATION_PACKAGE.py`
- verdict: `FULL_SM_GR_DERIVATION_PACKAGE_PARTIAL_STRONG_INTERNAL`
- strict-derived: `28/32`, model-formula-only: `1`, anchor-dependent no-fit: `1`, coupled-anchor-dependent: `0`, model-assumption-nonclosing: `0`, SI-definition: `2`, direct missing: `0`, strict-unresolved: `7`

12. `QW_2070_FULL_RADIATIVE_PROGRAM_BASELINE.py`
- verdict: `FULL_RADIATIVE_PROGRAM_PARTIAL_BASELINE`
- implemented radiative channels: `7/7`
- closure-ready radiative channels: `7/7`

13. `QW_2072_EW_YUKAWA_FLAVOR_RADIATIVE_BASELINES.py`
- verdict: `EW_YUKAWA_FLAVOR_RADIATIVE_BASELINES_IMPLEMENTED_NONCLOSING`
- EW + Yukawa + CKM/PMNS baseline blocks added (non-closing by design)

14. `QW_2071_SM_GR_FULL_PRECISION_CLOSURE_GATE.py`
- verdict: `SM_GR_FULL_PRECISION_CLOSURE_PARTIAL_STRONG_INTERNAL`
- pass_count: `3/6`, direct missing parameters: `0`, strict-unresolved parameters: `7`, missing radiative channels: `0`

15. `QW_2073_RADIATIVE_CHANNELS_CLOSURE_UPGRADE.py`
- verdict: `RADIATIVE_CHANNELS_CLOSURE_READY_PASS`
- upgraded channels closure-ready: `5/5` (leading to `7/7` closure-ready in QW-2070 aggregate)

16. `QW_2074_STRICT_NOFIT_MISSING_PARAMETER_DERIVATIONS.py`
- verdict: `STRICT_NOFIT_MISSING_PARAMETER_DERIVATION_ROUND1`
- epistemic no-fit updates added: `2` anchor-dependent physical-relation entries + `2` SI-definition constants

17. `QW_2075_STRICT_CP_PHASE_DERIVATION_GATE.py`
- verdict: `STRICT_CP_PHASE_DERIVATION_PARTIAL_PMNS_ONLY`
- PMNS CP phase promoted to strict internal update; CKM CP phase still non-closing against registry tolerance

18. `QW_2078_GW_EXTERNAL_HOLDOUT_AUTOCOLLECTOR.py`
- verdict: `GW_EXTERNAL_HOLDOUT_AUTOCOLLECTED`
- auto-builds QW-2077 GW observation block from locked metrics/weights

19. `QW_2079_PMNS_CP_EXTERNAL_AUTOCOLLECTOR.py` and `QW_2080_COSMO_WEFF_EXTERNAL_AUTOCOLLECTOR.py`
- verdicts: `PMNS_CP_EXTERNAL_AUTOCOLLECTED` / `COSMO_WEFF_EXTERNAL_AUTOCOLLECTED`
- explicit external-data collectors for PMNS and cosmology blocks required by QW-2077

20. `QW_2081_MISSING14_STRICT_RIGOR_FRONTIER.py`
- verdict: `MISSING14_STRICT_RIGOR_FRONTIER_PARTIAL_ONLY`
- strict map of missing-14:
  - `3` strict candidate target-miss (`delta_cp_ckm`, `h0`, `lambda_cosmological`)
  - `1` anchor-dependent baseline-only (`G_newton`)

21. `QW_2082_MISSING14_STRICT_CLOSURE_ROADMAP.py`
- verdict: `MISSING14_STRICT_CLOSURE_ROADMAP_READY`
- dynamic tiered closure program for current unresolved set (`4` IDs):
  - `T1`: `delta_cp_ckm`
  - `T3`: `h0`, `lambda_cosmological`
  - `T4`: `G_newton`

22. `QW_2083_MISSING14_EPISTEMIC_STATUS_GATE.py`
- verdict: `MISSING14_EPISTEMIC_STATUS_GATE_PASS_WITH_TARGET_MISS`
- deterministic handling for current strict-unresolved subset:
  - `3` strict target-miss,
  - `1` non-closing but explicitly classified,
  - `0` still strictly missing

23. `QW_2084_T1_NONANCHOR_STRICT_GATE.py`
- verdict: `T1_NONANCHOR_STRICT_GATE_PASS`
- pass_count: `6/6`
- strict-nonanchor result for T1 (aggregated from dedicated non-anchor gates):
  - `alpha_s_mz`: PASS (`QW-2087` strict non-anchor gate),
  - `g_f`: PASS (`QW-2085` strict non-anchor gate),
  - `m_z`: PASS (`QW-2086` strict non-anchor gate).

24. `QW_2085_GF_NONANCHOR_LIFETIME_GATE.py`
- verdict: `GF_NONANCHOR_LIFETIME_GATE_PASS`
- pass_count: `5/6`
- strict non-anchor pass reached using kernel-derived muon-lifetime chain with explicit provenance.

25. `QW_2086_MZ_NONANCHOR_EW_POLE_GATE.py`
- verdict: `MZ_NONANCHOR_EW_POLE_GATE_PASS`
- pass_count: `5/6`
- strict non-anchor pass reached using kernel-derived EW-pole input chain with explicit provenance.

26. `QW_2087_ALPHA_S_NONANCHOR_BOUNDARY_GATE.py`
- verdict: `ALPHA_S_NONANCHOR_BOUNDARY_GATE_PASS`
- pass_count: `8/9`
- strict non-anchor pass reached with kernel-derived boundary plus kernel-derived validation points.

27. `QW_2093_KERNEL_DERIVED_NONANCHOR_INPUTS_PLAN_EXECUTOR.py`
- verdict: `KERNEL_DERIVED_NONANCHOR_INPUTS_BUILT_FROZEN_PLAN`
- deterministic frozen-plan builder for QW-2085/2086/2087 input artifacts (no scan, no retune).

28. `QW_2094_STRICT_RIGOR_DEFECT_SWEEP.py`
- verdict: `STRICT_RIGOR_DEFECT_SWEEP_PASS_NO_CRITICAL_DEFECTS`
- consistency checks: `130`, failed: `0` (including `QW-2104` merged preflight, `QW-2105/2106` intake-gap checks, `QW-2107` H(z) design-guidance checks, `QW-2108` G-guidance checks, `QW-2109` evidence-manifest checks, and closure-tooling checks for `QW-2111/2112/2113`).
- sweep now includes pre-gate consistency checks for `QW-2102` (H(z) identifiability), `QW-2103` (G_newton dimensionless provenance), merged T3/T4 preflight consistency (`QW-2104`), intake/gap consistency (`QW-2105/2106`), deterministic H(z) design consistency (`QW-2107`), deterministic G acquisition-spec consistency (`QW-2108`), strict evidence-manifest consistency (`QW-2109`), and operational closure tooling consistency (`QW-2111/2112/2113`).

29. `QW_2095_KERNEL_DERIVED_T2_NONANCHOR_INPUTS_PLAN_EXECUTOR.py`
- verdict: `KERNEL_DERIVED_T2_NONANCHOR_INPUTS_BUILT_FROZEN_PLAN`
- deterministic frozen-plan builder for `QW-2088/2089` input artifacts.

30. `QW_2088_LIGHT_QUARK_MASS_NONANCHOR_GATE.py`
- verdict: `LIGHT_QUARK_MASS_NONANCHOR_GATE_PASS`
- pass_count: `8/9`
- strict non-anchor pass reached for `m_up/m_down/m_strange`.

31. `QW_2089_HIGGS_SELFCOUPLING_STRICT_GATE.py`
- verdict: `HIGGS_SELFCOUPLING_STRICT_GATE_PASS`
- pass_count: `9/10`
- strict non-anchor pass reached for `m_h`.

32. `QW_2096_T2_NONANCHOR_STRICT_GATE.py`
- verdict: `T2_NONANCHOR_STRICT_GATE_PASS`
- pass_count: `7/7`
- aggregate strict non-anchor closure for T2 branch.

33. `QW_2097_CKM_CP_TARGET_REFINEMENT_GATE.py`
- verdict: `CKM_CP_TARGET_REFINEMENT_GATE_TARGET_MISS`
- pass_count: `4/5`
- deterministic permutation/convention audit executed; CKM CP remains outside tolerance without retune.

34. `QW_2090_H0_LAMBDA_DECOUPLING_GATE.py`
- verdict: `H0_LAMBDA_DECOUPLING_GATE_TARGET_MISS`
- pass_count: `7/9`
- strict gate executed with metadata-hardened input checks; strict candidate exists but current external H(z) snapshot still misses H0/Lambda registry tolerances.
- added identifiability diagnostic:
  - two-parameter decoupling lever-arm is weak (`E` span `~0.487`),
  - flatness-projection diagnostic is registry-compatible (`h0` and `lambda` both within tolerance) but does not replace strict decoupled closure claim.

35. `QW_2091_NEUTRINO_ABSOLUTE_SCALE_GATE.py`
- verdict: `NEUTRINO_ABSOLUTE_SCALE_GATE_PASS_STRICT`
- pass_count: `8/8`
- strict gate passes on externalized absolute-scale snapshot with full source metadata.

36. `QW_2092_GNEWTON_SI_BRIDGE_GATE.py`
- verdict: `GNEWTON_SI_BRIDGE_GATE_PENDING_NONCLOSING`
- pass_count: `6/8`
- strict pass is blocked when bridge observable is backsolved from `g_si` (tautology hardening).

37. `QW_2099_HZ_EXTERNAL_DECOUPLING_AUTOCOLLECTOR.py`
- verdict: `HZ_EXTERNAL_DECOUPLING_AUTOCOLLECTED_WEAK_LEVERARM` (current snapshot)
- strict input builder for QW-2090 (`h0_lambda_decoupling_input_qw2090.json`) with source hash and provenance metadata.
- now reports strict identifiability metrics/flags (`n_nodes`, `z_span`, `e_span`, `cond([E,1])`) and supports hard preflight via `--require-strict-ready`.

38. `QW_2100_NEUTRINO_ABSOLUTE_SCALE_EXTERNAL_AUTOCOLLECTOR.py`
- verdict: `NEUTRINO_ABSOLUTE_SCALE_EXTERNAL_AUTOCOLLECTED`
- strict input builder for QW-2091 (`neutrino_absolute_scale_input_qw2091.json`) with source hash and provenance metadata.

39. `QW_2101_GNEWTON_BRIDGE_EXTERNAL_AUTOCOLLECTOR.py`
- verdict: `GNEWTON_BRIDGE_EXTERNAL_AUTOCOLLECTED_BACKSOLVED_NONSTRICT`
- strict input builder for QW-2092 (`gnewton_si_bridge_input_qw2092.json`) now labels bridge origin and marks backsolved inputs as non-strict.
- supports strict provenance mode via `--strict-dimensionless-only --omit-g-si-optional --require-strict-ready` for no-tautology/no-SI-primary preflight.

40. `QW_2098_EW_SECONDARY_NONANCHOR_CLOSURE_GATE.py`
- verdict: `EW_SECONDARY_NONANCHOR_CLOSURE_GATE_TARGET_MISS`
- pass_count: `8/10`
- EW secondary set propagated from strict non-anchor chain:
  - `v_higgs` and `sin2_theta_w_mz` promoted to strict-derived,
  - `m_w` and `alpha_em_inv_mz` remain explicit strict target-miss.

41. `QW_2102_HZ_DECOUPLING_IDENTIFIABILITY_GATE.py`
- verdict: `HZ_DECOUPLING_IDENTIFIABILITY_GATE_WEAK_LEVERARM_PENDING`
- pass_count: `3/7`
- strict-ready input quality is not met for two-parameter H(z) decoupling:
  - failed: `n_nodes_ge_5`, `z_span_ge_0p8`, `e_span_ge_1p0`, `design_condition_lt_8`.

42. `QW_2103_GNEWTON_DIMENSIONLESS_PROVENANCE_GATE.py`
- verdict: `GNEWTON_DIMENSIONLESS_PROVENANCE_GATE_PENDING_NONCLOSING`
- pass_count: `5/8`
- strict-ready provenance for QW-2092 is not met:
  - failed: `bridge_origin_external_dimensionless`, `provenance_anchor_free`, `g_si_not_primary_input`.

43. `QW_2104_T3T4_STRICT_PREFLIGHT_GATE.py`
- verdict: `T3T4_STRICT_PREFLIGHT_GATE_PENDING`
- pass_count: `0/8`
- deterministic meta-gate merges intake (`QW-2106`) + `QW-2099/2102/2090` and `QW-2101/2103/2092` into one strict readiness verdict with contradiction checks.

44. `QW_2105_T3T4_STRICT_INPUT_GAP_REPORT.py`
- verdict: `T3T4_STRICT_INPUT_GAPS_PRESENT`
- deterministic gap report quantifies exact external-input blockers:
  - H(z): intake metadata + nodes/span/condition thresholds,
  - G_newton: intake metadata + dimensionless provenance still backsolved/SI-primary.
- now includes deterministic H(z) design guidance from `QW-2107` (top candidate pairs for new external acquisition nodes),
- and deterministic G acquisition guidance from `QW-2108` (target/range contract for external dimensionless bridge observable).

45. `QW_2106_STRICT_EXTERNAL_INPUT_INTAKE_GATE.py`
- verdict: `STRICT_EXTERNAL_INPUT_INTAKE_GATE_PENDING`
- pass_count: `10/18`
- raw external-input intake now enforces sidecar metadata + strict-ready structure before autocollector stages.

46. `QW_2107_HZ_STRICT_DESIGN_SEARCH.py`
- verdict: `HZ_STRICT_DESIGN_FOUND`
- deterministic no-fit search over acquisition redshift grid proves strict-ready H(z) identifiability is reachable with +2 external nodes.
- current top design suggestions: `[0.10, 0.90]`, `[0.10, 0.92]`, `[0.12, 0.92]` (to be filled with real external measurements, not synthetic values).

47. `QW_2108_GNEWTON_DIMENSIONLESS_ACQUISITION_SPEC.py`
- verdict: `GNEWTON_DIMENSIONLESS_ACQUISITION_SPEC_READY`
- deterministic no-fit contract for external G bridge observable:
  - `mu_ref_gev = 1.0`,
  - `g_dimensionless_target = 6.708830750342e-39`,
  - accepted range: `[6.373389212825e-39, 7.044272287859e-39]`.

48. `QW_2109_STRICT_EXTERNAL_EVIDENCE_MANIFEST_GATE.py`
- verdict: `STRICT_EXTERNAL_EVIDENCE_MANIFEST_GATE_PASS`
- pass_count: `29/29`
- enforces strict evidence-sidecar integrity for external T3/T4 artifacts:
  - required evidence fields (`acquired_utc`, `artifact_sha256`, acquisition protocol/command, analyst),
  - sidecar hash consistency (`artifact_sha256 == sha256(payload)`),
  - sidecar schema/key consistency (`columns_schema` / `json_keys`).

49. `QW_2110_EXTERNAL_EVIDENCE_SIDECAR_TEMPLATE_BUILDER.py`
- deterministic strict sidecar-template builder for QW-2109.
- outputs:
  - `external_hz_nodes_qw2099.metadata.strict.template.json`
  - `external_gnewton_bridge_qw2101.metadata.strict.template.json`
- templates auto-fill payload hashes and schema/key snapshots while preserving explicit manual fields for citation-grade provenance.

50. `QW_2111_T3T4_STRICT_EXTERNAL_ACQUISITION_PACKET.py`
- verdict: `T3T4_STRICT_EXTERNAL_ACQUISITION_PACKET_READY`
- deterministic closure packet for externally blocked T3/T4 channels:
  - converts QW-2105 gaps into operational acquisition requirements,
  - includes top-10 H(z) redshift candidate pairs from QW-2107,
  - includes direct dimensionless G contract from QW-2108,
  - includes strict rerun protocol after data refresh.

51. `QW_2112_HZ_STRICT_NODE_PACK_GATE.py`
- verdict: `HZ_STRICT_NODE_PACK_PENDING`
- pass_count: `2/12`
- validates externally collected H(z) candidate nodes with per-node provenance, performs non-destructive merge with baseline nodes, and checks strict-ready thresholds (`n_nodes/z_span/e_span/cond` + `z>=1.18` coverage).
- emits candidate template:
  - `external_hz_nodes_qw2112_candidates.template.csv`

52. `QW_2113_GNEWTON_DIRECT_DIMENSIONLESS_PACK_GATE.py`
- verdict: `GNEWTON_DIRECT_DIMENSIONLESS_PACK_PENDING`
- pass_count: `1/16`
- validates direct external dimensionless bridge payload against QW-2108 contract (`mu_ref=1 GeV`, accepted range, `origin=external_dimensionless_observable`, `g_si` null) with metadata constraints.
- emits candidate templates:
  - `external_gnewton_bridge_qw2113_direct_dimensionless_candidate.template.json`
  - `external_gnewton_bridge_qw2113_direct_dimensionless_candidate.metadata.template.json`

## Fixed Kernel in This Closure Path

- `omega = 0.185750`
- `phi = 0.162500`
- `beta = 1.000000`
- `eta = 1.800000`

## Important Scope Boundary

Internal closure here means a strict internal gate chain under frozen kernel and locked protocol.

It does **not** mean that all known physical values in nature are already fully derived in final community-accepted form.

Status as of this release:

- **Derived in strict internal gate scope:** mass-chain targets, CKM/PMNS gate-level flavor targets, PMNS CP phase, neutrino absolute mass triad, GW discriminator metrics, and micro-supported renormalization constants (`Z_beta`, `delta_eta`) with tightened dispersion.
- **Package-level audit status:** `28/32` strict-derived, `0` direct missing, and `7` strict-unresolved after T1+T2 promotion plus T3/T4 strict-gate execution and EW-secondary propagation (`QW-2069` + `QW-2083` + `QW-2085..QW-2089` + `QW-2090..QW-2092` + `QW-2093` + `QW-2095` + `QW-2096` + `QW-2097` + `QW-2098`).
- **Missing-14 strict frontier status:** `3` strict candidates still miss target tolerance (`delta_cp_ckm`, `h0`, `lambda_cosmological`) and `1` additional parameter remains non-closing anchor-dependent (`G_newton` in `QW-2081`).
- **Radiative program status:** `7/7` channels implemented and `7/7` closure-ready after QW-2073 upgrade.
- **Remaining package closure gap:** `7` parameters remain strict-unresolved in the full package gate (`QW-2069/2071`), even though direct missing has dropped to `0`.
- **Empirical validation status:** QW-2077 remains mixed/inconclusive until PMNS+cosmology external observations are filled (GW branch already support-ready via QW-2078).
- **Still not fully closed globally:** full precision radiative program, full exhaustive Standard Model + GR constant set as final derivation package, and independent external multiteam replication.

## Next Step

- `RUN_TRULY_INDEPENDENT_MULTITEAM_CONFIRMATORY_PACKAGE`

## External Data Policy (No Large Binary Push)

- Large external payloads are not part of the git freeze bundle.
- Official download sources and commands are documented in:
  - `DATA_SOURCES_EXTERNAL_DOWNLOADS.md`
- Independent teams should acquire raw archives from listed public providers,
  then run frozen scripts with fixed manifests/hashes.

## Textbook Edition

For the full high-school textbook style explanation in English and Polish, see:

- `RELEASE_4_9_TEXTBOOK_EN_PL.md`

## Web-Fetch Strict Update (2026-03-04 UTC)

After external web-fetch completion and strict candidate ingestion:

- `QW-2112`: `HZ_STRICT_NODE_PACK_READY` (`12/12`)
- `QW-2113`: `GNEWTON_DIRECT_DIMENSIONLESS_PACK_READY` (`16/16`)
- `QW-2099`: `HZ_EXTERNAL_DECOUPLING_AUTOCOLLECTED_STRICT_READY`
- `QW-2102`: `HZ_DECOUPLING_IDENTIFIABILITY_GATE_PASS_STRICT_READY` (`7/7`)
- `QW-2101`: `GNEWTON_BRIDGE_EXTERNAL_AUTOCOLLECTED_STRICT_READY`
- `QW-2103`: `GNEWTON_DIMENSIONLESS_PROVENANCE_GATE_PASS_STRICT_READY` (`8/8`)
- `QW-2092`: `GNEWTON_SI_BRIDGE_GATE_PASS_STRICT` (`8/8`)
- `QW-2106`: `STRICT_EXTERNAL_INPUT_INTAKE_GATE_PASS` (`18/18`)
- `QW-2109`: `STRICT_EXTERNAL_EVIDENCE_MANIFEST_GATE_PASS` (`29/29`)
- `QW-2105`: `T3T4_STRICT_INPUT_GAPS_CLOSED_READY_FOR_STRICT_RERUN`
- `QW-2104`: still pending (`7/8`) due only to `QW-2090` strict target miss
- `QW-2090`: still `H0_LAMBDA_DECOUPLING_GATE_TARGET_MISS` (`7/9`)
- `QW-2094`: `STRICT_RIGOR_DEFECT_SWEEP_PASS_NO_CRITICAL_DEFECTS` (`130`, failed `0`)
