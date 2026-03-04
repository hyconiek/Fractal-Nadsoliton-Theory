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

## Fixed Kernel in This Closure Path

- `omega = 0.185750`
- `phi = 0.162500`
- `beta = 1.000000`
- `eta = 1.800000`

## Important Scope Boundary

Internal closure here means a strict internal gate chain under frozen kernel and locked protocol.

It does **not** mean that all known physical values in nature are already fully derived in final community-accepted form.

Status as of this release:

- **Derived in strict internal gate scope:** mass-chain targets, CKM/PMNS gate-level flavor targets, GW discriminator metrics, and micro-supported renormalization constants (`Z_beta`, `delta_eta`) with tightened dispersion.
- **Still not fully closed globally:** full precision radiative program, full exhaustive Standard Model constant set as final derivation package, and independent external multiteam replication.

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
