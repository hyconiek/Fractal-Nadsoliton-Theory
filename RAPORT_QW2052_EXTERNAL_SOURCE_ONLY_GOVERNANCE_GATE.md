# RAPORT QW-2052: EXTERNAL SOURCE-ONLY GOVERNANCE GATE

- Data UTC: 2026-03-04T01:00:06.777315+00:00
- Verdict: **EXTERNAL_SOURCE_ONLY_GOVERNANCE_PASS**
- Readiness: **SOURCE_ONLY_CONFIRMATORY_GOVERNANCE_READY**
- pass_count: 8/8
- size_threshold_bytes: 10485760

## Flags
- data_sources_doc_exists: True
- runbooks_reference_data_sources_doc: True
- release_docs_reference_data_sources_doc: True
- manifest_qw2033_declares_external_sources: True
- manifest_qw2050_declares_data_sources_doc: True
- manifest_qw2033_has_no_embedded_large_archives: True
- manifest_qw2050_has_no_embedded_large_archives: True
- independent_bundle_dirs_have_no_large_payloads: True

## Tracked Large Files (>= 10485760 B)
- count: 16
- `raw_strain_unfiltered/H1_unfiltered_1266965117_long.h5` (20711424 B)
- `raw_strain_unfiltered/L1_unfiltered_1266965117_long.h5` (20711424 B)
- `nadsoliton_scale10.gif` (20614806 B)
- `raw_strain_unfiltered/H1_unfiltered.h5` (16779264 B)
- `raw_strain_unfiltered/H1_unfiltered_1253326744.h5` (16779264 B)
- `raw_strain_unfiltered/H1_unfiltered_1266965117.h5` (16779264 B)
- `raw_strain_unfiltered/H1_unfiltered_1267051517.h5` (16779264 B)
- `raw_strain_unfiltered/H1_unfiltered_1267137917.h5` (16779264 B)
- `raw_strain_unfiltered/L1_unfiltered.h5` (16779264 B)
- `raw_strain_unfiltered/L1_unfiltered_1253326744.h5` (16779264 B)
- `raw_strain_unfiltered/L1_unfiltered_1266965117.h5` (16779264 B)
- `raw_strain_unfiltered/L1_unfiltered_1267051517.h5` (16779264 B)
- `raw_strain_unfiltered/L1_unfiltered_1267137917.h5` (16779264 B)
- `raw_strain_unfiltered/L1_unfiltered_1267569917.h5` (16779264 B)
- `raw_strain_unfiltered/V1_unfiltered_1266965117.h5` (16779264 B)
- `Kopia_notatnika_12_(2)podsumow.ipynb` (11042592 B)

## Independent Bundle Large Payload Check
- large files in independent bundle dirs: 0
- none

## Required Next Step
- Keep freeze bundles source-only and provide this data-source document with each handoff.

## Artifacts
- JSON: `report_qw2052_external_source_only_governance_gate.json`
