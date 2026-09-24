# FIN missing-laws — complete flat handoff, 2026-09-24

This archive consolidates the research chain from
`FIN_MISSING_LAWS_RESEARCH_HANDOFF_20260924.zip` through the current
`NL-11..NL-14` state.

## Design

- No ZIP file is stored inside this ZIP.
- Every nested predecessor/original-package ZIP was recursively expanded.
- Byte-identical payloads were deduplicated globally by SHA-256 and stored only once.
- `SOURCE_MAP.tsv` records every original archive/member and points duplicates to the
  one canonical stored copy.
- Different versions of a report are retained when their bytes differ.
- The older superseded standalone `FIN_MISSING_LAWS_CONTINUATION_HANDOFF_20260924.zip`
  is not included because the canonical continuation chain used here is
  `RESEARCH_HANDOFF -> NL01_NL05 -> NL06_NL10 -> NL11_NL14`.

## Stages

1. `01_RESEARCH_HANDOFF` — four missing-law hypotheses plus the recursively expanded
   post-PHY parent handoff and its unique original-package payloads.
2. `02_NL01_NL05` — first continuation.
3. `03_NL06_NL10` — phase/Fisher/cutoff/radial-coexistence continuation.
4. `04_NL11_NL14` — current radial-stiffness/refinement/static-wall state.

## Deduplication summary

- unique stored source payloads: 279
- duplicate payload occurrences skipped: 1471
- nested ZIP entries expanded (not stored as ZIPs): 38
- kept by stage: {'01_RESEARCH_HANDOFF': 174, '02_NL01_NL05': 38, '03_NL06_NL10': 37, '04_NL11_NL14': 30}
- duplicate occurrences by stage: {'01_RESEARCH_HANDOFF': 209, '02_NL01_NL05': 383, '03_NL06_NL10': 421, '04_NL11_NL14': 458}

`MANIFEST.sha256` verifies every stored file in this consolidated handoff except
itself. `SOURCE_ARCHIVES.sha256` records the hashes of the four source ZIPs.
