# P2807/S1757 pynauty canonical-label toolchain gate

Status: `P2807_PYNAUTY_COMPACT_CANONICAL_CERTIFICATE_AUDIT_NO_ROW_CSV_NO_SOURCE_LAW_NO_CLOSURE`

## Toolchain probe
- pynauty_available=True
- available_symbols=['Graph', 'certificate', 'canon_label', 'autgrp']
- missing_symbols=[]

## Compact pynauty audit counts
- decoded_graph_count=16828
- pynauty_certificate_hash_class_count=16828
- pynauty_duplicate_certificate_class_count=0
- pynauty_certificate_max_class_size=1
- import_error_type=None

## Decision
P2807 checks the pynauty canonical-label stack and preserves a compact-diff boundary.  If pynauty is unavailable, canonical labeling remains a toolchain blocker; if available, the compact certificate audit is still not a strict spectral source/coupling theorem or K/L_total closure.

## Recommendation
Pynauty is available and P2807 produced unique compact canonical certificate hashes, so the next honest step is P2808: a small canonical digest/automorphism manifest with any row-level CSV kept out of the review diff, followed only then by a separate strict spectral source-law/coupling audit.  Do not promote P2807 to K/L_total, bridge closure, role transfer, selector closure, or ToE closure.
