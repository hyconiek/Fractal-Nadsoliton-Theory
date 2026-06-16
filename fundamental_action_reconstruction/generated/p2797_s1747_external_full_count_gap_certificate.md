# P2797/S1747 external full-count gap certificate

Status: `P2797_EXTERNAL_FULL_COUNT_GAP_CERTIFICATE_NO_CLOSURE`

## External target count and current gap
- external_full_connected_16node_4regular_class_count=8037418
- p2791_distinct_class_lower_bound=8
- p2795_named_subclass_union_class_count=6
- p2795_new_classes_beyond_p2791=0
- uncovered_class_gap_after_p2791=8037410
- p2791_coverage_fraction_exact=4/4018709
- p2791_coverage_fraction_decimal=0.000000995345

## External source
- source=Markus Meringer Regular Graphs page
- url=https://www.mathe2.uni-bayreuth.de/markus/reggraphs.html
- retrieved_utc_date=2026-06-16

## Decision
P2797 sizes the missing full-class obligation against an external GENREG-lineage count: current witnesses cover 8 classes while the reported full connected 16-node 4-regular class count is 8,037,418.  The imported number is a target/gap certificate only, not graph-list provenance or a strict source law.

## Recommendation
Use the external 8,037,418 count only as a target-count/gap obligation.  The next proof-grade move is to import the actual graph-list artifact/toolchain behind that count with graph6/hash provenance and run the full quotient/charpoly/complement/orbit audit, or else export a strict nadsoliton spectral action/source law fixing the admissible class and K/L_total coupling before testing.  Without one of those objects, preserve the P2697-P2797 no-canonical-geometry/no-closure certificate.
