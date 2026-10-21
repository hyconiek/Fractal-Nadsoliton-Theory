# Verification notes

Date: 2026-10-21.

The first authoring run completed its scientific algebra but stopped in
the link checker: REPORT_PL.md already linked RESULTS.json, which --record
had not yet created. No result manifest was authored by that failed run.
The metadata repair moves link checking after metadata creation.

A subsequent read-only replay rejected RESULTS because multithreaded rg
returned identical matching lines in a different order. Exact input hashes
and every scientific control still agreed. The saved search now sorts
those lines before comparison. These two packaging/replay repairs change
no rate, analytical claim, fixture, old input or scientific bound; no second
scientific route or analytical repair was introduced. The final --check
run must validate the complete authored packet read-only.

The checker also uses fixed algebraic s fixtures, not a fitted parameter
scan. Its count of assertions includes repeated substitutions of defining
identities; it is not that many independent proofs. The actual A7 N3
argument is analytical, whereas the rational event number belongs only to
the separately declared integer-potential control.
