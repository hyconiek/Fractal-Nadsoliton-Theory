# Replaying the combined RB / RBCUT intake

From the repository root, without writing into any incoming source:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_rbcut_weighted_degree_review_20261016/finalize_review.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_rbcut_weighted_degree_review_20261016/verify_intake.py --replay --output /tmp/fin-rbcut-review-replay.json
```

The first command is read-only. The second checks original manifests before
replay, copies all five campaigns to a fresh /tmp directory, runs their
unmodified scientific scripts, records stdout/stderr and compares saved
scientific JSON bytes. It rechecks original incoming hashes afterward.
It writes only the specified output and temporary copies. Python standard
library is sufficient for this intake's finite controls and producers.

Four original package checkers pass. The repair package checker, after
checking its manifest, expects historical /mnt/data input files and fails
there. The verifier preserves that exact failure, checks the copied manifest
independently, resolves the available originals separately and then runs
the scientific script unchanged. It neither suppresses arbitrary failures
nor calls that original checker passing. Missing archived AGENTS and old
wrapper inputs remain unavailable, not synthesized from current files.

The new finite controls use the older versioned AC exact-count helper, a
different q4/N2 rational model with nonuniform weights and orbit stabilizers,
full-cut/orbit-cut enumeration, explicit lifted flows and negative controls.
q12 Fourier multiplicities and native lower-bound arithmetic are checked
separately. Small-N fixtures do not occupy the asymptotic MT domain.

No numerical script proves the analytical statements about every N or
solves the remaining actual cut inequality. Review proofs state their
inherited model, symmetry, moment and boundary hypotheses explicitly.
Full original ZIP/member binding is checked; absent original predecessor
wrappers and past checkout/reading cannot be reconstructed by replay.

`finalize_review.py` without --check is the mechanical authoring helper for
this **new review only**: it validates JSON/AST/links/input hashes and writes
PACKAGE_CHECKS.json and the review manifest. It is not a scientific reference
or a reason to regenerate archived source narratives. Future append-only
AGENTS changes are checked against the bound historical prefix.
