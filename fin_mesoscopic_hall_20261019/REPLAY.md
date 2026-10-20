# Bounded read-only replay

From the repository root:

```sh
prlimit --as=1073741824 --cpu=55 -- timeout 55s env PYTHONDONTWRITEBYTECODE=1 python3 fin_mesoscopic_hall_20261019/verify_research.py --check --output /tmp/fin-mesoscopic-hall-check.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
git diff --check
```

The verifier is standard-library-only. It checks exact Gram determinants,
cosets, finite residue periods, projected cell factors, all upward steps
on six exact N120 fixtures, rational exponent balances, a positive-series
bound, exact source cancellation and departure redistribution. It does
not enumerate the growing ten-dimensional lattice or solve a full maxflow.

`--record` was used only to author this new package's RESULTS, INPUTS and
MANIFEST. `--check` writes only the optional requested /tmp output and checks
the stored metadata without rewriting source evidence. Ten predecessor
manifests are checked with their declared excluded helper/metadata files
separately hashed. Current full AGENTS is bound, not falsely equated to an
old full-file input. Later append-only changes require a fresh context review.

No network, browser, PDF cache or old scientific writer is needed for replay.
The browser/PDF record attests actual reading, not a replay dependency or a
permission to reconstruct missing copyrighted files. The analytical shell
and adjacent theorems are not mechanically certified by these finite checks.

A no-op patch attempt encountered a capitalization mismatch before the real
paragraph was inserted; no prior file or scientific inequality changed.
All terminal verification must use the final sealed bytes.
