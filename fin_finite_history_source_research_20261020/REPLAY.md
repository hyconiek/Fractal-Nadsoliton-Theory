# Read-only replay

From the repository root:

```sh
timeout 55s env PYTHONDONTWRITEBYTECODE=1 python3 fin_finite_history_source_research_20261020/verify_research.py --check --output /tmp/fin-finite-history-check.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
git diff --check
```

The first command is standard-library-only. It writes only the optional
requested /tmp result. Its frozen enumerations are 3^7 candidate words in
one support fixture, 3^6 in the root-bias fixture, and separate 2^4/2^5
controls for periodic stabilizers and finite-record memory. This is not a
growing FIN lattice.
The log controls use exact Fraction series with a proved tail and outward
80-bit dyadic output, not floating near-zero tests.

`--record` authors only this new package's RESULTS, INPUTS and MANIFEST;
it was used when creating the package, not as an old-source replay. Do not
run it to disguise later input changes. `--check` binds exact current bytes
and requires a fresh scoped review if they change. Historical AGENTS-prefix
passes do not replace this package's current full-file context binding.

The literature cache is optional and ignored. Its absence is not a missing
FIN proof or a reason to reconstruct unobserved original PDF hashes. No
old writer/finalizer, network, heavy algebra or external action is required.

The first read-only replay found a JSON round-trip mismatch: integer keys
in the root-count diagnostic became strings on disk. The serializer now
uses explicit string keys before comparison. All 42 arithmetic assertions
had passed; no proof, old source or recorded input was changed by this
local implementation correction. The final read-only replay is required
after recording the corrected new package.
