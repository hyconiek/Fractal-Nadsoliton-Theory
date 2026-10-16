# Replay of the executed intrinsic preparation

Only Python standard library is required. From repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_activity_covariance_research_20261016/verify_research.py --output /tmp/fin-activity-covariance-replay.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_activity_covariance_research_20261016/finalize_research.py --check
```

The first command writes only its requested output. It does not modify
historical packages, physical records, AGENTS, or its original RESULTS.json.
Exact rational fixtures test the finite arguments. The native A7/g3 check
uses binary64 and is not an interval certificate or a critical-N witness.
Full asymptotic proofs are in THEOREMS.md, with their stated hypotheses;
flag counts do not mechanize them.

The package finalizer checks input hashes, source syntax, links, class scopes
and the next order. --check is read-only. Without it, only generated
INPUTS.sha256, PACKAGE_CHECKS.json and MANIFEST.sha256 are authored.
No external modules, simulations, installations or network are needed.
