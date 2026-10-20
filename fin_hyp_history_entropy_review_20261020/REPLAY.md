# Replay

From repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_hyp_history_entropy_review_20261020/verify_review.py --check --replay-source --output /tmp/fin-hyp-history-review-check.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
git diff --check
```

Before source replay the new helper validates the original manifest and
every saved source input hash. The source runs only with --check, in a
temporary cache-free copy with explicitly needed read-only dependencies.
Its 48 controls reproduce its saved JSON exactly. Its --record writer and
duplicate search are not run. The subprocess has 60-second/1-GiB limits.

The 68 new checks use separate code and stdlib exact fractions, prime-log
expressions, outward rational log/square-root bounds and a frozen modest
root bracket. Full short-word probabilities test the path identity; source
formulas are not merely repeated as their own independent reference.
These are finite diagnostics, not 68 probability-theorem proofs.

Authoring --record writes only this new review's verification, input list
and manifest. Later --check writes only the requested /tmp output. If
instructions change later, keep these actual historical hashes; do not
rewrite old evidence to make a full-file check green. The package is not
a new network acquisition or an independent-custody attestation.
