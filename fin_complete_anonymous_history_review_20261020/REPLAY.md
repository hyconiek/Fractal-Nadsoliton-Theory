# Replay and integrity

From repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_complete_anonymous_history_review_20261020/verify_intake.py --replay --output /tmp/fin-ah-intake-replay-20261020.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_complete_anonymous_history_review_20261020/finalize_review.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_post_rp004_history_review_20261020/finalize_review.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_rp004_flux_review_20261019/finalize_review.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
```

The intake verifier validates the incoming manifest, then runs the inspected
source package checker **before** its scientific writer on a fresh temporary
copy. Each child has a 60-second CPU/wall and 1-GiB address-space bound.
No historical absolute path rewrite is needed by this source. Source RESULTS
is byte-identical here, not universally promised across Python/libm versions.

The new checks use a separate exact rational implementation of kernels,
encodings, departure/edge counting, Euler-tour control and log-series
intervals. They do not import source science functions or simulate the huge
certified length. The certificate checker recomputes rational intervals and
rejects an undersized integer endpoint. Its analytical logarithm remainder
and the martingale proofs are supplied in the proof completion, not hidden
behind a PASS count. No new external literature theorem is applied unread.

`finalize_review.py` writes only this new review's certificate, input binding
and manifest when run without flags. `--check` is read-only. Old sources,
AGENTS, planning snapshot and prior status documents are never rewritten.
Hash/JSON/AST/link checks are reproducibility, not a formal verification of
the entire analytical or physical theory. The original uploaded order remains
absent; actual current inputs are bound separately without fabricated equality.
