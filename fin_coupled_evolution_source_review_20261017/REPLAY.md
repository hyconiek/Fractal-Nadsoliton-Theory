# Read-only checks and copy-only scientific replay

From the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_coupled_evolution_source_review_20261017/verify_intake.py --replay --output /tmp/fin-ue-intake-replay.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_coupled_evolution_source_review_20261017/finalize_review.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_unified_strict_evolution_audit_20261017/finalize_audit.py --check
git diff --check
```

REPLAY/verify.py is a writer: run it only in a temporary copied campaign.
The intake does that itself and preserves original bytes. Source11 labels
are compared semantically with runtime differences separated; byte equality
is reported, not assumed. The independent check uses the accepted safe
archive-expression parser but no incoming scientific producer functions.
Its matrix, full Hessian, trial-integral and Taylor controls have their
stated bounded scopes. Analytical proofs are not replaced by fixtures.

The finalizer authors only this new review's manifest in write mode;
--check is read-only. All required local scientific inputs must hash-match.
The supplied older Library AGENTS is not today's ledger and is not treated
as proof of historical access. No online access or custody is replayed.
