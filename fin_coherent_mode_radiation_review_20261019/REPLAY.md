# Replay and instruction contexts

From the actual repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_coherent_mode_radiation_review_20261019/verify_intake.py --replay --output /tmp/fin-rp-intake-replay.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_coherent_mode_radiation_review_20261019/finalize_review.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_agents_update_20261019/verify_update.py --check
git diff --check
```

The inspected producer script executes only in a temporary source copy.
Its16 flags use mpmath point quadrature/differentiation; no exact/outward
or PDE-proof status is inferred. The intake uses the accepted safe archive
parser, independent symbolic differential/Fourier identities, rational
margins and finite-interval quadrature. No scientific producer routine is
imported. Analytical completeness and flux limits are stated separately.

AGENTS is bound by the actual byte length and hash seen during this intake.
A later append can be recognized as a separate instruction tail; it is not
silently equated with this revision. Earlier source hashes are checked on
their verified historical prefix. The older instruction-update helper can
report a later tail outside its own frozen revision; this review supplies
the additional current context, not a rewrite of its contract.

Archived full-file AGENTS validators require the explicit historical-context
replay route, not direct execution against a changed ledger. The new
finalizer writes only this new review's metadata/manifest in write mode;
--check is read-only. No browser, external experiment or Git write is run.
