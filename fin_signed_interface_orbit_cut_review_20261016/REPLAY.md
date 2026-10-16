# Replay of the SC intake and separate new research

From repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_signed_interface_orbit_cut_review_20261016/finalize_review.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_signed_interface_orbit_cut_review_20261016/verify_intake.py --replay --output /tmp/fin-sc-intake-replay.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_balanced_phase_cut_research_20261016/finalize_research.py --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_balanced_phase_cut_research_20261016/verify_research.py --output /tmp/fin-phase-current-replay.json
```

The SC verifier is standard-library-only. It uses hash/check helpers from
the accepted RBCUT review, not imports of incoming producer scientific
code. Original package checker and scientific writer run unchanged only
on a new temporary copy. The source manifest is checked before and after;
the saved source scientific output is byte-identical in this environment.
Historical path references are resolved against actual matching extracted
files and actual ZIPs, with unavailable wrappers/snapshots kept missing.
One source alias duplicates the shell file instead of an MT domain file.

New research has separate exact finite and nominal numerical controls;
its analytical actual-count theorems are not proved by these fixtures.
The first and third commands check mechanical JSON/AST/links/input hashes
and review/research manifests without writing. Authoring-mode packaging
helpers write only their new package metadata, not historical sources.
Future append-only AGENTS changes are checked by the saved prefix hash.
