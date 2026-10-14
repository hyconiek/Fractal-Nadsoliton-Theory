# Source access / Git history limitations

Research date: 2026-10-14. The provided plan is copied unchanged to `SOURCE_PLAN_20261014.md` with local SHA-256 `3c7d51ffa6a3df3715e985bffa873972605fbb75b29352a9b4d3675f6d494f71`. The plan cites baseline **c24773ec**, but that SHA was not verified against GitHub objects here.

## Public remote objects read by web inspection

- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/AGENTS.md — confirms 2026-10-14 guardrails and post-326 reconciliation.
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_post_326_methodology_review/REVIEW_20261014.md
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_post_326_methodology_review/verification.json
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_physical_bridge_campaign_01/PHYS-001/DICTIONARY.md
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_physical_bridge_campaign_01/PHYS-003/PROOF.md
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_physical_bridge_campaign_01/PHYS-004/COUNTERMODELS.json
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_physical_bridge_campaign_01/PHYS-004/CALIBRATION_CONTRACT.json

## Access failures

The GitHub `commits/main` endpoint was restricted in the browsing interface; `git ls-remote` and file downloads in the execution container failed for lack of DNS resolution. Therefore **the latest commit IDs, order, change diff, main HEAD, and exact hash of remote source bytes were not established**. The GitHub homepage listed 1,539 total commits and its root directory was truncated; that does not establish the recent commit sequence. The independently reconstructed A7 implementation faithfully matches the reported 4→3 contrast, but is *not* a byte-verified copy of `safe_seed.py`.

The methodology review and repo guardrails do not show a recoverable FIN-M0 simulator. No Kaggle trajectories or hardware data were imported. This archive is a *new execution package*, not a patched Git checkout. The original repository was not modified.

## Follow-up for a commit-specific audit

From a network-enabled clone, run `git rev-parse HEAD`, `git log -n 20 --date=iso --name-status` and `sha256sum` on the source inputs. Bind the resulting full SHA and source hashes, compare the live source formulas and rerun the producer/verifier. If a formula changed, this handoff cannot be assumed to apply unchanged.
