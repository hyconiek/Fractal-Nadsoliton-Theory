# FIN generative campaign — execution handoff (2026-10-14)

## Scope / status

RG-C01 to RG-C05 are executed against an independently reconstructed q=12 rank-seven finite model using published formulas, not a cloned/verified checkout. RG-C06 remains reserved for the architect, RG-M00 is blocked on sources.

## Results

- RG-C01: supplied A7 Fourier construction, gain and per-copy leave-one-out CTMC frozen; source checkout SHA unavailable.
- RG-C02: FIN third contrast 0.271270584589 (earlier review ≈0.27127058458855), actual pair TV 0.023358855; mathematically valid tanh(osc/4) nominal certificate 0.051301996. Source interval checked separately where applicable.
- RG-C03: nonzero strong-lumpability defect 0.17191417, K(0)=0.0099520152 and explicit finite-time deviations of observed dynamics from averaged Markov model.
- RG-C04: trace-matched Potts and flat rank-seven numerical controls; see CSV.
- RG-C05: preregistered static computational holdout: FIN |error| 0.0174996 vs pair TV envelope 0.051301996.

## Blockers and next action

The latest commit SHA and full original source/input hashes could not be confirmed: raw files were viewable, `git clone` and GitHub commit history were unavailable. Independently sourced physical A7, FIN-M0 simulator/trajectories, an exact-kernel outward interval certificate for the full F3 oscillation and external experimental data remain absent. Do not promote any of these. RG-C06 requires an architect to decide whether the next branch retains memory, restores FIN-M0, or targets an engineered readout.

## Replay

`python run_campaign.py` recreates numerical outputs; `python independent_verify.py` independently checks claimed invariants. Versions and checksums in `ENVIRONMENT.json`, `INPUTS.sha256`, and `MANIFEST.sha256`.

