# MP7 clean replay

## Predecessor archives

The package contains immutable copies of the three input handoff archives under `inputs/archives/`. Verify `inputs/archives/SHA256SUMS.txt` before extraction/use.

## Existing clean replay

`results/MP7-045_clean_replay.json` records the completed fresh-tree replay. The fresh-tree input and output manifests are copied into `replay/` in the final handoff.

## Minimal executable sequence

Use Python 3 with the already available `mpmath` and `numpy`; no network/package installation is required by this handoff. Set roots to extracted predecessors and the current MP7 tree:

```bash
export R7P_ROOT=/path/to/extracted/fin_rank7_followup
export R7N_ROOT=/path/to/extracted/FIN_R7N_HANDOFF_20260920
export R7O3_ROOT=/path/to/extracted/FIN_R7O3_TARGETP_HANDOFF_20260920
export MP7_WORK_ROOT=$PWD

python scripts/mp7_016_global_exhaustion.py
python scripts/mp7_017_global_first_transition.py
python scripts/mp7_020_spectral_margin_fast.py
python scripts/mp7_021_local_spectral_sensitivity.py
python scripts/mp7_022_023_response_quant.py > results/MP7-022_023_quantitative_response.json
python scripts/mp7_026_controlled_fold.py
python scripts/mp7_034_035_local_weights.py
python scripts/mp7_035_cap_lipschitz_scan.py
python scripts/mp7_035_optimize_explicit_error.py
python scripts/mp7_036_validate_finiteN.py
python scripts/mp7_039_transverse_crossing.py
python scripts/mp7_041_phase_continuation.py
python scripts/mp7_043_check_dag.py
python scripts/mp7_044_negative_controls.py
python scripts/mp7_044_new_controls.py
```

The scripts intentionally expose input roots through environment variables so a replay cannot silently fall back to producer caches.

## Independence level

MP7-045 is a clean-environment replay of the **same mathematical implementations** against copied inputs. It verifies portability, explicit dependencies and deterministic scientific outputs. It does not claim that every theorem was reimplemented independently in a second codebase.
