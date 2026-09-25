# Reproducibility

The predecessor intake correctly criticized constant-printing replay scripts.
This continuation therefore includes real producer scripts for the key new
quantities that are cheap enough to reproduce locally.

## Included producers
- `REPLAYS/replay_strict_memory.py`
  - rebuilds strict q=12 X7 and hidden k=1,2 basis;
  - rebuilds uniform/saddle/localized p states from source roots;
  - recomputes residual hidden covariance eigenvalues and memory strengths.
- `REPLAYS/replay_hierarchy_entropy.py`
  - iterates the exact colored-tree recurrence and outputs h4/h8/h12.
- `REPLAYS/replay_tree_bulk.py`
  - checks an irregular-tree local DtN/Schur identity;
  - checks regular pole/residue scaling formulas.
- `REPLAYS/replay_kurtosis_formula.py`
  - rebuilds strict lambda5 and evaluates the closed k=5 asymptotic formula.
- `REPLAYS/replay_four_phase_partition.py`
  - rebuilds the localized orbit from strict inputs;
  - uses exact unordered count-vector partition recursion;
  - recomputes n=1,2,4,8 minimizer thresholds and alpha=5 collapse.
- `REPLAYS/replay_eight_phase_partial.py`
  - lightweight n=1/n=2 8-phase partition checks only.
- `REPLAYS/replay_inverse_noise.py`
  - verifies the analytic resistance perturbation/topology threshold algebra.

## Not fully reproduced here
The complete 495-composition 12-phase M=16 partition sweep is retained as a
numerical scout record.  A full producer should be exported in the next wave.
The predecessor reaction-valley numerical campaign also remains governed by the
repository intake warning and is not re-certified by this package.
