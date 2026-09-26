# Reproducibility

## Baseline
Repository baseline: `FIN_MISSING_LAWS_AFTER_POST_FLAT_FULL_HANDOFF_CONTINUATION_20260925`, commit `62b3e4f8571492a1f2ef33fa29f08d6a53641680`.

## Included artifacts
`ARTIFACTS/` preserves post-checkpoint source code and small/medium text outputs with their original relative paths. Compiled executables, NumPy cache arrays and very large regenerable state dumps are not duplicated; every local artifact, including excluded files, is listed with SHA-256 in `RAW_ARTIFACT_CATALOG.csv`.

## Main runtime families
- Python 3 with NumPy/SciPy/SymPy/mpmath for strict-kernel algebra, interval arithmetic, finite-N generator replays and analysis.
- C++17 producers for exact count-vector recursion, derivative jets, branch-and-bound, tensor enumeration and interval-like directed-rounding certificates.
- FFT-based M=64 hierarchy producers use the exact Gaussian-gauge convolution identity; FFT is an acceleration of the exact recurrence, not a different statistical model.
- Some M=64 depth runs use FFTW where noted by source filename.

## Reproduction strategy
Read each report together with its producer/output. Historical reports are retained even when later superseded; `CORRECTIONS_AND_SUPERSESSIONS.md` and `CURRENT_STATE.md` define the accepted final interpretation.

## Integrity
`MANIFEST.sha256` covers every file in this handoff except itself. `RAW_ARTIFACT_CATALOG.csv` records SHA-256 for every file in the original local post-checkpoint work directory, whether copied or intentionally excluded.
