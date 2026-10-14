# RG-F02 — replay

From campaign root: `OPENBLAS_NUM_THREADS=1 python run_calibration.py`; then `OPENBLAS_NUM_THREADS=1 python decision_rule.py` and `OPENBLAS_NUM_THREADS=1 python independent_verify.py`. Task-specific artifacts are under `RG-F02/`. No physical data is fetched or generated. `write_handoff.py` only writes reports; `finalize_package.py` hashes after these outputs.
