# Reproducible replay (local, NO device access)

From root with Python3, numpy, scipy and mpmath installed:

```
OPENBLAS_NUM_THREADS=1 python run_calibration.py
OPENBLAS_NUM_THREADS=1 python decision_rule.py
OPENBLAS_NUM_THREADS=1 python independent_verify.py
python write_handoff.py
python finalize_package.py
```

`validate_records.py some_file.csv` audits only syntax/required fields and never treats records as independent unless externally established. `independent_verify.py` uses the 144-state labelled generator and its own branch propagation; it does NOT import positive_process nor run_calibration. The derived reporting and manifest stages are deterministic for a fixed local environment. To compare saved calculations across environments, use tolerances stated in independent verifier rather than byte-identical locks. The source context is immutable.

No formal IEEE/mpmath all-operations proof, remote checkout, physical measurements or experimental preregistration is supplied.
