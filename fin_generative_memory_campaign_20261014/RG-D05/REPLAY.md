# RG-D05 replay

`OPENBLAS_NUM_THREADS=1 python verify_frozen.py` recomputes heldout source references only after verifying the original preregistration hash. `python independent_verify.py` separately verifies the original immutable predictions and scored full laws **without importing the generator**. This command overwrites only post-open validation outputs, not any pre-open file.
