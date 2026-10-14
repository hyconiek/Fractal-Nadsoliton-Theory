# RG-D04 replay

`OPENBLAS_NUM_THREADS=1 python freeze.py` recomputes frozen prediction files from training inputs and overwrites the original lock, therefore **do not rerun in an audited existing package**. Instead copy the whole research folder to a scratch directory, run the production sequence, compare the newly frozen SHA-256 to the archived immutable originals, then run `verify_frozen.py`.
