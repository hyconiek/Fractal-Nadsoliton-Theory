# Reproduction steps

Extract ZIP, `cd fin_generative_memory_campaign_20261014`.

1. Verify saved files and freeze without altering them: `python independent_verify.py` and `sha256sum -c MANIFEST.sha256`.
2. For a full regeneration, make a SECOND working-directory copy so you do not overwrite frozen preregistration. In that copy run (in this strict order):

```
export OPENBLAS_NUM_THREADS=1
python interval_inputs.py
python train.py
python memory_moments.py
python freeze.py           # freeze must PRECEDE any validation reference
python verify_frozen.py    # independently opens reserved t=.35/.8, delta_2
python independent_verify.py
```

The data-producing scripts are deterministic for fixed numpy/scipy/BLAS. Hashes of float results are environment-sensitive; scientific equivalence means numerical/envelope agreement, NOT automatically byte-for-byte frozen SHA across different versions. The original saved LOCK+MANIFEST bind the delivered artifact.

For exact recheck of Poisson-reference correctness, `model.uniform_tensor` is distinct from sparse `model.uniform_p`. A minimum floating/reference accounting disclosure is at `RG-D02/ERROR_CERTIFICATE.json`. Independent verifier reads raw frozen predictions and reference vectors and recomputes scores without model imports.
