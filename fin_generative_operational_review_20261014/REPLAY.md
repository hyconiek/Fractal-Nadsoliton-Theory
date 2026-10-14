# Reproduction and scope

Run from the repository root; dependencies: Python3, numpy, scipy, mpmath.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
python3 fin_generative_operational_review_20261014/verify_intake.py \
  --output /tmp/fin-operational-verification.json

python3 fin_generative_operational_review_20261014/check_gain_certificate.py \
  /tmp/GAIN_BOX_CERTIFICATE.json
```

The first command writes the verification and gain certificate beside the
specified output. It never writes the source package. It independently
reconstructs N2 quotient/parent checks, static intervals and the gain cover.
The second uses exact fractions and rejects a missing coverage row; probability
box validity still depends on the audited interval construction.

To check the immutable saved arithmetic without recomputing inputs:

```bash
python3 fin_generative_operational_review_20261014/check_gain_certificate.py \
  fin_generative_operational_review_20261014/GAIN_BOX_CERTIFICATE.json
```

Source replay ran on `/tmp/fin-operational-intake-agHnEJyY/campaign`, using
run_operational.py, interval_static.py and independent_verify.py in that
order. All three write outputs; run them on a **copy**, not the audited
original. make_reports.py writes authorial status/preregistration additions
and is not required for scientific replay. The scratch path is historical,
not a permanent dependency. Optional --replay attaches its verification.

Check the source manifest in its own directory. The review's numerical
generator/exponential checks do not certify all floating operations. Static
and gain-cover certificates assume mpmath.iv correctness. All source/model
inputs remain supplied, and no physical measurements are implied by PASS.
