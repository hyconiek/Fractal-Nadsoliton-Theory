# Replay and certificate boundaries

Run from the repository root with Python3, numpy, scipy and mpmath.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
python3 fin_generative_calibration_review_20261014/verify_intake.py \
  --replay /tmp/fin-calibration-intake-nUipS8ry/campaign \
  --output /tmp/fin-f-recheck.json
```

The --replay folder is a source **copy**. To recreate it, copy the F source
package to a new scratch directory and run run_calibration.py,
decision_rule.py and independent_verify.py using python3 and one BLAS thread.
Use that directory in --replay. Source writers must not run in the immutable
original. Reporting/archive stages are not needed for scientific replay;
finalize_package.py has a historical /mnt/data destination.

The intake imports F helper prototypes only for in-memory interface probes.
Scientific process reconstruction uses versioned safe_seed/admitted E
quotient code; threshold intervals use exact-decimal parameters. Synthetic
fixtures have zero physical evidence status. The command writes a verification
and a binary certificate beside --output, never in the source package.

Exact-fraction check of saved event coverage:

```bash
python3 fin_generative_calibration_review_20261014/check_binary_certificate.py \
  fin_generative_calibration_review_20261014/BINARY_EVENT_NO_GO_CERTIFICATE.json
```

All8190 proper events are represented by4096 valid-label masks, with INVALID
included/excluded as described in the proof. Omitted event coverage is a
negative control. Probability boxes inherit the audited E mpmath.iv input
construction; the stdlib checker does not verify transcendental arithmetic.

Check the source and review manifests in their own directories. Numerical
replay/hash identity is not physical calibration, formal IEEE proof, custody,
iid evidence or an independent source for the supplied interaction.
