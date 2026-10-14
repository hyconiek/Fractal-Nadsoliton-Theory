# Replay and verification scope

Run from the repository root with Python3, numpy, scipy and mpmath.
The source package needs no historical remote checkout for these checks.

## Immutable source integrity

```bash
(cd fin_generative_memory_campaign_20261014 && sha256sum -c MANIFEST.sha256)
```

The supplied independent_verify.py writes a result file. Run it, and all
producer scripts, in a **copy** of the package. Source REPLAY.md specifies
the six-script sequence. Use python3 where the `python` alias is absent.
Set OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1 and
PYTHONDONTWRITEBYTECODE=1. Freeze must precede validation in that copy.
Do not run finalize.py just to replay scientific outputs.

Actual scratch replay: /tmp/fin-memory-intake-WKoEBCYo/campaign. It is not
a permanent scientific dependency; saved verification.json records results
and versions. Different-library replay is expected to alter float hashes.

## Independent repository intake and fixed-witness recheck

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
python3 fin_generative_memory_review_20261014/verify_intake.py \
  --witness-input fin_generative_memory_review_20261014/moment_witness.json \
  --output /tmp/fin-memory-recheck.json
```

The new campaign's model.py/bounds.py are not imported by this checker.
The independent parent uses versioned safe_seed; uniformization uses rate4
and cutoff90, rather than source rate3/cutoff65. Reference comparisons remain
floating numerical results. For the moment certificate, all exact input
functions and weighted residuals are interval evaluated; the exact dyadic
candidate is specified by float.hex strings. The analytic gap is supplied
by the proof, not a floating eigenvalue.

The output path receives the verification and its accompanying witness.
The command does not modify either source package. To regenerate a new
candidate rather than validate the saved one, omit --witness-input.

The hash manifest and verifier checks establish integrity/replay, not original
independent custody, empirical validation or proof-assistant verification.

An additional intake negative control evaluated moment_certificate with an
identically zero candidate function and the exact interval inputs. It was
rejected: the residual cannot support positive tight M0/M1 enclosures. This
checks that a candidate/solver status alone is not sufficient for promotion.
