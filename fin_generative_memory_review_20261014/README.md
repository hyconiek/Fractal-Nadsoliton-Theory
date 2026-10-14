# Reviewed FIN memory results

Start with [INTAKE_20261014.md](INTAKE_20261014.md), then
[RG_D06_DECISION.md](RG_D06_DECISION.md).

The original [source campaign](../fin_generative_memory_campaign_20261014/HANDOFF.md)
is retained unchanged. This review admits scoped results and separately adds
[the hidden-refresh gap and moment proof](HIDDEN_REFRESH_GAP_AND_MOMENTS.md).
The next executable plan is
[FIN_GENERATIVE_OPERATIONAL_NEXT_CHAIN_20261014.md](../FIN_GENERATIVE_OPERATIONAL_NEXT_CHAIN_20261014.md).

Local independent verification:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
python3 fin_generative_memory_review_20261014/verify_intake.py \
  --witness-input fin_generative_memory_review_20261014/moment_witness.json \
  --output /tmp/fin-memory-verification.json
```

Run from the repository root. This reads source files only; it writes the
verification and a witness alongside the selected output path. Omit
`--witness-input` to independently regenerate a candidate. Optional `--replay`
points at a six-script scratch replay for environment-comparison diagnostics.
Historical scratch directories are not required for certificate replay.

[verification.json](verification.json) records the actual input hashes,
independent numerical comparisons and outward FIN moment intervals.
Status limitations remain binding even when every computational check passes.
