# Replay of intermediate-clock controls

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_intermediate_clock_bridge_research_20261017/verify_bridge.py --output /tmp/fin-intermediate-clock-replay.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_intermediate_clock_bridge_research_20261017/finalize_research.py --check
```

Exact Fraction fixtures test a rational q3/N2 labelled lift, reversibility,
count lumping, trace rates, clock transforms and barrier/reflection controls.
The separate versioned A7/N2 count calculation uses ordinary numpy/scipy
and is nominal, not outward or an asymptotic proof. The full-N trace and
clock statements are analytical in the theorem file, not proved by sampling.

BROWSER_READ_LOG.jsonl records actual new local browser visits. Full public
DOM text is only in ignored .firecrawl; optional cache absence is reported
without reconstructing original data. Replaying math does not repeat those
visits or attest original PDF hashes. The finalizer authors metadata/hashes
for this new package only and does not modify any predecessor/source.
