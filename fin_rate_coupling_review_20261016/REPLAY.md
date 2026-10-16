# Independent rate-coupling intake replay

From repository root; only Python standard library and ordinary Git/coreutils:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_rate_coupling_review_20261016/verify_intake.py --replay --output /tmp/fin-rc-intake-replay.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_rate_coupling_review_20261016/finalize_review.py --check
```

The verifier does not import producer code. It independently checks the
original rational certificate, exact finite fixtures and surviving controls.
With --replay it copies the producer plus its twelve original inputs into
a temporary miniature tree, runs the scientific writer there, and checks
both original output hashes. No archived source is rewritten.

All12 inputs matched before integration. Appending the new AGENTS guardrail
necessarily changes that live file. Replay binds its original snapshot by
`git show 761e0731:AGENTS.md`, with the original SHA in INPUTS.sha256;
the eleven other original files must still match. Do not rewrite the
historical manifest or advertise its live AGENTS entry as still identical.
The older E AGENTS snapshot is verified as an exact hash-bound prefix of
the append-only ledger. Alias input names are resolved explicitly.

verification.json records actual byte replay and its temporary path; the
path is a log, not a required persistent dependency. Exact Fraction bounds
are separate from ordinary math.exp/log diagnostics. No stochastic limit,
finite-data confidence or physical source is certified by replay.

finalize_review.py --check is read-only. Without --check it authors only
this review's PACKAGE_CHECKS.json and MANIFEST.sha256 after syntax, links,
input hashes and conditional next-plan checks.
