# Read-only replay

From the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_sqrt_local_regime_research_20261020/verify_research.py --check --output /tmp/fin-sqrt-local-check.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
git diff --check
```

The first command verifies this package and its exact finite controls,
bound live inputs, predecessor manifests with coverage, syntax and local
Markdown links. Its --record mode is only for authoring this fresh package;
it must not be used to rewrite predecessor evidence or historical hashes.
The check command writes only the explicitly supplied /tmp output file.
An output-order bug found in the first full replay and its local repair
are retained in [VERIFICATION_NOTES.md](VERIFICATION_NOTES.md).

No external package is needed: all new scientific controls use Python's
standard library and Fraction arithmetic. The q3/N2 fixture is a separate
rational A=I model, not the actual tiny A7/MT domain. The q12/N2 fixture is
at g=0. Neither small control proves an asymptotic stochastic theorem.

The old strategy's inputs-only check validates its 87 historical scientific
inputs and matching AGENTS prefix; it does not execute its old research
tasks or override live instructions. Protected manifests and the full live
AGENTS hash are separately bound by this packet. No historical scientific
writer/finalizer is run here and no Internet acquisition is claimed.
