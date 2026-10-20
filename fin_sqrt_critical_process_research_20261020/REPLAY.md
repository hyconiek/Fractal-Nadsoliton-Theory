# Small read-only replay

From the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_sqrt_critical_process_research_20261020/verify_research.py --check --output /tmp/fin-sqrt-critical-check.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
git diff --check
```

The new verifier uses only the standard library and exact Fraction
arithmetic. It checks coefficients, sparse finite jumps, signs, count
multiplicities, scaling budgets, actual input hashes, protected manifests,
JSON/AST and local links. It is not a stochastic simulation or a mechanical
proof of tightness/UI/uniqueness. No q12 A7 all-type graph is enumerated.

The q12 symmetry fixture is a separately identified rational circulant
rank-five matrix; its nonzero fast curvature is an adverse control to
dropping all fast terms. The q12/N2 Gibbs jump fixture uses A=I and
g=4 log(2). Neither is the asymptotic MT domain or nature's kernel.
The original unaugmented q3/N2 source is an exact c=0 control, not an
asymptotic counterexample to a uniform mixing estimate.

--record authors only this new package; --check is read-only except for
the explicit /tmp output. Do not regenerate old evidence, source hashes,
plans or AGENTS. The strategy inputs-only check preserves its historical
instruction prefix and does not execute its stale NOT_EXECUTED tasks.
No historical scientific writer/finalizer, network acquisition, experiment
or commit is part of these commands.
