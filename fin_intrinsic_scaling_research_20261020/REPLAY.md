# Read-only exact replay

From the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_intrinsic_scaling_research_20261020/verify_research.py --check --output /tmp/fin-intrinsic-scaling-check.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
git diff --check
```

The new standard-library verifier checks finite rational scaling/clique
identities, the exact colour family quotient, count multiplicities,
discrete variance correction, spectral positivity identities, source
integrity, JSON/AST and links. Small generated-source fixtures are not
actual A7 critical trajectories. No q^N enumeration is performed at large N.

--record authors only the fresh IS package; --check writes only the explicit
/tmp replay output. Do not rewrite predecessor evidence or AGENTS. The
old strategy inputs-only check preserves its historical prefix distinction
and does not rerun its completed/stale tasks. No historical scientific
writer/finalizer, Internet, external action, experiment or commit is used.
