# Exact replay and proof scope

From repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_spectator_limit_research_20261020/verify_research.py --check --output /tmp/fin-spectator-limit-check.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
git diff --check
```

The verifier uses stdlib Fraction arithmetic, exact bisection of the actual
two-state stationary-flow equation and outward alternating exponential
bounds. Its finite-c nonuniform pi is fixed across the named test cases;
probabilities are not redesigned at each c to make a scaling root rational.
It also checks a time-warp fixture, initial-law budget, joint-gauge completion,
full-path-TV common-mass countercontrol and an independent labelled/count
rate quotient. The 54 checks do not mechanize the general path theorem.

Authoring --record writes only this new packet's results, bounded duplicate
search, input list and manifest. Later --check is read-only except the
requested /tmp result. Old package manifests and actual live AGENTS are
bound and preserved. Do not run record modes against historical sources or
rewrite their hashes after an instruction update. No network, external
service, heavy computation or metadata writer is needed.

The documentation check excludes code blocks/spans from link parsing, so
mathematical matrix products are not treated as paths. One initial patch
context mismatch was corrected before the verifier file was created; no
mathematical postulate or adverse result was repaired.
