# Read-only checks and preservation

From repository root, with the already installed sympy/numpy/scipy:

```sh
timeout 55s env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_rp004_full_input_self_audit_20261020/verify_audit.py --check --output /tmp/fin-rp004-full-input-self-audit.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
git diff --check
```

The new checker verifies package coverage, 22 inputs, 8 protected manifests,
47 bounded symbolic/rational controls, saved actual replay scope and source
preservation. It writes only the optional /tmp output. It does not rerun a
PDE solver, independently certify uniform asymptotics, or use a web service.

The following read-only replay was actually run before sealing this audit:

```sh
timeout 55s env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_rp004_flux_review_20261019/verify_intake.py --replay --output /tmp/fin-rp004-current-self-audit-replay.json
timeout 55s env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_agents_causal_update_20261020/verify_update.py --check --replay-archives --output /tmp/fin-rp004-self-audit-historical-context.json
```

The first copies scientific writers to /tmp; the existing eight truncated
quadratures remain finite-point, ordinary-floating diagnostics. The second
uses ONLY the byte-verified historical AGENTS prefix in a temporary replay
tree. Seven historical package checks pass there, including RP's checker;
that is ARCHIVED_CONTEXT, not a direct full-live-AGENTS pass or scientific
independence. Originals and today's instructions are unchanged.

Actual outputs are archived here as SOURCE_REPLAY.json and
ARCHIVED_CONTEXT_REPLAY.json. Temporary directory names and floating values
are historical replay metadata, not portable byte-identity promises. A new
environment may differ in final numerical digits. Check scientific scope
before treating any drift as an error; never overwrite old evidence to match.

`--record --source-replay PATH --context-replay PATH` authored only the new
audit metadata/manifest from the observed logs. A syntax typo (one extra
parenthesis) was caught before any scientific execution and corrected.
No proof, source result or input was changed to force a PASS.
