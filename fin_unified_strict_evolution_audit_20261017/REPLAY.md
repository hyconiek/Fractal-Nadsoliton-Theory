# Replay and scope

From the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_unified_strict_evolution_audit_20261017/verify_evolution.py --output /tmp/fin-unified-evolution-replay.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_unified_strict_evolution_audit_20261017/finalize_audit.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_dense_orbit_fibre_transport_review_20261017/finalize_review.py --check
```

SymPy is the symbolic backend, version recorded in verification.json.
The archive reader admits only arithmetic/named functions/Derivative AST
nodes and rejects execution/attributes/floats. It does not import historical
QW writers or execute their Lean wrappers. Exact symbolic checks do not
formalize the analytical PDE proofs. No theorem is inferred from a trajectory.

Browser visits are historical evidence in BROWSER_READ_LOG; math replay
does not recreate browsing or original PDF bytes. Optional public-text/PDF
caches are untracked and may be absent. All required local FIN sources
are hash-bound independently. No external experiment or physical record.

The finalizer writes only this new package's own index/manifest when run
without --check. In check mode it is read-only and verifies required files,
saved checks, source hashes, claim evidence, local links and objective
coverage ledger. AGENTS and incoming science are never changed.
