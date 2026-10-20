# Exact replay

From repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_reference_transport_research_20261020/verify_research.py --check --output /tmp/fin-reference-transport-check.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
git diff --check
```

The new implementation uses stdlib exact Fraction and Q(sqrt(5)) arithmetic.
Its 55 bounded checks cover gauge changes, full neighbour counts, weighted
nested embeddings, a 50-state source, the 144->12/78->7 q12 routes, and
positive/negative Markov-closure controls. They are not 55 proofs. No
floating root, huge lattice, matrix simulation or physical record is used.

An early unmatched-parenthesis syntax error in the new verifier was fixed
before its first complete run; it was not a mathematical or model failure.
Scientific counterexamples remain part of the packet. No source model,
postulate, multiplicity, old proof or metadata was repaired.

Authoring --record writes only this fresh directory's result, bounded
duplicate search, input list and manifest. Later --check is read-only
except its requested /tmp output. Predecessor manifests are validated and
the actual live AGENTS hash is bound unchanged. Later instruction changes
must not be disguised by rewriting these historical hashes.
