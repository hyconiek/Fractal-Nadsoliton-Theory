# Exact read-only replay

From the repository root:

```sh
prlimit --as=1073741824 --cpu=55 -- timeout 55s env PYTHONDONTWRITEBYTECODE=1 python3 fin_order_free_planar_obstruction_20261021/verify_research.py --check --output /tmp/fin-order-free-20261021.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
git diff --check
```

The new script uses only the standard library, exact integer roots and
Fraction arithmetic. It rebuilds the certificate, tests all twelve
five-port orders, validates actual ordered circular minors, irredundancy,
the perturbation bound, an exact planar-star control and corrupted-sign
controls. It checks bound inputs and named predecessor manifests before
and after, without running historical writers or changing AGENTS hashes.

--record authors metadata only in this fresh directory. --check compares
against those actual bytes, with optional output restricted to /tmp. It
requires neither Internet nor a browser/PDF cache. The external necessity
theorem remains explicitly attributed; it is not proved by script flags.

The actual PDF read is hash-bound in the literature register. If its
temporary cache is absent, report ABSENT_OPTIONAL_CACHE, not a missing
FIN certificate or a reconstructed download. Never copy copyrighted full
papers into the science packet to make a cache-availability flag pass.
