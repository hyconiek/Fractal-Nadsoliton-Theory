# Read-only replay

From the repository root:

```sh
prlimit --as=1073741824 --cpu=55 -- timeout 55s env PYTHONDONTWRITEBYTECODE=1 python3 fin_channel_composition_hypothesis_20261021/verify_research.py --check --output /tmp/fin-cc-20261021.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
git diff --check
```

The new verifier uses only the standard library and exact Fraction/integer
arithmetic. It checks coefficient constraints, detailed balance, rational
conditional examples, the N3 event, parity-line/linear-response algebra,
negative controls, named predecessor manifests and bound actual inputs.
It imports no historical scientific writer and writes only the explicitly
requested /tmp result in --check mode. --record authors metadata only in
this new directory; it must not be used to rewrite predecessors.

The actual-A7 N3 theorem is analytical under the bound versioned symmetry
and nonzero-parity premises. The rational q3/q12 examples use another
declared matrix or direct weights; they are not substituted for A7 or MT.
No floating numerical equality, simulation or small graph proves a general
theorem. No browser or PDF cache is a proof dependency here.
