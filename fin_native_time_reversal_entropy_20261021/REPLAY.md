# Read-only replay

From the repository root:

```sh
prlimit --as=1073741824 --cpu=55 -- timeout 55s env PYTHONDONTWRITEBYTECODE=1 python3 fin_native_time_reversal_entropy_20261021/verify_research.py --check --output /tmp/fin-tr-20261021.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
git diff --check
```

The new standard-library verifier uses exact Fraction algebra, formal prime
logarithms and rational logarithm bounds. Tiny jump/finite-step examples
check endpoint identities, reversal preparations, support/nonreversibility
controls, sharp TV counterexamples, count geometry and equivariant records.
They are not critical-domain simulations or a proof of the inherited MT gap.
The tiny reversible reset control has the exactly embeddable skeleton
K=(I+Pi)/2 at time log(2); Pi^2=Pi pays that identity. Its sampled-path
calculations are not relabelled whole continuous paths. Holding-factor
cancellation and endpoint sufficiency are separately derived analytically.

The actual q12 g0 control tests only the declared microscopic normalization.
No numerical native N0 or entropy at an example N is certified. Input and
predecessor hashes are checked before/after; no historical writer runs.
--record authors only this new directory's metadata; --check is read-only
except for an explicitly requested /tmp result. No Internet is required.
Absent PDF cache is ABSENT_OPTIONAL_CACHE, never a fabricated acquisition.
