# Reproduction of *bounded* NT diagnostics

From the extracted campaign directory, with Python 3, NumPy, SciPy and SymPy installed:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 verify_nt.py
python3 package_check.py
```

`verify_nt.py` reconstructs A7 according to the archived MT source constants; checks a **tiny N=2 global Gibbs** q12 count Markov generator (78 count states), detailed balance, stationary measure, symmetric native transfer, nonnegative low spectrum, edge energy, a nontrivial complex two-time reflection Gram matrix, bounded `cos(0.7 m)` correlation identity; and exact symbolic critical operator/ground state identities. This is a sanity fixture **not within** the primary N∈12N local regime and does **not** prove the conditional t,2t stationary convergence or all cylinder theorem.

The only new convergence proof is the mathematical argument in NT-002/STRONG_LIFTED_TRANSFER_PROOF.md, conditional on the named MT/OBS hypotheses. Independent audit of these source hypotheses remains required. Script outputs are ordinary binary64 numerical results; precision is not an interval enclosure.
