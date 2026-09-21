# MP7-036 — deterministic validation of the finite-N extension

Scientific state: **NUMERICAL_EVIDENCE validating analytic identities**.

The finite-copy model and its main identities are analytic results of MP7-031--033.
This task supplies an independent deterministic implementation check; it does not
upgrade a numerical observation into a theorem.

Replay:

```bash
python scripts/mp7_036_validate_finiteN.py
```

The script uses no sampling and no random seed.  It verifies:

1. direct label enumeration against occupation-count enumeration for `N=1,2,3,4`
   at two nonzero gains; the maximum relative discrepancy is below `9e-16`;
2. the `N=1` constant-diagonal closed form at several gains;
3. invariance of occupation weights under all 24 D12 actions on two nontrivial
   occupation fixtures;
4. the simplex/mediator determinant identity on twelve deterministic interior
   probability vectors, with absolute discrepancies at floating roundoff scale;
5. the exact finite-N source response identity by comparing the analytic
   `N Cov(mu)` with a central finite-difference Jacobian of the independently
   enumerated tilted expectation for `N=2,3,4`; the largest entrywise discrepancy
   is below `8.5e-12`;
6. the auxiliary-field total-covariance identity by independently assembling the
   within- and between-condition covariance terms.

All declared gates pass in `results/MP7-036_validation.json`.

These checks test implementation consistency only.  They do not source a physical
meaning for `N`, `g`, the auxiliary field, temperature or time, and no failure to
sample between phases can arise because no Monte Carlo sampling is used.
