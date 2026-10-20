# Replay and limitations

Run from an environment with Python 3, NumPy and SciPy:

```bash
python3 verify.py
```

The script rebuilds the strict Z12 W/A, confirms its positive semigroup structure, integrates the Fisher kinetic energy of a frozen heat path, independently evaluates the closed Fisher geodesic lower bound, compares regularized discrete KL sums at six grids, and finite-differences the Fisher entropy-gradient tangent map. Also checks a fixed-escape comparison generator and evaluates its stationary path KL-rate.

All computations are deterministic. It uses float64, ordinary SciPy quadrature and matrix exponentials, without interval arithmetic. Numerical values support only the chosen example. Theorems in `PAPER_INTRINSIC_DYNAMICS_FIN.md` have separate analytic arguments. A verifier `PASS` is not a PDE or formal-math proof, scientific peer review or empirical evidence.
