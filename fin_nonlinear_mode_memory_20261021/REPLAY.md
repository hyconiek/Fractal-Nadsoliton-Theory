# Read-only bounded replay

From the repository root:

```sh
prlimit --as=1073741824 --cpu=55 -- timeout 55s env PYTHONDONTWRITEBYTECODE=1 python3 fin_nonlinear_mode_memory_20261021/verify_research.py --check --output /tmp/fin-nm-20261021.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
git diff --check
```

The new checker uses standard-library integer/Fraction and exact Q(i)
polynomial arithmetic (rational real/imaginary pairs, not noise sampling).
It checks exact profile/eigen identities, moments,
Fourier coefficient reduction, threshold/tail factors, causal signs,
reciprocal force/energy controls and residual time-order bookkeeping.
It does not numerically simulate or mechanically prove the nonlinear PDE,
spectral completeness, asymptotic stationary phase or the stated Gronwall
bootstrap. Those remain explicit analytical proofs/inherited premises.

The finite two-oscillator fixture checks coefficient algebra only; it is
not a discretization of the actual continuum or a new physical law.
Old listed manifest members and actual current inputs are checked before
and after. No historical writer/finalizer is executed. --record authors
only this fresh packet; --check writes only the optional /tmp result.
No external package or browser/PDF cache is a computational dependency.
