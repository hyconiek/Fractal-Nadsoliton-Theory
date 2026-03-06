# RAPORT QW-2140: KERNEL INVERSE FINITE-DOMAIN GATE

- Date UTC: 2026-03-04T19:49:01.376214+00:00
- Verdict: **KERNEL_INVERSE_FINITE_DOMAIN_GATE_PASS_PARTIAL**
- pass_count: `6/7`

## Tested grids
- N=32: cond~`647.074`, near-zero-frac=`0.000e+00`
  exact: err=`2.435e-17`, center/offmax=`1.012e+18`
  regularized: err=`1.489e-03`, center/offmax=`2.592e+03`
- N=40: cond~`835.886`, near-zero-frac=`0.000e+00`
  exact: err=`2.199e-17`, center/offmax=`1.482e+18`
  regularized: err=`2.490e-03`, center/offmax=`1.549e+03`
- N=48: cond~`603.629`, near-zero-frac=`0.000e+00`
  exact: err=`2.166e-17`, center/offmax=`1.873e+18`
  regularized: err=`1.304e-03`, center/offmax=`2.958e+03`

## Scope boundary
- Full continuum action-level Green operator proof: `False`
