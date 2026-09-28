# 320 — DIRECT-N11-HOLDOUT

Date: 2026-09-27

Status: PASS. Direct finite-N microscopic holdout evaluated only after task-319 freeze.

## Microscopic size

Exact count-space state count:

`705,432`.

Exact generator transition count:

`47,263,944` nonzero entries.

Stationarity residual from the exact generator build:

`2.3783e-15` in L1.

The conservative high-gap basin classifier leaves stationary mass

`0.0006154901411593956`

unlabelled. It is charged adversarially in the final certificate.

## Exact C12 representation reduction

A global eigensolve on 705k states is unnecessary. Because the microscopic generator commutes exactly with cyclic relabelling, the reversible generator decomposes into C12 Fourier sectors.

The reduction was first validated on opened N=10. It reproduced all six known slow eigenvalues with absolute errors of order `1e-15`.

For N=11 there are 58,786 cyclic orbits, all of size 12. Each sector therefore has dimension 58,786 instead of 705,432.

Exact slow eigenvalues:

- k1: `-0.001637878255016554`
- k2: `-0.001949662537803995`
- k3: `-0.001007642272593455`
- k4: `-0.0010996085944970186`
- k5: `-0.001432188154907444`
- k6: `-0.0014096115978799258`

Hence

`rho_11(exact) = 0.0010996085944970186`

and

`R_exact = [1.48951023411, 1.77305138170, 0.916364493363, 1, 1.30245267459, 1.28192122627]`.

Task-319 clock relative error:

`2.3561078562 %`.

The smallest first-fast-sector decay is

`gamma_fast = 0.462398829543`.

At burn-in t=24 the omitted-fast spectral contribution has certified TV bound at most

`1.1819e-5`;

after the long second interval it is below `1.5e-88`.

## Observed two-time process

Canonical model, core comparison:

`max joint TV = 0.01836456425`.

Mean joint TV:

`0.01637817761`.

Adversarial preparation ambiguity bound:

`0.00733621266` maximum.

Pair-label ambiguity bound:

`0.00123060145` maximum.

After adding all ambiguity and spectral-tail terms:

`certified upper = 0.02240485792`.

This is far below the frozen process envelope

`B_313 = 0.05582267457`.

Therefore task 320 passes without changing any predeclared threshold.

The optional task-316 correction improves the certified upper only modestly to

`0.02164595670`.

It is not necessary for the PASS.

## Error decomposition on the labelled core

- microscopic reduction/history: `0.00050340748`
- cross-N transition transfer: `0.00492934035`
- preparation transfer: `0.01471003512`

The preparation map remains the largest source of error.

## Comparator

Predicted minimum separation before opening:

`0.08493141529`.

After the full N=11 conservative certificate, triangle lower bound:

`0.06252655737`.

Thus the direct microscopic FIN process is still at least 6.25% TV away from this fixed comparator under the declared protocol.

## Boundary

This is a finite-N observed-process result. It does not prove N->infinity, global saddle exhaustion, or a fundamental physical identification of the state variables.
