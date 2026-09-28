# 322 — DIRECT-N12-HOLDOUT

Date: 2026-09-27

Status: PASS. Direct finite-N holdout evaluated only after task-321 freeze.

## Microscopic system

Exact count-space state count:

`1,352,078`.

The conservative basin classifier leaves

`345,170` states unlabelled, but their total stationary mass is only

`0.00034544952969208946`.

Thus the labelled localized core carries

`0.999654550470308`

of equilibrium mass.

## Exact C12 reduction

There are 112,720 cyclic orbits. Small periodic orbits are handled exactly by stabilizer-compatible character sectors; they are not discarded.

Exact slow eigenvalues:

- k1: `-0.0008798095680030068`
- k2: `-0.001060751810160859`
- k3: `-0.0005307867688677706`
- k4: `-0.0005813796077097058`
- k5: `-0.0007737319035035645`
- k6: `-0.0007606752905205955`

Hence

`rho_12(exact) = 0.0005813796077097058`

and

`R_exact = [1.51331342953, 1.82454251249, 0.912977961093, 1, 1.33085490658, 1.30839692420]`.

Task-321 clock error:

`2.6077127955 %`.

The smallest fast-sector gap is

`0.475154641814`.

All sector eigenpair residuals are below `1.4e-11`.

## Process holdout

Canonical labelled-core maximum two-time joint error:

`0.02103587617`.

Mean:

`0.02021439380`.

Worst preparation ambiguity bound:

`0.00412970172`.

Worst pair-label ambiguity:

`0.00069077972`.

Omitted-fast burn-in bound:

`8.42e-6`.

Full conservative certified upper:

`0.02288357865`.

This remains far below

`B_313 = 0.05582267457`.

Therefore task 322 is a blind PASS.

Optional task-316 correction improves the certified upper slightly to

`0.02232450920`, but remains unnecessary.

## Core error decomposition

- reduction/history: `0.00032732875`
- cross-N transition: `0.00543305550`
- preparation: `0.01833902508`

Preparation is again the dominant error channel.

## Comparator

Predicted separation before opening:

`0.09039884715`.

After subtracting the complete conservative holdout certificate:

`TV(microscopic FIN, comparator) >= 0.06751526850`.

The directly evaluated core separation is even larger, about `0.10552` TV.

## Barrier trend

The newly observed local clock exponent is

`beta_12 = -log(rho_12/rho_11) = 0.63730565923`.

It lies only

`B4 - beta_12 = 0.02491347790`

below the mapped B4 barrier. This supports the barrier-aware finite-N interpretation, but does not prove that B4 is the global asymptotic communication exponent.
