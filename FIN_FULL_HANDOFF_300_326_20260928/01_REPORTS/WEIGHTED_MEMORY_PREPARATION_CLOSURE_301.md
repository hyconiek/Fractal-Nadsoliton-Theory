# WEIGHTED-MEMORY-PREPARATION-CLOSURE-301
## Short memory plus vanishing initial-slip weight gives the strengthened closure criterion that report 299 was missing

Date: 2026-09-27

Status:
- exact reversible spectral statement;
- exact N=6 reconstruction from the microscopic generator;
- accepted finite-N N=3..8 slow-overlap and memory-moment inputs;
- no N->infinity theorem.

## 1. Why epsilon_mem alone is insufficient

Report 299 used

    epsilon_mem = rho M1/M0.

This controls the duration of the integrated memory relative to the slow clock.
It does not, by itself, force the initial-slip amplitude to approach one.

The strengthened object must track BOTH:

1. how much spectral weight remains outside the slow mode;
2. how fast that residual weight decays on the slow time scale.

## 2. Exact reversible spectral decomposition

For a normalized reversible resolved correlation,

    C_N(t) = sum_a w_a exp(-r_a t),
    w_a >= 0,
    sum_a w_a = 1.

Suppose the desired slow mode has rate r_N and total weight Z_N. Write

    C_N(t)=Z_N exp(-r_N t)+R_N(t).

If all residual modes carrying nonzero resolved weight satisfy gamma >= gamma_fast,N, then positivity gives the exact bound

    0 <= R_N(t) <= (1-Z_N) exp(-gamma_fast,N t).

This controls the missing preparation/amplitude layer directly.

A sufficient scalar closure programme is therefore not merely

    epsilon_mem -> 0,

but, for example,

    Z_N -> 1,
    r_N/rho_N -> 1,
    relevant residual spectral weight moves to gamma/rho_N -> infinity.

Under these conditions the normalized correlation converges on every fixed positive slow-time window to the single effective exponential.

This is still a statement about the selected resolved sector, not a theorem for the full path process.

## 3. Exact slow-mode weight improves together with epsilon_mem

Accepted values are:

| N | epsilon_mem | exact slow weight Z | exact slip deficit 1-Z |
|---:|---:|---:|---:|
| 3 | 0.0490898 | 0.9108542 | 0.0891458 |
| 4 | 0.0258437 | 0.9239994 | 0.0760006 |
| 5 | 0.0142576 | 0.9583781 | 0.0416219 |
| 6 | 0.00700031 | 0.9720520 | 0.0279480 |
| 7 | 0.00388037 | 0.9851609 | 0.0148391 |
| 8 | 0.00189304 | 0.9913530 | 0.00864705 |

Thus over N=3..8 BOTH diagnostics improve:

    relative memory duration: 4.91% -> 0.189%
    exact initial-slip deficit: 8.91% -> 0.865%.

This is substantially stronger finite-N evidence than report 299 alone.
It still does not prove either quantity tends to zero asymptotically.

## 4. The M1 preparation map predicts the exact slow weight increasingly well

The low-frequency preparation factor

    Z_M1 = 1/(1+M1)

has relative discrepancy from the exact slow-mode weight of approximately:

    N=3: 8.86e-3
    N=4: 4.61e-3
    N=5: 1.92e-3
    N=6: 6.82e-4
    N=7: 2.43e-4
    N=8: 7.18e-5.

So the same memory moment that corrects the clock also predicts the preparation/slip map with rapidly improving finite-N accuracy.

## 5. Exact N=6 process-level scalar certificate

The new direct N=6 diagonalization gives for the k=4 / Z3 resolved correlation:

    Z6 = 0.972052011478719
    r6 = 0.0226112751871153
    gamma_fast = 0.678817607845826.

The first slower full-generator modes outside the target k=4 representation are dark to this observable; their weights are at numerical symmetry-zero level. The first actually resolved fast mode occurs near 0.679.

Therefore

    |C6(t)-Z6 exp(-r6 t)|
      <= (1-Z6) exp(-0.6788176078 t).

Direct checks:

| t | actual absolute error | rigorous spectral-weight bound |
|---:|---:|---:|
| 0.5 | 8.53e-3 | 1.99e-2 |
| 1 | 4.21e-3 | 1.42e-2 |
| 2 | 1.43e-3 | 7.19e-3 |
| 4 | 2.40e-4 | 1.85e-3 |
| 8 | 1.14e-5 | 1.22e-4 |
| 16 | 4.59e-8 | 5.36e-7 |
| 32 | 9.00e-13 | 1.03e-11 |

This directly repairs the over-broad closure criterion identified after report 299.

## 6. Why spectral weight is better than the global hidden gap

The full microscopic chain contains modes slower than the target Z3 mode in other symmetry sectors.
At N=6 the pair near rate 0.02132466 is slower than the target rate 0.02261128, but its overlap with the k=4 resolved observable is about 1e-25 or smaller.

Likewise the N=7 and N=8 accepted computations show extremely slow hidden modes with coupling weights at numerical zero, while the first materially coupled hidden modes remain O(0.5).

So a theorem requiring a uniform gap for EVERY hidden mode would be unnecessarily strong.
The physically relevant target is weighted spectral separation: slow hidden modes are allowed if their coupling weight vanishes.

## 7. Scientific verdict

The strengthened finite-N picture is now:

    short memory duration
      + shrinking unresolved spectral weight
      + accurately predicted preparation/slip factor
      -> controlled late single-mode dynamics.

This is the correct replacement for the claim

    epsilon_mem -> 0 implies Markov closure.

What remains open is to extend the statement from one resolved correlation to the required family of observable finite-dimensional distributions / path statistics, and to prove an N-uniform or asymptotic version rather than infer it from N=3..8.
