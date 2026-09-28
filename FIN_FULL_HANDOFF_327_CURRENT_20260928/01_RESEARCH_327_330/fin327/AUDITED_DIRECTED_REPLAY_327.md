# FIN 327 — AUDITED-DIRECTED-REPLAY-OF-326
## Independent replay of the global B4 separator certificate with directed interval decisions

Date: 2026-09-28

Status: **PASS_DIRECTED_REPLAY / EXTERNAL-REPOSITORY-ACCEPTANCE-PENDING**.

This task was opened because the external review of reports 300–326 correctly identified a proof-status gap in 326: the original global cover used conservative long-double energy margins, but some feasibility / `root_possible` exclusions were not accompanied by an equally explicit directed-rounding error budget. This report does not merely rerun the original code. It replaces the global exclusion logic by interval-valued decisions and independently replays the critical boxes and local saddle/path certificates on a high-precision mathematical kernel.

## 1. What was wrong with treating the original 326 as final

The original 326 source is present in its ZIP, but the report did not surface the full exclusion chain. Inspection of `interval_exhaust_326_final.cpp` shows two qualitatively different rejection mechanisms:

1. an energy lower bound, protected by a deliberately conservative subtraction and a positive acceptance margin;
2. feasibility / fixed-point (`root_possible`) rejection using ordinary long-double evaluations with small hand-added pads.

The first mechanism had an explicit numerical buffer. The second did not have a complete propagated rounding certificate. Therefore the original statement `Gamma=B4 is proved` was stronger than what the surfaced evidence justified.

## 2. Directed interval replay of the six-dimensional half-wall

The exact reduction from 325/326 is retained. The representative separator is

`P0 = P1 >= P2`,

and after eliminating the normal multiplier the constrained stationary problem is a six-dimensional fixed-point problem for the separator dual.

The new checker uses Boost interval arithmetic with directed transcendental rounding policy (`rounded_transc_std<long double>`) throughout the box evaluations. The robust replay additionally encloses every coefficient of the six-dimensional feature matrix by an interval of radius

`2e-15`,

which covers the observed differences between the old float64 construction and an independently evaluated high-precision trigonometric kernel.

A box is rejected only if an interval statement itself proves one of:

- separator infeasibility;
- fixed-point impossibility;
- energy lower bound strictly above the d4 target level;
- membership in a separately certified local d4 neighborhood.

No `+1e-6` acceptance threshold is needed for the directed replay.

### Full-cover result

Robust directed replay:

- processed boxes: **3,402,039**;
- feasibility / fixed-point exclusions: **431,636**;
- energy exclusions: **1,265,741**;
- local d4 boxes: **3,643**;
- unresolved boxes: **0**;
- maximum subdivision depth: **60**;
- minimum positive feasibility exclusion gap: **2.236822820353475e-08**;
- minimum positive energy exclusion margin: **1.1141814885935524e-09**.

Thus the directed implementation exhausts the complete representative six-dimensional half-wall without an unresolved cell.

## 3. Independent high-precision check of the critical global boxes

The two boxes realizing the smallest directed margins were replayed independently with 70-digit `mpmath.iv` interval arithmetic and a high-precision reconstruction of the mathematical FIN kernel rather than the saved float matrix.

For the critical feasibility box:

`gap_exact_kernel >= 4.0806488482909208e-4`.

For the critical energy box:

`energy_lower_bound - V_d4 >= 2.2346782835220542e-7`.

These margins are much larger than the corresponding long-double-directed minima. Hence the global replay's tightest recorded decisions are not artifacts of a float64 representation of the kernel.

## 4. Triple-junction boundary replay

The boundary

`P0=P1=P2=1/3`

must be certified separately because a boundary-constrained minimum need not satisfy the interior six-dimensional fixed-point equation. The exact boundary problem reduces to five dimensions.

Robust directed replay:

- processed boxes: **59,271**;
- feasibility exclusions: **6,740**;
- energy exclusions: **22,896**;
- unresolved boxes: **0**;
- maximum depth: **20**;
- minimum feasibility margin: **1.1752879963434983e-05**;
- minimum energy margin above the d4 level: **3.6046418241732507e-06**.

Therefore the triple junction cannot contain a state at or below the d4 separator energy in this replay.

## 5. High-precision local saddle certificate

The d4 separator stationary point was reconstructed using a high-precision mathematical kernel. Its six reduced coordinates are

`(2.5703048563159231754, 0, 1.2992516656905943475, 0.6817926532921314071, -1.1808995157291637697, 1.9759863787908110659)`

(up to the displayed precision).

A high-precision Krawczyk test gives strict inclusion. On the declared local cube of radius `0.02`, the certified Hessian lower bound is

`lambda_min >= 0.05510885939961514 > 0`.

Thus d4 is a strict constrained local minimum of the separator potential in that neighborhood.

The independently recomputed energies are

`V_localized = -0.80531946214230935605557878998...`,

`V_d4       = -0.14310032501484691985100395454...`,

hence

`B4 = V_d4 - V_localized`

`   = 0.6622191371274624362045748354471982765...`.

This differs slightly in the last digits from earlier float64 summaries and should be preferred when a high-precision value is required.

## 6. Explicit upper-bound path replay

Both explicit segments

`localized minimum -> d4`

and

`d4 -> translated localized minimum`

were replayed on the same high-precision kernel. The interval derivative has the correct strict sign throughout the middle of each segment, and endpoint curvature intervals have the required strict signs.

The worst middle derivative magnitude certificate is approximately

`0.0023927747496585843`,

with opposite signs on the two path halves.

Therefore the explicit path reaches no energy above `V_d4`, establishing the constructive upper bound

`Gamma <= B4`.

The exhaustive separator replay establishes the matching lower bound in this implementation,

`Gamma >= B4`.

## 7. Combined conclusion with report 324

Within the declared finite FIN model and the source/target sectors used by 324, this directed replay supports

`Gamma = B4`.

Together with the theorem of 324,

`lim_{N->infinity} -(1/N) log cap_N(A,B) = Gamma`,

this yields

`lim_{N->infinity} -(1/N) log cap_N(A,B) = B4`

**provided the directed-rounding replay is accepted as the repository's required proof standard**.

## 8. Remaining epistemic caveat

This is materially stronger than the original 326 evidence, because every global rejection decision is interval-valued and the critical decisions were replayed independently at high precision. However:

- the exhaustive 3.4-million-box cover uses the platform's `long double` plus Boost directed transcendental interval policy;
- it is not a complete Arb/MPFR replay of every box;
- no external referee or repository maintainer acceptance is asserted here.

Therefore the correct status is:

**PASS_DIRECTED_REPLAY; external/repository acceptance pending.**

A full Arb/MPFR replay remains a valuable formal-audit lane, not a reason to discard the present result.

## 9. Reproducibility artifacts

Primary sources and outputs in this packet:

- `interval_exhaust_326_robust.cpp` — robust directed six-dimensional cover;
- `MAIN_DIRECTED_ROBUST_STDOUT.txt` / `STDERR.txt`;
- `triple_boundary_exhaust_326_robust.cpp` — directed five-dimensional boundary cover;
- `TRIPLE_DIRECTED_ROBUST_STDOUT.txt` / `STDERR.txt`;
- `critical_mp_interval_326.py` and `CRITICAL_MP_INTERVAL_326.json`;
- `verify_local_path_exactkernel_327.py` and `EXACT_KERNEL_LOCAL_PATH_327.json`;
- original 326 sources retained for comparison.
