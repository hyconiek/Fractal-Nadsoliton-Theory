# K6-METASTABLE-ISOLATION-132
## The ±k6 pair has a bounded candidate autonomy window and a concrete competing escape saddle

Date: 2026-09-26

Microscopic process:
exact leave-one-out Gibbs heat bath from report 61.

Status:
- exact k6 invariant-line barrier inherited from report 128;
- competing saddle is a genuine full-X7 index-one stationary point;
- actual heat-bath unstable-manifold adjacency checked numerically;
- local capacity exponents follow by the report-63 reversible local-tube argument;
- global capacity comparison remains open.

## 1. Full-X7 stability window

The ±k6 daughters are full-X7 local minima only in

    g6 < g < g_(3|6),

with

    g6
      = 5.123427551398616,

    g_(3|6)
      ≈ 5.180490619637444.

So the binary interpretation must live inside this finite interval.

## 2. Direct ± switching channel

Along the exact k6 invariant line the separating saddle is the uniform state.

The direct local barrier is

    B_pair(g)
      =
      Phi_uniform-Phi_k6
      =
      -Phi_k6(g).

This is the exact barrier from report 128.

## 3. Competing escape channel to the localized phase

The known main index-one saddle branch also exists in this interval.

Using the actual deterministic heat-bath equation

    p_dot
      =
      softmax(g A7 p)-p,

the one-dimensional unstable manifold of this saddle was integrated in both
directions.

At g=5.13:
- one unstable direction converges to +k6 to distance about 2.5e-8 in p;
- the other converges to the localized minimum to machine precision.

The same adjacency is reproduced at g=5.1452287 and g=5.15.

Thus this is not merely a same-energy stationary point:
it is a direct dynamical basin-boundary saddle between k6 and the main
localized phase in the declared mean heat bath.

Its escape barrier is

    B_out(g)
      =
      Phi_main_saddle(g)-Phi_k6(g).

## 4. Barrier comparison

Since the uniform saddle has Phi=0,

    B_out-B_pair
      =
      Phi_main_saddle.

The main-saddle energy crosses zero at

    boxed:
    g_x
      ≈ 5.150374750802444.

Therefore the known-channel ordering is:

### g < g_x
    B_pair < B_out.

Direct ±k6 switching is cheaper than the known localized escape channel.

### g > g_x
    B_out < B_pair.

Escape from k6 toward the localized phase is cheaper than the direct
within-pair switch through the uniform saddle.

Representative values:

| g | B_pair | B_out | B_out-B_pair |
|---:|---:|---:|---:|
| 5.130000 | 1.23233e-6 | 5.22818e-5 | +5.10495e-5 |
| 5.140000 | 7.81682e-6 | 3.45735e-5 | +2.67567e-5 |
| 5.145228719 | 1.35110e-5 | 2.70219e-5 | +1.35110e-5 |
| 5.150000 | 2.00497e-5 | 2.10513e-5 | +1.00166e-6 |
| 5.150374751 | 2.06174e-5 | 2.06174e-5 | ~0 |
| 5.160000 | 3.78914e-5 | 1.10591e-5 | -2.68324e-5 |

The point

    g_bal
      ≈ 5.145228719489144

balances the direct pair barrier with the barrier-margin to the known escape
channel:

    B_pair
      ≈ Phi_main_saddle
      ≈ 1.3510972867e-5,

and

    B_out
      ≈ 2.7021945734e-5.

## 5. Intrawell relaxation

For the pure k6 order parameter, write

    c=g/g6,
    J=c tanh J.

The deterministic local decay rate inside either k6 well is

    r_well
      =
      1-c sech^2 J.

At g_bal:

    r_well
      ≈ 0.00846711747,

so the dimensionless local relaxation time is about

    tau_well
      ≈ 118.1.

Near birth,

    r_well ~ 2 epsilon,
    epsilon=(g-g6)/g6.

Thus a genuine metastable bit also requires

    N B_pair >> log(tau_well),

not merely B_pair>0.

## 6. Stationary-saddle discovery check at g=5.13

A 1500-start full-7D stationary discovery scan found 18 energy/index classes
at the chosen numerical tolerance, including five index-one energy classes.

Their unstable heat-bath endpoints showed:
- three low-energy saddle classes connect localized minima to localized minima;
- the uniform saddle connects the ±k6 pair;
- the small positive-energy main saddle connects k6 to a localized minimum.

This is strong numerical support for the two relevant local channels above,
but it is not a stationary-exhaustion proof.

## 7. Strengthened verdict

The earlier statement

    "±k6 supplies a binary fiber"

must be narrowed.

A necessary candidate autonomy interval from the known channels is only

    boxed:
    g6 < g < g_x
    =
    5.1503747508...

not the whole full-X7 stability interval.

For

    g_x < g < g_(3|6),

the known localized escape channel is already cheaper than the direct
within-pair switch.

Even for g<g_x, GLOBAL two-state autonomy is not yet proved:
one must exclude any additional lower-capacity exit channel.

That is now a precise capacity problem rather than a branch-search problem.
