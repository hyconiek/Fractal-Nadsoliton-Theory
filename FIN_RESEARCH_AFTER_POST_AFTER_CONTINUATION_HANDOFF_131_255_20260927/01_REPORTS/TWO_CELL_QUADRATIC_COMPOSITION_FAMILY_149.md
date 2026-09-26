# TWO-CELL-QUADRATIC-COMPOSITION-FAMILY-149
## Cell symmetry and diagonal consistency still leave one unsourced dimensionless interaction ratio

Date: 2026-09-26

Status:
exact classification within the minimal symmetric quadratic A7 grammar.

Consider two equal candidate FIN cells with compositions p1,p2.

Require:
1. exchange symmetry 1<->2;
2. the same internal A7 tensor in every quadratic term;
3. diagonal fusion consistency:
       when p1=p2=p,
       recover the original single-cell interaction
       -(g/2) p^T A7 p.

The most general quadratic interaction per total system is

    V_int
      =
      -(1/4)[
        a p1^T A7 p1
        +a p2^T A7 p2
        +2b p1^T A7 p2
      ].

Diagonal consistency gives only

    boxed:
    a+b=g.

So one dimensionless ratio remains free.

Write

    eta=b/g,

then

    a/g=1-eta.

## 1. Positivity / attractive-cross domain

The two-cell block kernel is

    [[a A7, b A7],
     [b A7, a A7]].

Since A7>=0, positivity requires

    a+b>=0,
    a-b>=0.

With g=a+b>0 this gives

    a>=b.

If the cross coupling is also required nonnegative,

    b>=0,

then

    boxed:
    0 <= eta <= 1/2.

Two natural endpoints are:

### independent cells

    eta=0,
    a=g,
    b=0.

Each cell retains the full original interaction independently.

### exchangeable global population split

    eta=1/2,
    a=b=g/2.

Only the average pbar feels A7; the relative mode has zero A7 coupling.

Both satisfy the same diagonal single-cell law.

## 2. Common and relative modes

Linear combinations separate exactly.

The common mode sees

    g_common
      =
      a+b
      =
      g.

The relative mode sees

    boxed:
    g_relative
      =
      a-b
      =
      g(1-2 eta).

Therefore all single-cell diagonal observations are blind to eta, while the
inter-cell relative dynamics depends directly on it.

This is the composition analogue of the earlier fiber-identifiability problem.

## 3. Sensitivity near the current working gain

At

    g=5.145228719489142

the first k6 threshold of a single cell is

    g6=5.123427551398616.

A relative k6 instability would require

    g_relative>g6.

Thus

    eta
      <
      (1-g6/g)/2
      ≈ 0.002118581047.

So at this working gain, a cross coupling larger than only about

    0.2119 %

of g already suppresses the relative k6 instability.

The predicted collective phase structure is therefore extremely sensitive to
the unsourced composition ratio eta.

## 4. No-go

The following conditions do NOT select a unique two-cell law:

- same A7 tensor;
- cell exchange symmetry;
- positive quadratic kernel;
- nonnegative cross coupling;
- exact recovery of the original single-cell law on p1=p2.

A continuous eta-family survives.

Therefore the next physical step cannot be obtained from single-cell
consistency alone.

A genuinely new FIN principle must source the inter-cell coupling ratio or
replace this quadratic grammar with a more fundamental composition rule.
