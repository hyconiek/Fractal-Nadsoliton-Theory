# CERTIFIED-FOLD-DYNAMIC-SCALING-59
## The interval-certified FIN fold implies square-root critical slowing in the declared maximum-entropy heat-bath dynamics

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Inputs:
- `fin_rank7_followup/certificates/R7P-031_simple_fold.json`;
- the exact heat-bath dual-coordinate dynamics derived in reports 54-58.

Status:
- local asymptotic theorem using already interval-certified fold coefficients;
- numerical coefficient interval propagated from the accepted certificate;
- not a physical-time prediction and not a global first-transition theorem.

## 1. Heat-bath dynamics in dual retained coordinates

For fixed gain g let the reflection-even dual field be s and let

    m(s)=C4^T softmax(C4 s).

The dual stationary gradient is

    F(s,g)=grad Phi_g(s)=s/g-m(s).

The retained mean heat-bath equation is

    dot mu = m(g mu)-mu.

With

    s=g mu

and constant g,

    boxed:
    dot s
      = g m(s)-s
      = -g F(s,g).

So the accepted stationary fold is automatically a dynamical saddle-node of
the declared mean heat-bath flow.

## 2. Certified fold data

R7P-031 certifies one simple fold near

    g_f = 3.51564471683959...

with normalized kernel vector v and the interval signs

    a := v . F_g
       in
       [-0.212571640084014,
        -0.212571626020459],

    b := v . F_ss[v,v]
       = D^3 Phi[v,v,v]
       in
       [0.118982527450499,
        0.118982667103224].

Therefore

    a<0,
    b>0

throughout the certified fold box.

## 3. Center-manifold normal form

Let

    delta = g-g_f

and y be the coordinate along the fold kernel v.

The simple-fold expansion gives

    v.F
      = a delta + (b/2)y^2
        + higher-order terms.

Hence the heat-bath flow has leading center equation

    boxed:
    dot y
      = -g_f [
          a delta +(b/2)y^2
        ]
        + higher-order terms.

## 4. Which side contains the two branches?

Stationarity requires

    a delta +(b/2)y^2=0,

so

    y^2
      ~ -2a delta/b.

Because

    a<0,
    b>0,

real nearby branches occur for

    boxed:
    delta>0,
    i.e. g>g_f.

Thus the local saddle-node creates/annihilates its pair on the high-g side of
the certified fold.

This is a local statement; it does not identify the fold as the first global
event.

## 5. Branch displacement

The leading branch amplitudes are

    y_+/- =
      +/- sqrt(-2a/b) sqrt(g-g_f)
      +O(g-g_f).

Propagating the certified coefficient intervals gives

    boxed:
    1.89027850
      < sqrt(-2a/b)
      < 1.89027968.

So the local branch separation has a tightly determined square-root amplitude.

## 6. Stable versus unstable branch

Differentiate the center flow with respect to y:

    partial_y dot y
      ~ -g_f b y.

Therefore:

- y>0 has negative linear eigenvalue and is the locally stable branch;
- y<0 has positive linear eigenvalue and is the locally unstable branch.

This orientation uses the accepted sign convention for v in R7P-031.

Flipping v flips the labels +/- but not the stable/unstable physical content.

## 7. Critical slowing coefficient

On the stable branch,

    y_+
      ~ sqrt(-2a/b) sqrt(delta).

The magnitude of the slow relaxation eigenvalue is therefore

    lambda_slow
      ~ g_f b y_+

      = g_f sqrt(-2ab) sqrt(delta).

Using the full certified intervals,

    boxed:
    0.790704515
      < g_f sqrt(-2ab)
      < 0.790705010.

Hence

    boxed:
    lambda_slow
      ~ 0.7907048 sqrt(g-g_f)

in the declared unit-rate heat-bath clock.

Equivalently the slow relaxation time diverges as

    tau_slow
      ~ [1/0.7907048] (g-g_f)^(-1/2).

## 8. Gaussian fluctuation implication

From the fluctuation-dissipation theorem of report 57, the local Gaussian
variance is controlled by the inverse static curvature.

Since the soft curvature vanishes as

    O(sqrt(g-g_f))

on the stable saddle-node branch, the corresponding Gaussian variance diverges
as

    O((g-g_f)^(-1/2))

inside the local Gaussian approximation.

Thus the same certified fold predicts both:

    relaxation time ~ (g-g_f)^(-1/2),

and

    soft-mode variance ~ (g-g_f)^(-1/2),

for the declared heat-bath dynamics.

## 9. Exact scope boundary

This result does NOT establish:

- a physical second;
- a physical critical exponent measured in nature;
- that g is temperature;
- that R7P-031 is the first/global physical transition;
- validity of the Gaussian approximation arbitrarily close to the finite-N
  fold.

It is a mathematically forced local dynamic consequence of combining:

1. the already certified static simple fold;
2. the maximum-entropy heat-bath evolution.

## 10. Significance for the FIN programme

This is the first point in the present continuation where an already
interval-certified static singularity yields a quantitative dynamic scaling
law without fitting an additional kinetic mode.

The remaining free factor is the overall refresh/event clock rate.
