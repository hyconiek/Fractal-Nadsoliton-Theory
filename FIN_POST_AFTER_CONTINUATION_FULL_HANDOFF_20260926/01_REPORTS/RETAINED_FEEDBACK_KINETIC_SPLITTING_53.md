# RETAINED-FEEDBACK-KINETIC-SPLITTING-53
## Isotropic microscopic refresh plus retained-only feedback generates strict spectral kinetics without direct hidden coupling

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
- exact linearization inside the declared empirical heat-bath model;
- exact uniform Fourier-sector rate formula;
- exact stationary visible/hidden decoupling formula;
- conditional on the supplied heat-bath refresh law.

## 1. Declared target dependence

The heat-bath target depends on the retained sufficient statistic only:

    mu = X^T p.

For a stationary state p*, perturb retained coordinates by x.

The accepted target expansion gives

    delta q
      = g S X x,

where

    S=diag(p*)-p* p*^T.

Therefore

    delta mu_dot
      = (gF-I)x,

with

    F=X^T S X.

## 2. Hidden raw drift

Let Y span the four hidden k=1,2 directions and define raw hidden fluctuation y.

Its linear drift is

    dot y = -y + g H x,

where

    H=Y^T S X.

The target has no independent hidden argument; hidden forcing appears only
through correlation with the retained perturbation.

## 3. Conditional residual decoupling

Define

    K=H F^{-1},
    z=y-Kx.

Then

    dot z
      = (-y+gHx)
        -K(gF-I)x
      = -y+Kx
      = -z.

Thus at every admitted stationary state with F invertible,

    boxed:
    dot x = (gF-I)x,
    dot z = -z.

This is the deterministic part of the accepted Gaussian stationary theorem.

The hidden residual has one universal dimensionless decay rate even though the
visible retained sector has a nontrivial spectrum.

## 4. Uniform state

At p*=1/12,

    S=P_C/12,
    H=0,
    F=X^T X/12.

Since X^T X has retained strict eigenvalues lambda_3,...,lambda_6,

    visible decay rates:
      gamma_k = 1-g lambda_k/12,
      k=3,4,5,6;

    hidden decay rates:
      gamma_1=gamma_2=1.

In label space the centered linear drift is exactly

    boxed:
    B_lin = -P_C + (g/12) A7,

where

    A7=X X^T.

So the rank-seven projected strict operator, not the full strict A, carries the
feedback.

## 5. Numerical rates at coexistence

At

    g_eq=3.7183448981203875,

the retained decay rates are

    k3: 0.392234400136
    k4: 0.318437032585
    k5: 0.287749091286
    k6: 0.274246613070,

while the two hidden Fourier planes retain rate 1.

These reproduce the rates already used in `EDGEWORTH-STATIONARY-15`.

## 6. Direct A-coupling cannot mimic this pattern

A direct generator

    -P_C-rho A

would assign

    1+rho lambda_k

to ALL six nonzero Fourier sectors.

To preserve equal hidden rates at exactly 1 requires

    rho=0

because lambda1 and lambda2 are nonzero/distinct.

Then all retained rates also equal 1.

Therefore no nonzero direct-A kinetic coupling reproduces the heat-bath
retained-only splitting.

The two mechanisms are structurally different:

    direct strict mobility:
        A acts on hidden + retained sectors;

    heat-bath feedback:
        A7 acts only through retained sufficient statistics.

## 7. Emergent-role interpretation

The declared heat-bath therefore contains a clean two-layer mechanism:

    isotropic microscopic forgetting
        -I

    plus retained self-consistency
        +g F
        / +g A7/12 at uniformity.

Hidden directions are forgotten directly.
Retained directions are forgotten and simultaneously regenerated through the
macrostate-dependent target.

This supplies a precise sense in which information selected as "retained" can
acquire slower effective dynamics without assigning a separate microscopic
clock to each mode.

## 8. Source boundary

The mechanism is mathematically exact inside the declared heat-bath model.

Still supplied rather than derived:
- why microscopic updates are complete refresh;
- why the target depends exactly on X7^T p;
- the physical meaning/value of g;
- the overall event rate / seconds.

So this is a role theorem, not a fundamental source theorem.

## 9. Next atom

### MAXIMUM-ENTROPY-REFRESH-54

Test whether complete refresh can be selected information-theoretically:

Among one-label transition kernels with a fixed target/stationary distribution
q, maximize the conditional entropy of the new label given the old label.

If the unique maximizer is the independent reset kernel

    K_ij=q_j,

then isotropic microscopic forgetting gains a target-blind naturality
principle rather than being merely postulated.
