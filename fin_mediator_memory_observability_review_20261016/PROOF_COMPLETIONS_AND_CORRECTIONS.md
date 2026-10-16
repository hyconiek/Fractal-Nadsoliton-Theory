# MM proof audit: actual fast process, labelled records and reference gates

Date:2026-10-16. Analytical model arguments under inherited exact spectral,
MT local-domain and stationary-moment premises; not mechanized or peer reviewed.
Primary FIN is unchanged. The mixed OU sampler remains a separate diagnostic.

## 1. Actual inherited inputs, not retrospective custody

The real HM completion supplies the conservative reversible diagnostic,
its stationary count marginal pi_D and conditional normal rho_n. MT supplies
all fixed polynomial moments and pi_D(B_D)<=poly(N)exp(-bN). KO supplies
the joint conditional finite-population CLT; DKO supplies the native fast
process and boundary/UI method. These actual local files are hash-bound in
verification.json. Their availability now does not prove historical access
to the producer's missing wrappers/reviews.

The new source's generic three-pole factor is useful but is not a discovery
of the projection-memory principle or old report298's spectral curvature.
The new contribution is the diagnostic fast-readout theorem and a constraint
which removes the private observation pole. The two kappas are {1,4}, not
the single fraction1/4 suggested by ambiguous source shorthand.

## 2. Full labelled characteristics and stationary initialization

The instrument reads a predetermined block's label counts, not theta and
not merely the total count vector. Let alpha=k/N and
a_l=N/k-1 inside the block, -1 outside. For r=1,3 put

    F_r=N^(-1/2)sum_l f_r(x_l), U_r=N^(-1/2)sum_l a_l f_r(x_l).

The diagnostic uses the same q(theta)=softmax(X theta) for each label slot,
so off the departure layer LF_r=sqrt(N)f_r.q-F_r and LU_r=-U_r exactly.
Every joint count bracket is the source's weighted finite jump sum. Cross
jump/field brackets vanish exactly; count cross brackets do not vanish
statewise. Our repeated-label, nonuniform block and nonzero-field fixture
checks all sixteen entries, not only equilibrium diagonal values.

Set c=Lambda3/Lambda6, gamma=1-c and
H=sqrt(N)*sqrt(Lambda3/12)*theta_3. At the exact g_c and lambda0,

    LH=kappa(cF_3-H), Gamma(H)=2kappa c.

Under Pi_D=pi_D rho_n, H=cF_3+sqrt(c)G exactly at this gain, with G a
standard normal independent of counts. MT's global fast Gaussian law gives
Cov(F_3,H)=gamma^(-1)[[1,c],[c,c]]. Conditional on counts, the block is
hypergeometric independently of this Gaussian residual. KO's conditional
characteristic-function limit gives private normal variance
c_alpha=(1-alpha)/alpha independently of global/field/quartic variables.
The null feature is A7-null, not null for the unprojected Laplacian. Its
global variance is1; its private variance is c_alpha. Product Gaussian
initialization follows from these conditional limits, not from covariance alone.

## 3. Integrated errors, maximum tightness and joint filtration

Write X theta=h sigma+xi, with xi transverse. For every fixed moment order,
HM/MT give E|h|^q=O(N^(-q/4)) and E||xi||^q=O(N^(-q/2)). Both f1 and f3
have zero softmax expectation on the entire parity axis. The first
transverse derivative at h=0 is H/sqrt(N) for f3 and zero for f1.
Bounded softmax derivatives therefore bound the integrated scaled remainder by

    C_T sqrt(N) E[|h| ||xi||+||xi||^2]=O_T(N^(-1/4)+N^(-1/2)).

This global derivative argument includes the unbounded Gaussian field; a
pointwise small-field assertion is not needed. At fixed alpha>0, conditional
hypergeometric moments, p->u in moments and q(theta)->u similarly give
integrated bracket errors tending to zero. For example the nonzero
weighted-old-label terms cancel after conditional block averaging, with
fluctuations O(N^(-1/2)), not identically in each configuration.

The departure rejection changes F/U drifts by at most
C_alpha sqrt(N)1_B and brackets by at most C_alpha 1_B. Their stationary
integrated expectations vanish by the actual MT exponential upper bound.
Field drift/noise is not reflected. No allowed-label renormalization is used.

Let X_N=(F3,H3,U3,F1,U1). Its exact semimartingale equation has fixed linear
matrix J, residual of vanishing expected total variation and a martingale
with deterministic bracket-density upper bound C_alpha,kappa. Variation
of constants or Gronwall bounds the maximum on [0,T] by tight initial data,
the martingale maximum and the residual variation. BDG gives bounded
martingale second and fourth maxima: jumps are O(N^(-1/2)), their fourth
compensator is bounded by the maximal jump squared times the second
compensator, and the field Brownian bracket is constant. Thus square
martingales are uniformly integrable, not merely bounded in expectation.
Localized drift increments give Aldous tightness, and vanishing jumps
make each limit continuous.

Pass the martingale and square/cross-martingale identities jointly with
the initial variables in their natural limiting filtration. The bracket is
diag(2,2kappa c,2c_alpha,2,2c_alpha)t. Levy characterization gives independent
Brownian coordinates with future increments independent of the initial
variables. The linear SDE has a unique law; the covariance above is its
stationary covariance. In particular,

    dF3=(-F3+H3)dt+sqrt(2)dB,
    dH3=kappa(cF3-H3)dt+sqrt(2kappa c)dW,
    dU3=-U3dt+sqrt(2c_alpha)dV.

F1,U1 are independent stationary rate-one OU processes, independent of
this active system. The count parity scaled by N^(1/4) is frozen on this
fast horizon: its drift expectation is at most O(N^(-1/4)) from h-m and
the nonlinear field terms; its scaled bracket is O(N^(-1/2)). Doob and
the integrated drift bound pay the maximum frozen-coordinate statement.

This completes the source's abbreviated step5. Scope: stationary reflected
diagnostic, fixed alpha in(0,1), fixed T, lambda0, N in12N, kappa1 or4.
For alpha1 the private variable vanishes and is omitted. Native late-iid
bounded-path transfer uses its own existing coupling; no diagnostic iid,
unreflected, shrinking-block or growing-horizon theorem is imported.

## 4. Filtered determinant and independent arithmetic

For C_j=A x^j+B y^j+D q^j,

    d_j=C_(j+1)-q C_j,
    d1*d3-d2^2=A B x y (x-q)(y-q)(x-y)^2.

For the diagnostic beta_s<gamma<1<beta_f, so x>q>y>0 and both A,B>0;
the determinant is strictly negative. The native q-filter leaves one
exponential and gives zero even with the private pole. Ordinary unfiltered
native curvature is positive, so is an essential false-positive control.
Normalization by C0 changes magnitude, not sign.

The Gaussian arcsin identity converts seven bounded ternary sign-pair means
into these normalized correlations. The independent certificate freshly
binds exact-decimal spectral values to inherited rational boxes and evaluates
the factored expression, rather than subtracting two nearly equal interval
products as the source does. Public decimals are rounded outward by exact
Fractions of raw dyadic interval endpoints. It retains the mpmath.iv
correctness assumption and is not an IEEE/libm or proof-assistant certificate.

The original gamma box is now paid locally. The very long producer decimal
strings and float previews are not all endorsed as outward endpoints.
Ordinary binary64 matrix-exponential diagnostics differ from the narrow
certificate by up to about3.35e-16 and need not lie inside it: the verifier
records that failed inclusion explicitly and tests numerical agreement
separately. Those point values do not certify arithmetic.

## 5. Correct MM-004: score convergence is not exact-gate convergence

The ten-sign law on {-1,0,+1}^10 converges in TV on this fixed finite
alphabet, because its Gaussian coordinate boundaries have probability zero.
Finite zero outcomes are retained, not dropped. The seven means, sine
transforms and ratio q_N=R_null,2/R_null,1 converge, with denominator tending
to exp(-1)>0 and q_N eventually in(0,1). Thus the **raw** determinant H_N
converges to0 or the certified negative value. The reference residual

    Delta_N=R_null,3-R_null,1*q_N^2

also tends to zero. A normalized single-rate reference additionally satisfies
R_null,1=q; that residual may be retained as well. A finite check of a few
lags does not identify an arbitrary instrument's full noise law.

However, the source JSON's exact test Delta_N==0 is discontinuous. For a
valid stationary Gaussian reference with correlations
R_j=(1-epsilon)q^j+epsilon r^j, q!=r in(0,1),

    R1*R3-R2^2=epsilon(1-epsilon)q r(q-r)^2>0.

As epsilon->0 these laws approach the desired reference, while every exact
equality gate fails. Therefore do not claim that MM's boolean guarded
selection becomes supported at finite N. We admit the joint population-law
limit (H_N,Delta_N), with an open denominator guard, and the ideal-law
constraint, not an operational finite-N acceptance rule. A tolerance sequence
exists abstractly after convergence, but its prescribed rate, uniform error
and finite-data meaning are unpaid. No guessed threshold is introduced.

Distinct limiting finite-alphabet laws imply a strictly positive asymptotic
TV separation for each fixed comparison. This does not supply its useful
numeric TV lower bound, N0, sample size or dependent-data confidence rule.

## 6. Repair the reset control at the labelled level

Equality of total count paths alone does not imply equality of block
histories. The source MM-004/NEGATIVE_CONTROLS.md item3 makes that invalid
inference. For instance (x1,x2)=(0,1) and (1,0) have equal counts and
different first-slot f3 observations. A label-swap process can change such
histories without moving counts at all.

A correct auxiliary lift exists explicitly. Let Q_D^lab(x,x') be the
original labelled native generator, pi_D^lab its stationary lift, n(x) its
total counts and rho_n the Gaussian posterior. On (x,theta) define

    G_eta f(x,theta)=sum_(x'!=x) Q_D^lab(x,x')
       [integral f(x',theta')rho_n(x')(dtheta')-f(x,theta)]
       +eta[integral f(x,theta')rho_n(x)(dtheta')-f(x,theta)].

The product conditional law pi_D^lab(x)rho_n(x)(dtheta) is reversible:
each off-diagonal edge is paired by labelled detailed balance; reset is
conditionally reversible. On f(x) the generator is exactly Q_D^lab f,
independently of theta and eta. Hence all labelled paths, and thus every
predetermined block record, are preserved at every eta. This **explicit
labelled extension** pays the intended control; it is not a consequence
of the old count-only statement and is not a new adopted FIN law.
The small three-atom residual fixture checks intertwining and reversibility,
not Gaussian discretization or critical trapping at N2.

## 7. Meaning and nonconclusions

Fixed per-channel nonzero real gains leave the sign autocorrelations
unchanged. A shared positive time-unit factor changes all poles consistently;
the zero/negative sign remains, but its magnitude and certified numbers are
not clock-independent. No uniform margin at alpha->0 or q->0,1 follows.

The result distinguishes the licensed native partial-observation class from
the specified feedback class. It does not identify a unique hidden medium,
exclude coloured readout memory, infer nature's kappa, source A7/kinetics,
calibrate an apparatus or produce quantum interference. Negative determinant
is a classical score, not a negative probability. All physical gates remain.
