# KO audit: comparator, finite-population CLT and probability arithmetic

Date: 2026-10-15. Analytical model audit under the inherited exact spectral
and MT/UF/OBS premises. Not mechanized or externally peer reviewed.

## 1. The comparator's local regime is paid without a global phase transfer

For A_epsilon the trace/diagonal/parity identities are exact. For both frozen
operators the largest transverse eigenvalue remains Lambda5, as checked
against the inherited rational spectral boxes. Therefore the source local
Hessian estimate is valid on the same cylinder.

At exact g_c its rational lower bound is `4458/49375>0.09`. Retaining the
earlier extra gain-box reserve gives `8679/98750>0.08`. The older certificate
contains that reserve; the two constants must not be equated without it.
The entropy derivative bound, parity-line potential and transverse first-
derivative cancellation are unchanged. Bounded derivatives, positive interior,
real count moves, root-path energy excess and type-weight estimates thus
remain those in the actual MT proof completion, not merely a matching Hessian.

PSD and the unchanged positive gain pay invariant-weight domination for the
comparator too. The count-graph path proof gives polynomial reflected mixing,
and the same positive boundary height gives the parent exit estimate. The
Gaussian local measure uses the same scaled lattice/tail argument with its
changed fast precision. No global minimum census or B4 theorem is transferred.

A common deterministic origin can exceed both existential polynomial burn-ins
while remaining below both exponential exit windows. No useful numeric N0,
burn-in or laboratory clock follows. The real MT/OBS/LI review completions
are available here and hash-bound, resolving the producer's conditional
input-access gate locally, not retrospectively attesting its access.

## 2. Complete the Fourier-inversion CLT tail estimate

The source's Bernoulli-conditioning formula is correct. To justify its tail,
write a=k/N, centered bounded population vectors b_l, and z_l=t.b_l/sqrt(N).
For a bounded away from0 and1,

    |1-a+a exp(iu)|^2 = 1-4a(1-a)sin^2(u/2)
       <= exp[-c dist(u,2pi Z)^2].

For theta in [-pi,pi], |z_l|<=C_t/sqrt(N). If
`|theta|>=2C_t/sqrt(N)`, distance from theta+z_l to2pi Z is at least
|theta|/2. The product of all N factors is consequently at most
`exp(-c' N theta^2)`. This is a bound on the **perturbed** product; it is
not an unjustified domination by the unperturbed binomial product.

Set x=sqrt(N)theta. On |x|<=2C_t the product modulus is at most1; elsewhere
it is bounded by exp(-c' x^2). This gives an N-independent integrable
dominator for the rescaled integral. For fixed x the logarithmic expansion is

    i a sum_l (x+t.b_l)/sqrt(N)
      -a(1-a)/2 sum_l (x+t.b_l)^2/N + o(1).

The external exp(-ik theta) cancels the first term exactly because k=Na
and sum b_l=0. The quadratic cross term vanishes for the same reason.
The denominator is the same integral at t=0; dominated convergence pays its
positive Gaussian limit. Thus the centered sampled sum divided by sqrt(N)
has covariance a(1-a)Sigma. Multiplication by1/a gives the source's U
covariance `(1-a)Sigma/a`.

This works for every deterministic count sequence approaching u. Applying
the sequential/subsequence criterion to random counts makes the conditional
characteristic functions converge in probability to the same deterministic
Gaussian characteristic function. Bounded conditioning then gives joint
independence from any convergent count-measurable statistic. This pays the
claimed stable Gaussian noise, not only a covariance calculation.

## 3. Population and block normalizations

The exact conditional sampling covariance is
`N(N-k)/(k(N-1)) Sigma_n`. At k=N/2 its limit is I for the declared normalized
features. Their Euclidean norm squared is12, so MT precision
`12-g_c Lambda_r` gives global variance `12/(12-g_c Lambda_r)=1/gamma_r`.
The A7 harmonic1 eigenvalue is zero; the original Laplacian harmonic1 value
must not be inserted there. The null observed variance is therefore2.

The full Gaussian count law and the conditional sampling CLT give the joint
snapshot theorem. Actual iid-start event probabilities transfer by the
exchangeable labelled TV lift. Unbounded moment convergence is admitted for
the local comparison only unless separately paid for the actual parent.

## 4. Ideal bounded contrast: validated public arithmetic

The event is fixed in the repository plan, so its model probabilities are
not tuned by this review. Its continuity boundary has zero Gaussian mass.
The strict sign is analytical: larger Lambda3 at fixed Lambda6 increases
`v=1+1/gamma3` and `p=erfc(1/sqrt(2v))`.

[certify_event.py](certify_event.py) independently uses the accepted rational
spectral boxes, freshly binds the exact-decimal kernel by an interval DFT,
and evaluates erf by a **positive**, rather than alternating, series:

    erf(x)=2/sqrt(pi) exp(-x^2)
             sum_(n>=0) 2^n x^(2n+1)/(2n+1)!!.

It uses18 terms. If T_n is the nth term, `T_(n+1)/T_n=2x^2/(2n+3)`.
For the enclosing x<1/3, the omitted positive tail is bounded by
`T_18/(1-r_18)`. Interval operations enclose all retained terms and this
tail. Public decimal endpoints are rounded outward from **exact raw dyadic
endpoints** using Fraction floor/ceiling, not ordinary endpoint printing.

The certificate assumes correctness of mpmath.iv and the inherited input
binding; it is not proof-assistant formalization. Its intervals are wider
than the source's65-digit printout because they use inherited rational input
boxes. No last-digit validity of every source printed endpoint is inferred.

The source's lower-gap bookkeeping uses `gamma4+.01` where the eigenvalue
shift implies `gamma4+.01/Lambda6`. Its reported minimum still equals gamma5;
the independent maximum-mode check pays the correct comparison, so the main
result is unaffected. Do not reuse the incorrect intermediate formula.

## 5. Scope corrections and remaining gates

The preserved drift and bracket are not protected by D12 symmetry **alone**.
Their numerical values also use the common parity restriction, entropy,
critical normalization, microscopic heat bath and unit attempt clock. They
are family-preserved under those joint premises, not constants of nature.

The common finite quantized parity law gives only the stated finite-alphabet
model-ID no-go. The optional stronger full-parity model-index likelihood
argument was not executed; OPEN_CROSS_MODEL remains open. LI's within-model
future-sign risk equality cannot replace it.

KO-004 is genuinely deferred: no fast-time two-OU law or dynamic kernel
contrast is admitted. Static success alone does not identify a process.
The source's physical-input activation gate is correct for empirical
promotion. It is not a theorem that no new consequence of the existing
microscopic law can be studied without an external experiment.
