# OBS intake: proof completions and corrections

Date: 2026-10-15. Status: analytical model audit, not mechanized proof.
Primary law: the supplied labelled unit-attempt LOO heat bath; exact spectral
critical gain, fixed bounded lambda, sufficiently large N through multiples
of 12. Inherited MT/UF premises are bound in [verification.json](verification.json).

## 1. Correct the dependency, not its historical custody

`source_context/MT_INTAKE_USER.md` is **not** the accepted MT intake. Its
normalized content is the MT *research plan*, including its unexecuted-task
status. This is a substantive misidentification, not a Markdown-only
difference from the intake. The source reports' claims that this attachment
contains the accepted intake must not be inherited.

The real local [MT intake](../fin_critical_metastability_review_20261015/INTAKE_20261015.md)
and its [proof completions](../fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md)
are present and hash-bound here. They supply the needed theorem. This repairs
the *current dependency chain*, not the absent proof of what the executor
actually read. The preserved MT producer ZIP is not the later review.

The attached OBS plan agrees after the recorded formatting transformations,
including the observed single-emphasis artifact `sqrt(N)<<k*<<N*`.
Do not generalize that normalization to deletion of scientific symbols.

## 2. Exchangeable lift and genuine future paths

For each count n let U_n be the uniform law on its labelled orbit. Generator
covariance and an exchangeable initial law give that conditional law exactly.
The disjoint supports imply

    TV(sum nu(n) U_n, sum eta(n) U_n) = TV(nu,eta).

At a deterministic t_N in the MT window, couple the entire labelled law to
the lift of pi_D,N, then use the *same full labelled* transition law. Entire
future-path TV is at most the initial coupling error. No independence across
times is inferred from the one-time orbit statement.

For the stationary reflected comparison, suppression affects only jumps
leaving the count domain. Its lifted Gibbs law is invariant by detailed
balance, whether or not one separately proves every labelled orbit communicates.
The stationary expected number of suppressed attempts in sqrt(N)S is at
most a polynomial times exp(-Nb), using the real boundary layer and the
local partition lower bound. Couple this process to the full parent until
the first suppressed jump. This pays the boundary term without pretending
the full-generator martingale is exact for the reflected chain.

MT supplies all fixed polynomial moments of Z0 and sqrt(N)w0 under pi_D,N.
UF's energy proof therefore extends to that random initial law, with initial
large-radius events explicitly retained. Compact soft stopping, negligible
fast exit and then nonexplosion localization give the stationary quartic
path law on each fixed scaled horizon. This is not substitution of a random
initial law into a theorem whose only stated start was Z0=0 without checking
its proof.

Transfer to the actual iid-start process is in path distribution/TV coupling.
An unspecified o(1) TV error alone does **not** transfer unbounded scaled
moments, entropy, or variances.

## 3. Pay the conditional finite-population CLT directly

Condition on a parity population of mean m. Draw k signs X_r without replacement
and put T_r=sum_(l<=r)(X_l-m). For r=0,...,k<N the Doob martingale of T_k is

    D_r = (N-k)/(N-r) T_r.

Its increment is `(N-k)/(N-r) * (X_r-E[X_r|F_(r-1)])`. Dividing by sqrt(k)
gives jumps at most 2/sqrt(k). If k/N->0 and m->0, the coefficient tends to
one uniformly. The remaining-population mean differs from m by at most
2k/(N-k), so the sum of conditional increment variances divided by k tends
uniformly to one. The bounded-jump martingale CLT gives the required normal
limit, uniformly when |N^(1/4)m|<=K. MT tightness removes the central cutoff.
This direct payment does not rely on an unread formulation of Hajek's theorem.

For two disjoint critical blocks, after the first draw the remaining mean
changes by `-k1(m_B1-m)/(N-k1)`. Its effect on the standardized second-block
mean is O_P(sqrt(k1*k2)/N)=o_P(1). Apply the same conditional CLT to the second
block. This establishes joint independent Gaussian errors in the limit;
zero limiting covariance alone would not prove their joint normality.

Thus for k/sqrt(N)->alpha>0 the joint snapshot channel is

    (Z_N,N^(1/4)m_B) => (Z,Z+G/sqrt(alpha)),

where Z has the local quartic prior and G is an independent standard normal.
Two readouts share Z; they are not unconditionally independent.

For k=o(sqrt(N)), the source coupling/Pinsker bounds establish the precise
statement

    TV(Law(Z_N,signs_B), Law(Z_N) x uniform_signs_k) -> 0.

They do not assert TV convergence of the lattice-valued Z_N to its continuous
limit. This is a parity-snapshot statement, not automatically an audit of
every conceivable full-label readout.

For sqrt(N)<<k<<N, conditional variance gives

    E_piD[(N^(1/4)(m_B-m))^2] <= C sqrt(N)/k -> 0.

Under the actual post-burn-in iid law this transfers in probability, not
automatically in L2. For the bounded sign rule, decision boundaries Z=0 and
Z+G/sqrt(alpha)=0 have zero limiting probability, so the risk converges.
Bayes optimality is proved for the *limiting* symmetric-prior Gaussian
experiment, not for every finite-N instrument.

## 4. Raw-path no-go requires a distribution, not a variance assertion

The source private-contrast identity and brackets pay the stationary OU
limit on each fixed interval of model time. At any fixed interior scaled
time s, with h_N=N^(-1/2), the actual post-window readout increment converges
in law to

    Normal(0, 2(1-exp(-1))/alpha).

The collective contribution vanishes on that small interval. If a J1-tight
subsequence existed, choose a deterministic s at which its cadlag limit is
almost surely continuous. Evaluation at s and s+h_N would force increments
to vanish in probability, contradicting the displayed nondegenerate limit.

The displayed number is the **variance of the limiting increment**. The
source's `lim Var(...)` wording for the actual iid-start law additionally
needs moment uniform integrability, which the bare MT TV error does not
supply. That extra finite-N variance convergence is not admitted here and
is unnecessary for the no-go.

## 5. Joint record martingale and initial independence

For fixed k let c_a=1{a in B}/k-1/N and v=m_B-m. With the actual LOO refresh,

    L v=-v+e_N,  |e_N|<=2 |g_N| max|A_ij|/N.

Indeed the gradient L1 norm of a softmax parity mean is `1-q^2<=1`, and
two deleted-copy fields differ in sup norm by at most the displayed bound.
The integrated endpoint payment is at most 4N^(-1/4), and the LOO forcing
payment is at most C S N^(-3/4). No fast-mode averaging assumption is needed
for this uniform heat-bath identity.

Use h_a=2(1-S_a q_(xi_a)) in the source's exact three bracket formulas. Globally,

    Gamma(m)<=4/N,
    Gamma(v)<=4(1/k-1/N),
    |Gamma(m,v)|<=8/N.

The last bound uses sum|c_a|<=2. Consequently the rescaled cross bracket is
at most 8S/sqrt(N). Occupation convergence to u pays the two diagonal limits
2s and 2s/k. Jumps of the two rescaled martingales are at most
2N^(-3/4) and 2N^(-1/4)/k, respectively.

For each scalar linear combination, bounded predictable bracket and the
bound on jumps imply a uniform fourth-moment maximum bound on finite
horizons. One can use the fourth-power jump formula or the jump BDG bound:
the fourth-jump compensator is at most maximal-jump-squared times the
quadratic compensator. This supplies UI of both the martingale and its
square-minus-bracket, not just an L2 stochastic-boundedness statement.

Pass their identities jointly with Z0 and the soft drift in the *natural
past filtration*. [Whitt Theorem 2.1(ii) and Theorem 5.3 with its proof](https://arxiv.org/html/0712.1929v2)
were checked for this use. The predictable brackets here are continuous;
their jump condition is automatic. Vanishing martingale jumps and the paid
UI give a two-dimensional Brownian motion in the joint limiting filtration,
with covariance diag(2,2/k). Its coordinates are independent and its future
increments are independent of Z0. Cubic dissipativity and local Lipschitzness
give unique, nonexplosive state dynamics. Thus the joint limit is

    dZ=(chi lambda Z-Z^3/3)ds+sqrt(2)dB,
    dY=Z ds+sqrt(2/k)dW,
    Z0~rho_lambda, Y0=0, B independent of W and Z0.

Z and Y themselves are not independent. No quantitative finite-N process-TV
rate, physical detector or physical stored memory is obtained.

## 6. Correct the exact helper, without rewriting the archive

The source `joint_hypergeom_pmf(12,2,4,3,3,0)` raises ValueError. The first
block asks for three plus signs from a population containing only two, so
the correct answer is zero. The implementation enters a nonexistent second
conditional population after a zero first factor. The formula is correct;
the general-purpose implementation is incomplete at its support boundary.

[observation_kernels.py](observation_kernels.py) is the separately versioned
exact-Fraction correction. It validates the contract and returns zero before
forming that nonexistent conditional law. The checker tests edge populations,
empty second blocks, exhaustive finite supports and a separate multinomial
allocation identity. Historical source/check results remain unchanged.

## 7. What the weak record theorem leaves unpaid

Joint weak convergence does not imply convergence of filters or information.
For example, let C be an equiprobable sign and Y_N=C/N. Then
`(C,Y_N)=>(C,0)`, yet C is exactly recoverable from each Y_N and not from
its limit. This is a logical counterexample, not a demonstrated FIN defect.

In FIN, the entire ideal integral path retains its right derivative, hence
the block parity history away from a null fixed-terminal-jump event. An
amplitude/path limit can hide fine observational structure. Whether the
derived diffusion experiment is predictively sufficient requires its own
likelihood argument, or a fixed finite-resolution observation contract.
It cannot be supplied merely by naming W an innovation.
