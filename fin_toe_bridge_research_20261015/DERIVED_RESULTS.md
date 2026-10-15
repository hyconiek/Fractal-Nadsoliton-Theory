# New bridge calculations and conditional analytical propositions

Date:2026-10-15. Status: **RESEARCH_SUBMISSION**, not an accepted physical
bridge or independent peer review. Historical source files are unchanged.
Arithmetic certificate: [verification.json](verification.json).

## 1. Three operators, not one physical object

Let W12 be the supplied positive q12 shell matrix with zero diagonal and

    w_d=cos(.18575 d+.1625)/(1+d^(9/5)), d=1,...,6.

Define L12=diag(W12 1)-W12, with Fourier eigenvalues Lambda_k. All six
weights are positive at this finite q. The admitted A7 instead retains only
harmonics3,4,5,6 of L12. It is a covariance/interaction matrix in its current
FIN role, not automatically an electrical stiffness or spatial propagator.

L_cycle denotes the ordinary nearest-neighbour C12 Laplacian. Its eigenvalues
are mu_k=4 sin^2(pi k/12). The internal cyclic label order is supplied
mathematical structure, not physical position. Infinite absolute-kernel,
signed-kernel and q-refinement models are different objects; none is silently
identified with this finite positive W12.

## 2. Scalar local-string completion: an explicit obstruction

For a Bernstein function

    f(s)=a+bs+integral_(0,infinity)(1-exp(-st))nu(dt), a,b>=0, nu>=0,

the derivative is nonnegative and the function is concave. A nonconstant
such function is strictly increasing on positive s. Complete Bernstein
functions are a subclass. Their local-string extension is characterized in
[Kwaśnicki–Mucha, Theorem3.1 and §2.1](https://arxiv.org/pdf/1707.02475).
For a finite graph the scalar separated extension reduces mode by mode to
the same one-dimensional string ODE; this is a finite spectral-calculus
application, not a claim that the paper directly studies FIN.

### 2.1 A7 itself

A7 has eigenvalue0 at the distinct positive mu_1,mu_2 and a positive
eigenvalue at mu_3. Thus A7 cannot equal f(L_cycle) for any Bernstein f.
Adding a common mass aI or rescaling by a positive clock/weight does not
remove the forbidden plateau. This is an exact support argument, not a
floating rank test. It does not exclude constrained boundary observables,
matrix-valued fields or more general auxiliary structures.

### 2.2 The unprojected L12 does not repair this class

The three upper base modes have mu_4=3,mu_5=2+sqrt(3),mu_6=4. Concavity
requires the successive slope difference

    D=(Lambda6-Lambda5)/(2-sqrt(3))
                    -(Lambda5-Lambda4)/(sqrt(3)-1)

to be <=0. Independent exact-rational interval propagation from the inherited
spectral boxes gives

    .027339348741403146 <= D <= .027339348741607070.

Hence L12 is not f(L_cycle) for any Bernstein f either. The certificate
checks the rational sqrt(3) enclosure by squaring, binds the exact-decimal
kernel back into the inherited spectral boxes, and pays public rounding.
The fresh transcendental binding assumes correctness of mpmath.iv; exact
rational checking does not mechanize that backend.

This excludes **scalar positive separated local extensions over the fixed
cycle base**, not all local hidden realizations, all possible base operators,
or every completion of the strict kernel. It is distinct from the older
legacy-to-strict mode-order obstruction and clock-subordination rigidity.
No legacy role, selector or physical scale transfers.

## 3. Circular-planar completion: a recovered small certificate

This question is **not new**: it was REF-004 in FIN_POST_PHY_TASK_CATALOG.md.
The later ledger records a terminal disposition with unavailable original
per-task bytes; fin_post_phy_review does not independently certify that old
missing result. Here one explicit independently replayable witness is exported.
This is evidence restoration for an old question, not a new campaign to scan
thousands of minors.

In the response convention M=-L12, a circular-planar network with positive
conductances and every port on a disk boundary has nonnegative circular
minors. The exact ordering and convention are in
[Alman–Lian–Tran, Definitions2.2.1/2.2.4 and Theorem2.2.6](https://arxiv.org/pdf/1309.2697).
For P=(0,1),Q=(6,4), the boundary order is p1,p2,q2,q1=(0,1,4,6), and

    det M[P,Q]=w6 w3-w4 w5.

Fresh exact-decimal interval input and rational endpoint arithmetic give

    -.000122681888983639 <= det M[P,Q] <= -.000122681888983638.

Thus the natural-order circular-planar passive realization fails. This does
not exclude a different port order, multiple boundary components, nonplanar
graphs, constrained/vector media or other notions of locality. It certainly
does not establish a physical dimension. Uniqueness of conductances for a
**given critical network** is not selection of a unique microscopic graph.

The original Curtis–Ingerman–Morrow PDF scan was not read successfully;
the precise readable statement used here is the primary Alman–Lian–Tran
paper, which cites it. No original-proof reading is invented.

## 4. Actual interacting histories supply temporal reflection positivity

This is the positive constructive result. Use the actual stationary reversible
finite FIN chain, or its already admitted stationary local reflected comparison.
Denote its full Markov state by X, invariant law by pi, and time reflection
by vartheta X(t)=X(-t). For a bounded functional F of finitely many future
times, the stationary Markov property gives conditional independence of
past and future given X(0). Reversibility identifies their conditional laws.
Consequently

    E[conjugate(F(vartheta X)) F(X)]
                 =sum_x pi(x)|E[F(X)|X(0)=x]|^2 >=0.

The mixed-cylinder version is a Gram form, so the assertion holds for all
finite linear combinations of bounded future cylinders, not just a Gaussian
two-point covariance. Quotient by zero conditional expectation identifies
the reconstructed temporal Hilbert space with L2(pi): functions of X(0)
already span it. Time translation induces the original reversible semigroup.

In row-generator convention,

    H_N=-D_pi^(1/2) L_N D_pi^(-1/2)>=0,
    H_N sqrt(pi)=0.

This resembles the known RK/stochastic-matrix construction; the new FIN
object is the **actual interacting hierarchy and its critical reduction**,
not a claim to have discovered the general mapping. Henley's
[§§3.4–4](https://arxiv.org/pdf/cond-mat/0311345) separates probabilities
from amplitudes explicitly.

Exact reflection positivity belongs to the stationary law. An actual late
iid-start parent inherits bounded finite-window predictions through the paid
MT/OBS path coupling, not exact stationarity or a new physical wall. Temporal
reflection positivity alone is not the full spatial/relativistic OS axioms.

## 5. Reversibility lifts the admitted path limit to strong slow transfer

This is a conditional analytical consequence of already accepted premises,
not a numerical large-N assertion. Put

    H_N=L2(pi_D,N), H=L2(rho_lambda),
    Z_N=N^(1/4) sigma.p, J_N f=f(Z_N),
    T_N(t)=exp(sqrt(N)t L_D,N).

The stationary reflected soft-path limit is paid in the OBS proof completion,
§2, by MT moments and the UF random-initial-law extension. It is not merely
UF's original deterministic Z0=0 theorem. In particular, for real f,g in C_b,

    <J_N g,T_N(t)J_N f>_N -> <g,T(t)f>_rho,
    (Z_N)_# pi_D,N -> rho_lambda.

The quartic diffusion has a bounded-continuous Feller transition action;
local Lipschitz coefficients plus localization/nonexplosion justify continuity
in the initial state. Do not assume uniform-norm strong continuity on all
of C_b. For fixed t, T(t)f is bounded continuous. Now expand exactly:

    ||T_N(t)J_N f-J_N T(t)f||_N^2
      =<J_N f,T_N(2t)J_N f>_N
       -2<T_N(t)J_N f,J_N T(t)f>_N+||J_N T(t)f||_N^2 ->0.

Each term follows from the corresponding stationary two-time law or marginal
law; limiting reversibility identifies the first with ||T(t)f||^2. For complex
functions use conjugation and the real part of the cross term. Neither exact
finite-N Markov closure of Z_N nor independence of hidden modes is assumed.

For r>0, integrate with exp(-rt). Contraction gives domination by
2||f||_infinity exp(-rt), proving the analogous strong **lifted resolvent**
convergence. This is already a substantive improvement over generator balances.
The varying-Hilbert/Mosco interpretation requires the exact dense-section
definitions; [Kuwae–Shioya, §2 and Theorem2.4](https://www.intlpress.com/api/bgcloud-front/resource/pdf/volume/1805786134762446850-1805786134762446850-e776617e92e2a310c0026f43fb47be64.pdf)
is the appropriate framework. The next plan explicitly audits that interface.

No full-space common-Hilbert norm, compact convergence, every eigenvalue,
spectral gap, zero-parameter resolvent, trace, uniform-in-time error or numeric
N0 is inferred. No stationary operator theorem transfers to an actual iid
law merely because its bounded path predictions converge.

## 6. The critical transfer operator is sextic, not an inserted quantum law

Let a=chi lambda, U(z)=z^4/12-a z^2/2, rho_lambda proportional to exp(-U).
The accepted diffusion generator is

    G=partial_z^2-U'(z) partial_z.

The unitary change of Hilbert representation f->sqrt(rho)f gives

    H=-partial_z^2+U'^2/4-U''/2
     =-partial_z^2+z^6/36-a z^4/6+(a^2/4-1/2)z^2+a/2.

Use the nonnegative closed-form/Friedrichs realization; the factorization
Q^*Q with Q=partial_z+U'/2 specifies its positivity. At lambda0,

    Hcrit=-partial_z^2+z^6/36-z^2/2,
    psi0 proportional to exp(-z^4/24), Hcrit psi0=0.

The symbolic identities were independently checked. This is a derived
interacting0+1-dimensional transfer representation. The coordinate z is
collective amplitude, not spatial position. The leading normalized operator
is shared by the admitted matched critical family; its numerical spectrum
would not by itself select A7 or particles.

Physical real-time replacement exp(-tH)->exp(-itH) changes predictions.
Even q12 uniform reset gives a return at t=2pi for the formal unitary but
not for the stochastic process. This is a known boundary, included as a
negative control rather than counted as a new quantum theorem. No physical
complex phase, Born instrument, conserved energy, hbar or Lorentzian time
has been derived.

## 7. Scientific status

The new scalar-completion obstruction is a scoped exact-input certificate.
The circular-planar witness restores evidence for REF-004. The temporal,
lifted-transfer and sextic propositions are analytical deductions under the
hash-bound stationary/model premises, submitted for independent integration.
The ten computational checks are not proofs of the all-N path theorem or
formal verification of the new analysis. AGENTS is not changed by this research.
