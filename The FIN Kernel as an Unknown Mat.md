The FIN Kernel as an Unknown Mathematical
Object
An independent referee report on the classification, spectral structure,
and adaptive-information interpretation of the operators
K legacy_ont ​
 and Kstrict_gate
 ​
of the Fractal-Nadsoliton-Theory repository
Research-level referee report with independent numerical verification
Prepared on the basis of repository state HEAD (2026-07-17)
July 2026
Contents
Abstract and verdict table
1 Provenance, exact definitions, and the convention trap
2 Part 1 — Identification of the mathematical object
3 Part 2 — Verification and falsification of the previous classification
4 Part 3 — Can the operator become an adaptive learning operator?
5 Part 4 — The self-modifying kernel on a dynamic graph
6 Part 5 — Deep spectral analysis
7 Part 6 — The mathematically smallest extension
8 Part 7 — Literature map: has FIN rediscovered a known operator?
9 Part 8 — What is genuinely novel
10 Part 9 — A publishable contribution without physics
11 Part 10 — Blind classification: top-10 fields, ranked
12 Falsification log, counterarguments, and limitations
Appendix A — Reproducibility of all numerical claims
References
Abstract and verdict table
Scope. This report analyzes the mathematical core of the Fractal-Nadsoliton-Theory (FIN) repository — two
frozen kernel functions Klegacy_ont​
 and Kstrict_gate ​
 used as coupling operators on a 12-node cycle and on
random Euclidean point clouds — stripped of all physical interpretation. The guiding question is whether the
kernel is more naturally interpreted as an adaptive information-propagation operator than as a physical field
theory. All spectral and numerical claims in the repository's own audits, and in a previous independent
analysis (repo file Analiza_operatora_jadra_FIN.md , 2026-07-16), were re-derived from scratch;
nothing was taken on trust.
Principal results. (i) On its operational domain Z12 ​
, Kstrict_gate ​
 is, to a Frobenius residual of 2.9%, a
normal-ordered discrete Yukawa (screened-Poisson) Green function N [ A(L + m 2 I) −1 ] of the cycle
Laplacian, with A = 2.053, m2 = 0.749 — not a "smoothed Helmholtz operator," as the previous
analysis conjectured, and not "the operator rather than its Green function," as that analysis concluded; the
correct statement is the reverse. (ii) Klegacy_ont​
 is a band-pass spectral-shell filter — a smeared, non-
idempotent analogue of the spectral measure δ(μ − μ⋆ ) ​
 of the Laplacian — confirming, with one
correction, the "on-shell" reading of the previous analysis. (iii) Neither kernel is a reproducing (Mercer)
kernel, a Gaussian-process covariance, a transition operator, or a pure resolvent; the classification battery is
exhaustive at the 7-mode resolution of the problem. (iv) As an adaptive object, the coupled system ˙ ψ​
 =
iK t ​
ψ, ˙ K = η(ψψ⊤ − γK) is well-posed and stable only with decay; the kernel is provably not a fixed
point of self-driven Hebbian learning (correlation of the learned operator with the time-averaged covariance
of the system's own excitations is 1.000; with the kernel itself, −0.07). (v) The four kernel parameters are
not jointly identifiable from the kernel's action: the spectral Fisher matrix has a numerically exact null
2
1 Provenance, exact definitions, and the convention trap
direction, and a manifestly distinct parameter tuple reproduces the operational profile to 6 × 10−5 . (vi) On
random Euclidean clouds the kernel matrices are Euclidean random matrices with Wigner–Dyson level
statistics (⟨r⟩ = 0.51), bulk modes localized on ∼4–5 sites, and — for the strict kernel — a heat-kernel
spectral dimension d s​
 ≈ 3 at intermediate times, recovering the embedding dimension. (vii) FIN has,
unintentionally, rediscovered a known object: the Chung–Yau discrete Green function (screened variant) on
the cycle, resp. a spectral-shell filter; the novelty is not the object but an unstudied combination: a normal-
ordered oscillatory Green kernel used simultaneously as Hamiltonian, as graph filter, and as the slow variable
of a self-referential Hebbian plasticity loop. That combination supports one honest publication; it does not
support a new paradigm.
Table 0. Verdict table. Confidence: H ≥ 95%, M 80–95%, L 50–80%.
Claim
 Verdict
 Confidence
Kstrict
 ​
 =
 normal-ordered
 discrete
 True,
 ∥R∥F / ​
 ∥K∥F ​
 = 0.029 ,
 m 2 =
H (99%)
Yukawa Green function on Z12
 ​
 0.749
Klegacy ​
 = smeared spectral-shell (band- True; corr. with shell filter 0.950; not
H (97%)
pass) filter, on-shell object
 idempotent, hence not a projector
Previous analysis: "Kstrict ​
 ≈ c(z −
 Refuted: rational-cutoff/resolvent+contact fit is
L)e−tL , closer to the operator than to 5×
 better (nRMSE 0.025 vs
 0.121);
 H (97%)
its Green function"
 operator/Green assignment inverted
False in every tested regime (residuals 4.8–6.7,
Repo claim: "K(d) looks like a
repo's own QW-2139; continuum FT sign-
 H (98%)
damped-Helmholtz Green function"
indefinite)
Repo
 claim:
 Hebbian
 learning
"generates"
 the
 kernel
 False as stated: rule learns the drive covariance
H (99%)
( nadsoliton_neural_analysis.p
 (corr 0.9998 with bare cosine), not the kernel
y )
Adaptive
 interpretation
 of
 K
 is
 True, with decay/homeostasis; fixed points are
H (95%)
mathematically coherent
 projector mixtures, never K itself
Four
 parameters
 (ω, ϕ, β, η)
 False:
 Fisher
 matrix
 singular;
 effective
H (96%)
identifiable from kernel action
 dimension ≤ 3 on the domain
FIN unintentionally rediscovers a known
 Yes: Chung–Yau screened discrete Green
H (95%)
operator
 function (strict); spectral-shell filter (legacy)
Physics-free
 publishable
 contribution
 Yes, one paper (operator/bootstrap analysis);
M (85%)
exists
 journal targets in §10
Spectral graph theory / graph signal processing;
Natural home of the bare kernel
 H (95%)
see §11 ranking
1 Provenance, exact definitions, and the convention trap
1.1 What the repository actually contains
The FIN repository (21,348 files, 564 MB) is a monorepo of an automated theory-development program. Its
mathematical content reduces, for the present purpose, to two frozen parametric kernels and a dynamics built on
3
1 Provenance, exact definitions, and the convention trap
them. The definitions below were extracted from the primary sources, not from secondary summaries:
cos(ωd + ϕ)
 π πKlegacy_ont (d) ​
 = αgeo
 ​
 ​
 ,
 αgeo ​
 = 4 ln 2 = 2.772589, ω = 4 ​
 , ϕ = 6 ​
 , β tors ​
 =
1 + βtors ​
 d
from
 TOE_FINAL_DOCUMENTATION_RELEASE_4_4_LEGACY_FULL.tex
 (Eq.
 408
 therein)
 and(1-1)
RELEASE_6_2_SHANNON_VOID_ASYMMETRY_EN_PL.md , §4.1; and
cos(ωd + ϕ)
Kstrict_gate ​
(d) =
 ​
 ,
 ω = 0.18575, ϕ = 0.1625, β = 1.0, η = 1.8,
1 + β dη
from the machine-readable gate artifact report_qw2118_ktotal_spectral_tripartition_gate.json (1-2)
( kernel: omega=0.18575, phi=0.1625, beta=1.0, eta=1.8 ), which also fixes the matrix specification:
n_octaves=12 , diagonal: zero_self_coupling , distance: cyclic_octave_distance . The strict kernel is
therefore the operator Kij ​
 = K strict ​
 ( dij ) ​
 on the cycle Z12 ​
 with
dij ​
 = min (∣i − j∣, 12 − ∣i − j∣) ∈ {1, … , 6},
 Kii ​
 = 0.
The dynamics used across the repository (e.g. QW-593_Information_Unity.py ) is the linear flow
 (1-3)
˙ ψ​
 = i K ψ,
 ψ(t) = eiKt ψ0 ,
 ​
i.e. K plays the role of a Hamiltonian (a continuous-time quantum walk generator), not of a propagator; the(1-4)
propagator is eiKt.
1.2 The convention trap — four different operators hide under one name
A referee's first duty is to note that "the FIN kernel" is not one object. The repository uses, interchangeably and
often silently, at least four distinct operators built from Eq. (1-1)–(1-3):
Table 1-1. The four operator variants found in the repository and their consequences.
Variant
 Where used
 Structure
 Consequences
Real,
 cyclic
 QW-2118 gate, strict-
 Symmetric
 circulant,
 Real spectrum, degenerate doublets; the object
distance, d ≤ 6
 core chain
 diagonalized by DFT
 analyzed here
Real,
 absolute
 nadsoliton_neur
 Symmetric Toeplitz (path
 No circulant symbol; row sums not constant;
distance, d ≤ 11
 al_analysis.py
 graph)
 boundary effects
Non-Hermitian Toeplitz
 iK not anti-Hermitian ⇒ dynamics (1-4) is not
Complex, ei(ωd+ϕ)
 QW-593
on random clouds
 unitary; authors renormalize ψ by hand each step
Real, on random
 Euclidean random matrix
QW-540–593 family
 Spectrum is a random variable; see §6
Euclidean clouds
 (ERM)
Unless stated otherwise, every statement below refers to variant 1 on Z 12 ​
, which is the only one with a frozen,
gate-certified definition. I verified the strict spectrum against the repository's own JSON artifact: all twelve
eigenvalues agree to nine significant digits (λmin ​
 = −0.68187476, λ max​
 = 1.66030728 ; artifact:
report_qw2118…json ).
4
2 Part 1 — Identification of the mathematical object
1.3 Method of this report
Every quantity was recomputed independently: kernels were rebuilt from Eqs. (1-1)–(1-3); spectra, fits, learning
experiments and random-matrix statistics were produced by fresh code (protocols in Appendix A). Where the
previous analysis (repo file Analiza_operatora_jadra_FIN.md ) is quoted, it is treated as a claim to be falsified,
per the task. Confidence levels are attached to every major conclusion (H ≥ 95%, M 80–95%, L 50–80%).
2 Part 1 — Identification of the mathematical object
2.1 Exact spectral symbols
Both kernels are circulants on Z12 ​
 , hence diagonalized by the discrete Fourier basis vm (j) ​
 = e2πimj/12 /
 12
 ​
. Writing the cycle-Laplacian eigenvalues as μ m ​
 = 2 − 2 cos(2πm/12) , the kernel symbols are
5
λm ​
 = 2 ∑ ​
 K(d) cos( 2πmd
 12
 ​
 ) + K(6)(−1)m ,
 m = 0, … , 6,
d=1
with multiplicities 1, 2, 2, 2, 2, 2, 1. Table 2-1 gives the exact symbols (recomputed; they coincide with(2-1)
the previous analysis's table to all printed digits).
Table 2-1. Exact symbols of the two kernels on Z12 ​
, against cycle-Laplacian eigenvalues μm ​
 and wavenumbers km ​
 =
2πm/12.
m
 km
 ​
 μm
 ​
 λm ​
 (strict)
 λm ​
 (legacy)
 sign pattern
0
 0.000
 0.000
 +1.6603
 −11.1741
 + / −
1
 0.524
 0.268
 +0.9062
 +2.0562
 + / +
2
 1.047
 1.000
 +0.0833
 +10.2030
 + / +
3
 1.571
 2.000
 −0.3011
 −3.2072
 − / −
4
 2.094
 3.000
 −0.5393
 −0.2516
 − / −
5
 2.618
 3.732
 −0.6383
 −2.7725
 − / −
6
 3.142
 4.000
 −0.6819
 −0.8819
 − / −
Three structural facts are immediate from the definitions, before any fitting:
5
2 Part 1 — Identification of the mathematical object
Proposition 2.1 (elementary obstructions).
(a) Neither kernel is positive semidefinite: λmin ​
(Kstrict ​
) = −0.6819 , λmin (Klegacy ​
 ​
) = −11.1741.
Hence neither is a Mercer kernel, an RKHS kernel, or a Gaussian-process covariance matrix. The correct
habitat for indefinite kernels is a reproducing-kernel Kreĭn space[26][27].
(b) Klegacy ​
 cannot be a monotone matrix function g(L) of the cycle Laplacian for any monotone g : its
symbol is non-monotone in μ (two sign changes, − + + − − − −). In particular it is not a resolvent
(zI − L)−1 , not a heat kernel e−tL , not any completely monotone function of L.
(c) K strict ​
 is elementwise non-negative on its domain (minimum off-diagonal entry +0.01107): the first
zero of cos(ωd + ϕ) occurs at d⋆ ​
 = ( π
 2 ​
 − ϕ)/ω = 7.58 > 6, outside the operational domain. Hence,
despite spectral indefiniteness, Kstrict ​
 is a legal weighted-graph adjacency matrix, and Perron–Frobenius
applies: λmax ​
 = 1.6603 is simple, with the strictly positive uniform eigenvector v0 ​
 = 1/ 12​
 (constant
row sums). [H, exact]
2.2 The classification battery
I fitted every standard operator family to the two symbols. The honest quality metric is the normalized RMS error
nRMSE = ∥ ^ λ − λ∥ 2 /∥λ ​
 − ˉ λ∥2 ​
 over the seven distinct modes (the previous analysis's "max relative error
8.8%" turns out to be measured relative to ∣λ∣max ​
, which flattens small-mode errors; on the worst small mode
their own fit errs by 77%).
Table 2-2. Family classification battery. Best-fit parameters by least squares; nRMSE and correlation over the 7 distinct modes.
Winners in bold.
Candidate family
 strict: nRMSE /
 legacy: nRMSE /
Interpretation
g(μ)
 corr
 corr
a/(z − μ)
 pure resolvent (Green function of L − z )
 1.000
 −0.91
 0.762
 0.76
a e−tμ
 heat/diffusion kernel
 0.514
 0.97
 0.718
 0.71
a/(μ + z)
 Yukawa-type (shifted resolvent)
 0.593
 0.92
 —
 —
aμs
 fractional power of L
 0.658
 0.78
 —
 —
Helmholtz operator + heat cutoff (prev.
c(z − μ)e−tμ
 0.121
 0.993
 0.999
 0.08
analysis)
c(z − μ)e−tμ
2
 Helmholtz + Gaussian cutoff
 0.221
 0.975
 —
 —
c(z − μ)/(1 + tμ)
 Helmholtz + rational cutoff
 0.025
 0.9997
 —
 —
A/(μ + B) + C
 resolvent + contact (identity) term
 0.025
 0.9997
 —
 —
a e−((μ−b)/s)
2
 smeared spectral shell (band-pass)
 —
 —
 0.763
 0.81
polynomials deg 2–3
 generic spectral filter
 0.195
 0.98
 0.60
 0.80
The two rows in bold are not competitors but the same object: the Möbius identity
c
 1 z − + tμ
 μ
 ​
 =
 A
 μ + B
 ​
 + C,
 A = t
 c
 (z + 1
 t
 ​
 ),
 B =
 1
 t
 ​
 ,
 C = − c
 t
 ​
 ,
(2-2)
6
2 Part 1 — Identification of the mathematical object
shows that "Helmholtz operator with rational cutoff" is "shifted resolvent plus a constant (contact) term." The
best-fit parameters are
A = 2.0528,
 B ≡ m2 = 0.7492,
 C = −1.0896,
 nRMSE = 0.025.
(2-3)
2.3 The strict kernel is a normal-ordered Yukawa Green function
The constant C in Eq. (2-3) is not a nuisance: it is exactly accounted for by the FIN convention of zeroing the
diagonal. On a vertex-transitive graph the diagonal of any resolvent is a constant, Gii ​
 = 1 12​
 ∑m ​
 A/(μm ​
 +
B) = 1.0884 ; the FIN convention Kii ​
 = 0 therefore subtracts precisely a multiple of the identity.
Numerically ∣C∣ = 1.0896 matches A G 00 ​
 = 1.0884 to 0.1%. Hence:
Theorem 2.2 (identification of Kstrict_gate; ​
 verified numerically, H 99%).
On Z12 ​
, with L the combinatorial cycle Laplacian and N the operation of deleting the diagonal (normal
ordering),
2 −1 2 ∥R∥F
 ​
K strict ​
 = N [ A (L + m I) ] + R,
 A = 2.053, m = 0.749,
 ​
 = 0.029,
∥K∥F
 ​
i.e. the strict kernel is, to 97% in Frobenius norm, the normal-ordered discrete screened-Poisson(2-4)
(Yukawa) Green function of the cycle — a screened instance of the Chung–Yau discrete Green function[1].
Equivalently, in spectral space λm ​
 ≈ A/(μm ​
 + m2) − A G 00 . ​
 The residual R is concentrated at large
distances (relative errors 0.1%, 5.6%, 3.8%, 17%, 20%, 27% for d = 1, … , 6; absolute errors ≤ 0.008).
Proof. Direct exact computation: build L; compute G = A(L + m2 I)−1; delete diagonal; compare entrywise
with Eq. (1-2). The real-space profile correlation is 0.9993; see Figure 1(a,b). The identification is unique in the 3-
parameter class (A, B, C): the least-squares problem has a single minimizer (residual surface unimodal over the
box A ∈ [0.5, 5], B ∈ [0.1, 3], C ∈ [−3, 0]). □
This resolves the repository's central open question — "of which operator is K the Green function?"
( Recenzja_adwersarialna…md , M11 roadmap item; RAPORT_QW2139 ) — with a precise answer and, remarkably,
in the affirmative direction for once: on the cycle, Kstrict ​
 is a Green function, namely of the screened-Poisson
operator L + m 2
 I , up to normal ordering. The repository's own test of this hypothesis (QW-2139) returned
residuals of 4.8–6.7 and rejected it, because it tested continuum radial Laplace/Helmholtz ansätze against the raw
kernel — the wrong operator class and the wrong substrate. The correct class is discrete, graph-native resolvents.
(For the continuum R 3 reading, the 3D Fourier–Hankel transform of Kstrict ​
(r) is not of Yukawa form — best-
fit correlation 0.07 — because d−1.8 decay is algebraic, not exponential; the identification is intrinsically discrete.
See §3.3.)
2.4 The legacy kernel is a spectral-shell (band-pass) filter
For Klegacy​
 no differential-operator family fits (Table 2-2, right column). What fits is the opposite kind of object:
the symbol is concentrated on modes m = 1, 2 (λ 1 ​
 = +2.06, λ2 ​
 = +10.20), suppressed elsewhere except
the uniform mode (λ0 ​
 = −11.17 ). Correlations with canonical masks:
Table 2-3. Legacy kernel vs canonical spectral masks (off-diagonal correlation, Z12 ).
 ​
7
2 Part 1 — Identification of the mathematical object
Canonical object
Yukawa resolvent (L + 0.749)−1
Heat kernel e−0.45L
Spectral projector onto {m ≤ 2}
Smeared shell exp[−((μ − 1)/0.5)2 ]
Random-walk transition P = D−1 A
corr with K legacy
 ​
0.471
0.575
0.824
0.950
0.601
corr with Kstrict
0.9993
0.984
0.872
0.307
0.931
Proposition 2.3 (identification of K legacy_ont ; ​
 H 97%).
K legacy ​
 is a band-pass graph filter localized on the spectral shell around k ≈ π/3 (peak at m = 2,
k2 ​
 = π/3 = 1.047 ; the kernel frequency ω = π/4 = 0.785 lies between k1 ​
 and k2 ), ​
 plus a strong
negative uniform-mode component. It is a smeared, non-idempotent analogue of the spectral measure δ(L −
~
 ~
μ⋆ ) ​
 : the normalized idempotency defect mean ∣ λ2 − λ∣ = 0.43 (a true projector has 0). It is therefore an
on-shell object (a resonance/shell selector), not an off-shell propagation operator. In the continuum reading
(§3.2) it is a standing-wave mixture of the principal-value and spectral-measure parts of the 3D Helmholtz
resolvent at q = ω — again on-shell.
2.5 Quantified similarity to every family in the task list
Table 2-4. Quantified family membership. “—” = excluded by Proposition 2.1.
Family
 Kstrict
 ​
 Klegacy
 ​
 Basis
Green
 operator
 (discrete,
0.9993 — member
 0.47 — no
 Thm 2.2
screened)
Reproducing (Mercer) kernel /
 Prop
 2.1(a):
excluded
 excluded
GP covariance
 λmin ​
 < 0
yes — low-pass with notch at μ ≈
 yes — band-pass at shell
Graph filter (GSP sense)
 exact symbols[7]
1.18
 k ≈ π/3
Diffusion operator / heat kernel
 0.98 approx only
 —
 Table 2-2
Helmholtz
 operator
 (z − L) 0.93–0.99 approx; but see §3: it is the
—
 Table 2-2
with cutoff
 resolvent
inverse relationship: K ≈ its Green
Screened-Poisson operator
 —
 Thm 2.2
function
Graph
 wavelet
 /
 spectral
 0.82 crude; shell filter
0.87 crude
 Table 2-3[8]
projectors
 0.95 better
Neural (NTK) / random-feature
excluded (indefinite)
 excluded
 Prop 2.1(a)[34][40]
kernel
Message-passing
 operator
 (1
yes, trivially: ψ ↦ Kψ
 yes, with signed weights
 definition[42]
linear round)
static
 analogue
 only
 (no
 data-
Attention kernel
 no
 §4.4[17]
dependence, no softmax)
Reservoir
 (fixed
 recurrent
 yes, unstable one (ρ =
yes, linear reservoir without readout
 §5.2[20]
operator)
 11.2 )
8
3 Part 2 — Verification and falsification of the previous classification
Integral transform / Nyström
 on point clouds:
on point clouds: exactly
 §6[51]
discretization
 exactly
Spectral
 projector
no (defect 0.29)
 no (defect 0.43)
 Prop 2.3
(idempotent)
Perron–Frobenius / transfer
 shares PF structure (positive entries, uniform PF mode) but
 no (73% negative
 Prop
operator
 is not stochastic (row sums 1.66  =
 1)
 weights)
 2.1(c)
of its own linear flow:
 eiKt is a unitary Koopman
Koopman operator
 same
 §5.4[35]
semigroup on observables ψ
The headline of Part 1: the two kernels belong to two different mathematical genera. Kstrict ​
 is an off-shell
object (a resolvent/Green function — a medium); Klegacy ​
 is an on-shell object (a shell selector — a resonance).
Any discussion that treats "the FIN kernel" as one thing conflates a propagator with a spectral measure. [H]
Figure 1. (a) Kernel profiles on Z12 ​
 with the Yukawa fit of Theorem 2.2. (b) Spectral symbols with family fits: the
resolvent+contact model (nRMSE 0.025) versus the previous analysis's cutoff-Helmholtz model (nRMSE 0.121). (c) Unfolded
level-spacing distribution of the legacy-kernel Euclidean random matrix against Poisson and Wigner–Dyson benchmarks. (d)
Self-driven Oja plasticity: correlation of the evolving operator with the frozen strict kernel decays to ≈ 0.6 — the kernel is not
a fixed point (§4.3).
3 Part 2 — Verification and falsification of the previous classification
The previous analysis (repo file Analiza_operatora_jadra_FIN.md ) concluded: (i) K legacy ​
 "resembles an on-
shell spectral projector"; (ii) Kstrict ​
 "resembles a smoothed Helmholtz operator," specifically c(z − L)e−tL
with z = 1.179, t = 0.456, correlation 0.993, and is "closer to the operator itself than to its Green function
— the Green function would be K −1 , not K." I attempted to falsify both.
9
3 Part 2 — Verification and falsification of the previous classification
3.1 Claim (i): legacy as on-shell projector — CONFIRMED with one correction
The on-shell reading survives every test: band-pass symbol (Table 2-1), failure of all off-shell families (Table 2-2),
0.95 correlation with a smeared shell mask (Table 2-3), and the continuum decomposition below (§3.2). The
correction: it is not a projector. The idempotency defect is 0.43; the object is a filter, not a projection — it
amplifies the shell modes rather than selecting them (ratio λ2 ​
 /λ1 ​
 = 5.0 , not 1). Calling it a "spectral projector"
overstates the structure by one idempotency equation. The precise term is spectral-shell band-pass filter. [H
97%]
3.2 Continuum reading of legacy: standing wave on the mass shell (exact decomposition)
On R3 , with G+ (r) = eiωr /(4πr) the outgoing 3D Helmholtz Green function at frequency ω , the legacy
kernel for r ≫ β −1 = 100 satisfies
Klegacy ​
(r) ∼
 α β
 geo ​
 ​
 cos(ωr r
 + ϕ)
​
 =
 β
 4παgeo
 ​
 ​
 [ cos ϕ Re G+ (r) − sin ϕ Im G+ (r)],
i.e. with ϕ = π/6 : 86.6% principal-value (off-shell reactive) part and 50% spectral-measure (on-shell(3-1)
radiative) part — a standing wave, never the outgoing Green function (which requires the quadrature
cos +i sin with relative phase i; K legacy ​
 is real, so it cannot). Note the envelope: on the whole operational
domain d ≤ 11 the damping 1/(1 + 0.01d) varies by less than 9% (from 0.943 at d = 6 to 0.901 at d =
11 ); the 1/r Green-like envelope exists only asymptotically for r ≳ 100, far outside any domain on which the
repository ever evaluates the kernel. So even the cosmetic 1/r resemblance is physically inert. [H]
3.3 Claim (ii): strict as smoothed Helmholtz operator — REFUTED, and the correct
statement is the reverse
The previous analysis's fit is reproducible (z = 1.1787 , t = 0.4564 , c = 1.3055 , corr 0.9927 — I get the
same numbers). But it is not the best 3-parameter fit in its own neighborhood: the rational cutoff c(z −
μ)/(1 + tμ) achieves nRMSE 0.025 vs 0.121 — a factor 5 — and by the Möbius identity (2-2) that family is
the resolvent-plus-contact family. The conclusion "K is closer to the operator than to its Green function" is
therefore an artifact of truncating the search at exponential cutoffs. With the search extended by one rational
family, the assignment inverts:
Falsification result.Kstrict ​
 is not best described as "a smoothed Helmholtz operator whose Green function
−1would be K ." It is best described as the (normal-ordered) Green function of the screened-Poisson
operator L + m2 I itself (Theorem 2.2). The repository's intuition "K looks like a Green function" was
right all along; the previous analysis's correction of it was wrong; and the repository's own QW-2139
rejection of the Green reading was wrong for a different reason (testing continuum radial ansätze instead of
graph resolvents). [H 97%]
Two honest qualifications. (a) On R3 the identification fails (§2.3): the algebraic tail d−1.8 is not a 3D Yukawa
tail. The phrase "effective dimension deff ​
 = 4.6" from the previous analysis (matching r −(d−1)/2 to r −1.8 ) is
arithmetically correct and physically empty; a more useful continuum reading is the Riesz kernel of the
10
4 Part 3 — Can the operator become an adaptive learning operator?
fractional Laplacian (−Δ)s/2 in R3 with 3 − s = 1.8 ⇒ s = 1.2[38] — but the oscillation breaks the
sign-definiteness that Riesz kernels have, so this, too, is envelope-level only. The Yukawa identification is a
property of the discrete object. (b) The previous analysis's headline ordering — "legacy on-shell, strict off-shell"
— is confirmed and sharpened by Thm 2.2/Prop 2.3; only the operator-vs-Green assignment is reversed.
3.4 Score sheet of Part 2
Table 3-1. Verification score sheet.
Prior claim
 Result
 Correct statement
Confirmed
 w/
Legacy ~ on-shell spectral projector
 on-shell spectral-shell filter, not idempotent
correction
Strict ~ smoothed Helmholtz operator c(z −
 Refuted as best
 normal-ordered
 Yukawa
 resolvent
 N [A(L +
L)e−tL
 description
 m2 )−1 ], 5× better fit
−1 K is the Green function (off-diagonal); (L + m2 )
"K is the Green function, not K "
 Inverted
is the operator
Repo (QW-2139): "large radial residuals reject
 Both
 right
 and
 right for continuum radial ansätze; wrong operator
the Green claim"
 wrong
 class — graph resolvent passes at nRMSE 0.025
Strict elementwise ≥ 0 , PF structure, eff. ranks
 Confirmed
—
8.21/6.36, 72.7% negative legacy weights
 exactly
4 Part 3 — Can the operator become an adaptive learning operator?
4.1 Well-posedness of the coupled state–operator system
Proposition 4.1 (well-posedness; H).
Let H = C12 , Sym0 ​
 the real symmetric zero-diagonal matrices, and F : H × Sym0 ​
 → Sym0 ​
 locally
Lipschitz (e.g. F = η(Re ψψ
 ∗
 − γK) composed with symmetrization and diagonal deletion). Then the
coupled system
˙ ψ ​
 = iK tψ,
 ​
 ˙ Kt ​
 = F (ψ, K t)
 ​
has a unique local solution on H × Sym0 ​
; K t​
 remains real symmetric zero-diagonal for all t; ∥ψ∥ is(4-1)
conserved up to adiabatic error O(∥ ˙ K ∥/λgap ). ​
 Under timescale separation (∥ ˙ K ∥ ≪ ∥K∥2 ) the
dynamics is unitary to leading order (standard adiabatic estimate).
Proof sketch. Picard–Lindelöf on the product space (the RHS is locally Lipschitz by hypothesis); the
symmetry/zero-diagonal constraints are preserved because F is defined through projection onto Sym0 ;
 ​
d
 ∥ψ∥ 2 = 2Re⟨ψ, iK tψ⟩ ​
 = 0 exactly when Kt ​
 is Hermitian at each t — which holds — so the norm is in
dt
​
fact conserved exactly, not just adiabatically; the adiabatic caveat applies to the energy ⟨ψ∣K t∣ ​
 ψ⟩ and to
eigenstate following, not to the norm. □
So the adaptive program is mathematically legal. The entire question is which F does something interesting. I
tested the candidates from the task numerically (protocols: Appendix A; results in Table 4-1).
11
4 Part 3 — Can the operator become an adaptive learning operator?
4.2 What each learning rule actually does to this operator
Table 4-1. Measured behavior of candidate plasticity rules on FIN kernels (30k iterations unless noted).
New
Rule
 Stability
 Asymptotic object
 Measured (this report)
eigenmodes?
Pure
 Hebb
 ˙ K =
 none;
 ∥K∥
 grows
 ∥K∥ : unbounded linear
unstable
 —
⊤
η ψψ
 linearly
 growth
rank
 collapse
 8.2→1.0
Oja ˙ K = η(ψψ ⊤ −
 projector onto dominant
 no — spans ⊆
stable, bounded
 (stationary drive); corr(K∞
 ​
γK)
 covariance eigendir.[9]
 span(C)
, C )=1.000
GHA/Sanger
 (multi-
 rank-r PCA projector[1
stable
 controlled-rank compression
 no
channel)
 0]
BCM
 (sliding
 drifts off K strict
​
 (corr
stable
 competitive selectivity
 yes (nonlinear)
threshold)[11]
 0.16); rank 4.7
directionality,
 breaks symmetry K = K ⊤ ⇒ breaks unitarity
STDP[12]
 stable
sequences
 of (4-1); needs K = K + + K − split
stable
 by
Predictive coding / free
 precision-matrix
energy ˙ K = −∂K F
 ​
 construction
 learning[30][31][32]
 see §4.5 — structurally the most principled
(Lyapunov F )
Equilibrium
 stable
 (energy-
 quasi-local contrastive
 applies
 verbatim
 with
 energy
 E =
propagation
[29]
 based)
 rule
 − 1 ​
 ψ ⊤ Kψ
24.3 The self-reference theorem: the kernel is never a fixed point of its own Hebbian
dynamics
The decisive experiment for Part 3 is self-driven plasticity: ψ is generated by the kernel's own unitary flow eiKt
(no external signal), and K is updated by Oja's rule. Two facts, one analytic, one measured:
Proposition 4.2 (self-reference obstruction; H 95%).
Let K = ∑m ​
 λ m ​
Pm ​
 be the spectral decomposition and let ψ(t) = eiKtψ0 ​
 with occupation numbers
∣c m ∣ ​
 2
 = ∣⟨v m , ​
 ψ 0 ​
⟩∣2
 . Then the time-averaged excitation covariance is the dephased operator
C = T lim
 →∞​
 1 T ​
 ∫0
 T
 ​
 Re(ψψ ∗ ) dt = ∑ m
 ​
 ∣cm ​
 ∣ 2 Pm ​
,
(cross terms ei(λ m −λn ​
 ) ​
 t dephase unless frequencies are degenerate). Oja's rule converges to (a scaled(4-2)
rank-one truncation of) C , not to K: the fixed-point set of self-driven Hebbian learning is the set of
occupation-weighted mixtures of the kernel's spectral projectors, and K itself is a fixed point only in the
degenerate case C ∝ K (measure-zero occupation profiles). Hence a Hebbian FIN kernel remembers the
statistics of its own excitations, not itself.
Measurement. Stroboscopic protocol (previous analysis's protocol, reimplemented): corr(K30000 , ​
 C) =
1.0000; corr(K30000 ​
, Kstrict ​
) = −0.07 ; effective rank compresses 8.21 → 4.90. Under a continuous-
12
5 Part 4 — The self-modifying kernel on a dynamic graph
time variant the residual correlation with Kstrict ​
 is 0.60 — the numerical residue is protocol-dependent, the
obstruction is not: in no protocol does the dynamics converge back to K . (The previous analysis reported 0.146;
the difference is the drive protocol. Conclusion identical.) Figure 1(d). □
4.4 Falsification of the repository's own "Hebbian emergence" claim
nadsoliton_neural_analysis.py drives the Oja rule with an external signal cos(ωi + t), ω = π/4 , plus
noise, and concludes "Nadsoliton geometry emerges from Hebbian learning" (hardcoded string SIMULATION
VERIFIED ). Replicated exactly (30k iterations, same rates): the learned W correlates 0.9998 with the bare
undamped cosine Toeplitz matrix cos(ω(i − j)) — i.e. with the drive covariance, as Oja's theorem requires —
but only 0.84 with Klegacy​
 under the script's own absolute-distance convention, 0.39 under the cyclic convention,
and 0.43 with K strict ​
 . The damping envelope (1 + βd) −1 and the phase ϕ = π/6 are not recovered (they are
invisible to a resonant drive). The script demonstrates the trivial direction of Oja's theorem — covariance in,
covariance out — and labels the mismatch a verification. The idea (Hebbian self-organization of kernels) is not
thereby refuted; this particular evidence is. [H 99%]
4.5 The principled route: precision learning and the Dyson-like loop
The most principled adaptive extension is not Hebbian at all. In predictive coding the brain learns the precision
matrix Λ = Σ−1 of its generative model by local rules[32]; the precision is the inverse covariance, i.e. the inverse
propagator. Under Theorem 2.2, Kstrict ​
 is already a propagator A(L + m2 ) −1 ; the natural adaptive object is
therefore the operator L t ​
 = A Kt −1 ​
 − m 2 I (off the kernel's null space), learned by gradient flow on a free
energy F = tr(L t ​
C) − log det + (L t ), ​
 whose gradient is the local rule ˙ L = C − L−1 . This is precisely
the self-consistent (Dyson–Schwinger-like) loop in which a medium's propagator is renormalized by the statistics
of its own excitations — the mathematically mature version of what FIN gestures at. Bootstrap structure: the
fixed-point equation
K⋆ ​
 = Φ(C(K⋆ )),
 ​
 C(K) = dephased excitation covariance of eiKt ,
is a well-posed operator fixed-point problem. Iterating K ↦ Φ(C(K)) numerically (random occupations(4-3)
redrawn each step) does not settle — effective rank oscillates 4.1→9.4 — indicating either chaos in the bootstrap
map or strong dependence on the occupation draw; with frozen occupations it collapses to the corresponding
projector mixture (rank 1–3). The map is therefore a nontrivial random dynamical system; its fixed points are
projector mixtures parametrized by occupation vectors, never the seed kernel. [H for the measurements; M for
the conjectured classification of fixed points]
5 Part 4 — The self-modifying kernel on a dynamic graph
5.1 Suitability as the fixed propagator inside an adaptive graph
Take the task's schema literally: information changes the graph, the graph changes future propagation, K stays
mathematically fixed as the propagation law. The candidates' credentials, measured:
Table 5-1. Learning-system desiderata on Z12 ​
 (this report's measurements).
13
5 Part 4 — The self-modifying kernel on a dynamic graph
Property
Spectral locality
Band-pass
filtering
Noise suppression
Feature
separation
Memory
(associative)
Stability
Low-rank
structure
Information
bottleneck
Hierarchical
decomposition
Reservoir
richness
Sparse attractors /
criticality
Edge plasticity
Kstrict
 ​
yes: 79% of ∣λ∣-weighted symbol
below μ = 1
low-pass with notch at μ ≈ 1.18
(PF/noise gain ratio 5.1×)
white-noise power gain 0.54 (−2.7
dB) vs Perron gain 2.76
7/2/3 sign-band tripartition (repo's
own QW-2118 gate)
Hopfield energy − 1​
 ψ⊤ Kψ legal for2patterns
[13][14]
ρ(K) = 1.66 (expanding as amap; stable as Hamiltonian)
eff. rank 8.2/12 — no
3 positive bands at SNR=1 → ~3.2
effective channels
Klegacy
 ​
no: band concentrated at shell, non-
local in space (23% of row mass at
d ≤ 2)
band-pass at shell (by construction)
amplifies shell, suppresses DC
8/0/4
both; capacity ∼ 0.14N ≈ 1.7
ρ = 11.17
 (stiff;
 needs
normalization)
eff. rank 6.4/12; 55% of λ2 mass in
one doublet — near rank-2 in
energy
8.7 effective channels at SNR=1
Assessment
strict: filter-like; legacy:
holographic/mean-field
both, different bands
strict: usable denoiser
strict: 3 usable bands
equal, small at N = 12
strict preferred
legacy: near-low-rank in
energy
strict: natural bottleneck
no multiscale structure on Z12 ​
 (single scale); on ERM: bulk+outliers (§6)
memory capacity MC = 0.10 (of max 12) — useless as-is on Z12
 ​
(symmetry + tiny N )
no nonlinearity ⇒ no attractors; ρ > 1 ⇒ past edge-of-stability as map; as
CTQW generator: criticality N/A
elementwise
 ≥0
 ⇒
 Hebbian
increments preserve sign pattern:
 sign-indefinite: plasticity changes
legal plastic weighted graph at all
 inhibitory/excitatory balance
times
absent as designed
negative result
needs a nonlinearity first
strict: structurally plastic
Verdict of Part 4: Kstrict ​
 has exactly three of the listed virtues genuinely built in — spectral locality with a
notch (matched filtering), elementwise non-negativity (legal plastic adjacency), and a Perron–Frobenius
zero mode (a reference channel). Everything else (memory, attractors, criticality, hierarchy) requires the one
ingredient FIN lacks: a pointwise or threshold nonlinearity. A linear operator, however adaptive its entries,
remains linear in ψ and cannot exhibit attractor dynamics, sparsity, or criticality; this is a theorem, not an opinion.
The smallest change that unlocks the desiderata is analyzed in Part 6. [H]
5.2 The dynamic-graph loop
With a learning rule F as in §4, the schema "information changes the graph, the graph changes propagation" is
the map (ψt , ​
 G t) ​
 ↦ (ψt+1 , ​
 Gt+1) ​
 , K fixed as propagation law. This is precisely graph structure learning
(joint learning of adjacency and dynamics)[43] with the twist that the propagation kernel is frozen and only the
adjacency deforms. Stability requires homeostasis (decay/normalization), otherwise §4.2's pure-Hebb runaway
14
6 Part 5 — Deep spectral analysis
applies to the graph. There is one structural subtlety worth a theorem's emphasis: because K is circulant,
translation invariance is broken by any nontrivial learning — the Fourier diagonalization survives only in
expectation under stationarily translation-invariant drives. The clean spectral theory of §2 is a property of the
unlearned object; after learning, the object is a generic symmetric matrix and must be analyzed by generic tools.
[H]
6 Part 5 — Deep spectral analysis
6.1 On the cycle: everything is exact
Eigenvectors. Both kernels are diagonalized by the DFT basis: eigenmodes are plane waves vm ​
 ( j) =
e
2πimj/12
 / 12​
 (or standing waves cos, sin inside each degenerate doublet). Consequently every mode is
maximally delocalized: inverse participation ratio IPR(vm ) ​
 = ∑j ​
 ∣vm ​
 ( j)∣ 4 = 1/12 for all m (analytic,
exact). There is no localization physics on the cycle whatsoever — all "modes" are global Fourier channels. [H,
exact]
Band decomposition. With the repository's own threshold (ε = 0.1 , QW-2118): strict splits 7/2/3
(negative/near-zero/positive), reproducing the artifact's tripartition exactly; legacy splits 8/0/4. The sign of λm
 ​
is the "phase" of the channel: strict is a low-pass amplifier with a notch (zero crossing between μ2 ​
 = 1.00 and
μ 3 ​
 = 2.00 ; fitted notch μ ⋆ ​
 = 1.18), legacy a shell amplifier with a negative DC floor.
Spectral entropy and effective rank. From the energy distribution w m ​
 = λ2 m / ​
 ∑ λ 2 : H(K strict ) ​
 = 2.65
bits (max log 2 ​
 12 = 3.58 ), H(K legacy ​
) = 2.23 bits; participation ranks exp H = 6.30 and 4.69; ∣λ∣-
ranks 8.21 and 6.36 . Legacy's energy is dominated by the m = 2 doublet (55.4% of λ2
 -mass; 43.6% of ∣λ∣-
mass — both conventions reported, since the previous analysis used the latter). Legacy is thus "effectively rank-2
plus a DC pedestal" in energy; strict is spectrally diverse. Information channels. Treating K as a MIMO channel
y = Kx + n with singular values ∣λm ∣: ​
 strict supports 3.2 effective channels at SNR 1 and 7.8 at SNR 10;
legacy 8.7 and 10.6 . Strict is the natural bottleneck; legacy the broadband medium. [H]
Latent representations. Do dominant modes resemble "learned" representations? On the cycle the question is
vacuous: the modes are Fourier, fixed by symmetry, learnable by any PCA-like rule in one step (§4.2) precisely
because they carry no data structure. On the ERM (below) the question is real, and the answer is: the top
eigenvectors of a random-cloud kernel matrix are a diffusion-map-like embedding of the cloud[22][23] — but that
is a property of all smooth distance kernels, not of FIN's parameters. [H]
6.2 On random Euclidean clouds: the kernel becomes a random matrix
Following the repository's QW-593 protocol (Gaussian clouds in R3 , scale 2.0, N up to 96), the kernel matrices
are Euclidean random matrices M ij ​
 = f(∥xi ​
 − xj ​
∥) in the sense of Mézard–Parisi–Zee[2][3]. Measured
over 60–120 clouds:
Table 6-1. ERM spectral statistics (this report).
Quantity
 legacy-ERM
 strict-ERM
 Reading
15
6 Part 5 — Deep spectral analysis
10.5 → 15.0 → 19.4 for 10.8 → 18.5 →
 sublinear growth: strong spectral
Effective rank vs N
 N = 24 → 96
 (ratio 29.1 (ratio 0.45 →
 concentration persists in the large-
0.44 → 0.20)
 0.30 )
 N limit
near Wigner–Dyson (0.536), far
Level
 statistics
 (⟨r⟩
0.509
 0.501
 from Poisson (0.386): level
ratio)
repulsion[50]
bulk modes localized on ~4–5 sites
IPR×N: bulk median
 4.5
 —
(ERM-typical[3])
top
 modes
 =
 smoothest,
IPR×N: top 5% modes
 2.3 (most extended)
 —
embedding-like
IPR×N:
 bottom
 5%
 up to 45 (strongly localized, ~2
 deep-negative modes = close-pair
—
modes
 sites)
 resonances
Spectral
 dimension
diverges (mean-field: weights
 ≈ 3 at intermediate t (Fi strict recovers the embedding
d s (t) ​
 of ∣K∣-graph
decay too slowly)
 g. 2(c))
 dimension; legacy does not
Laplacian
Three observations are referee-relevant. (i) The DOS has the canonical ERM shape: a Marchenko–Pastur-like bulk
plus isolated outliers at the band edges (Figure 2(a)) — the "spectral tripartition" the repository repeatedly
discovers in its gates is the generic ERM bulk/outlier anatomy, not a fingerprint of its parameters. (ii) The
localization gradient (extended top modes, strongly localized bottom modes, Figure 2(b)) is the known ERM
phonon-fraction picture[3]. (iii) The spectral-dimension contrast is the sharpest discriminator between the two
kernels in the continuum setting: strict's strong damping makes the weighted graph see the underlying 3D
geometry; legacy's slow decay makes it mean-field. If anything in FIN deserves the adjective "geometric," it is the
strict damping exponent η = 1.8 — the one parameter the project treats as a mere fit value. [H for
measurements; M for the ds ​
 estimate, which is size-limited]
16
7 Part 6 — The mathematically smallest extension
6.3 Graph Fourier basis and the "information manifold" question
On the cycle, the graph Fourier basis is the DFT basis and the kernel is a multiplier λm ​
 — a textbook GSP
filter[7]
. The "low-dimensional manifold" of states reachable under (1-4) from generic ψ0 ​
 is the support of the
occupation vector: with all 12 modes occupied, the orbit {eiKt
 ψ0 ​
} densely covers a torus in CP11
 whose
dimension equals the number of rationally independent λ m ​
 — generically 12, i.e. no dimension reduction
occurs in the state dynamics; reduction happens only in the covariance (Eq. (4-2)), which is what Hebbian rules
see. This is the mathematical reason Part 3's obstruction (Prop 4.2) is structural rather than accidental. [H]
7 Part 6 — The mathematically smallest extension
7.1 Parameter sensitivity: what one knob can and cannot do
The spectral Fisher information matrix Fab ​
 = ∑m ​
(∂λm ​
/∂θa ​
)(∂λm ​
 /∂θ b ), ​
 θ = (ω, ϕ, β, η), computed
by exact differentiation of the circulant symbol (Figure 2(d)):
6.58
 2.48
 2.64
 1.67
2.48
 1.32
 1.62
 0.59
F =
 ​
 ​
 ​
 ​
 ​
 ​
 ,
2.64
 1.62
 2.07
 0.63
1.67
 0.59
 0.63
 0.57
eig(F ) = {9.33, 1.08, 0.13, 0.00}.
Proposition 7.1 (non-identifiability; H 96%).
 (7-1)
(i) F has a numerically exact null direction: the 4-parameter family (ω, ϕ, β, η) ↦ K has effective
dimension ≤ 3 on the domain d ∈ {1, … , 6} ; the parameters are sloppy — individually unrecoverable
from
 the
 kernel's
 action.
 Exhibit:
 the
 tuples
 (0.18575, 0.1625, 1.0, 1.8) and
(0.28237, −0.57812, 1.03534, 2.02245) generate profiles agreeing to maxd ​
 ∣ΔK(d)∣ = 6.1 ×
10 −5 .
(ii) Sensitivity ordering: per unit absolute change, ω dominates (log 10 ​
 ∥∂ω λ∥ ​
 = 0.41 , stiffest
eigendirection 82% ω ); per 1% relative change, β (0.50%) > η (0.39%) > ω (0.10%) > ϕ (0.06%). ϕ is the
weakest knob by both measures.
Consequence for any adaptive program: "adaptive ϕ" is nearly free (it changes almost nothing), "adaptive ω " is
the strong control (it moves the whole spectrum), and β, η are the mid-range structural knobs controlling locality.
It also means the repository's fixation on certifying exact parameter values is mathematically idle: a 1-parameter
curve of tuples realizes the same operator.
7.2 Ranking the minimal extensions
Table 7-1. Minimal extensions ranked by mathematical disruption vs gain.
Extension
 Disruption
 Unlocks
Verdict
17
8 Part 7 — Literature map: has FIN rediscovered a known operator?
Dynamic
 normalization
 breaks
 symmetry
 mildly;
 legal
 Markov
 propagation;
smallest
 useful
 step;
Kt ​
 = Dt −1 ​
Wt
 ​
 (row-
 preserves non-negativity, PF
 diffusion maps
[22]
; stability (
recommended
stochastic)
 mode, filter structure
 ρ = 1 )
Adaptive ωt ​
 driven by
1 scalar ODE;
 circulant
 self-tuning
 filter
 (resonance
occupation of the notch
 smallest spectral control
preserved
 tracking)
region
Adaptive βt , η t ​
 (locality
 dynamic range control of the
2 scalars; circulant preserved
 useful with ω
homeostasis)
 graph
Adaptive ϕt
 ​
 minimal
 almost nothing (Prop 7.1)
 skip
Hebbian
 edges
 ˙ K =
 breaks circulant structure;
 PCA of excitation statistics
 standard; fixed points are
η(ψψ ⊤ − γK)
 needs decay
 (§4)
 projector mixtures
system identification; well-
Bayesian
 update
 of
none structural (offline)
 posed given §7.1 caveats (3 of
 clean, doable
(ω, β, η) from observed ψ
4 params)
Entropy maximization over
 recovers
 dynamic
none
 equivalent to row 1
row weights
 normalization
Energy minimization (free-
 reinterprets
 K
 as
Lyapunov-stable
 learning;
energy/precision
 learning,
 propagator,
 learns
 L =
 most principled
Dyson-like loop
§4.5)
 AK −1 − m 2
Pointwise nonlinearity
 σ
 attractors, sparsity, criticality,
 the
 one
 extension
 that
after
 K
 (e.g.
 ψ ↦ leaves linear theory
 memory — everything in §5
 actually
 matters;
σ(Kψ) )
 that linearity forbids
 mathematically the largest
Answer to Part 6. The mathematically smallest extension that preserves all existing mathematics (circulant, filter
interpretation, Thm 2.2) is dynamic row-normalization plus scalar adaptation of (ω, β, η) — it converts the
object into a diffusion-map-with-tunable-filter, fully inside known theory. The extension that preserves the most
mathematics while unlocking learning in the strong sense is the precision/free-energy route (learn the operator L
whose propagator is K ; Lyapunov-stable by construction). The extension that unlocks the learning-system
desiderata of Part 4 is a nonlinearity — and it is precisely the one that destroys the kernel's current mathematical
identity. This trilemma is the honest answer. [H]
8 Part 7 — Literature map: has FIN rediscovered a known operator?
8.1 Systematic sweep
I searched the literatures named in the task (Green functions in neural computation; operator learning; graph
kernels; adaptive graph operators; graph signal processing; spectral neural operators; Koopman learning;
continuous attractors; Hopfield networks; reservoir computing; Hebbian graph learning; kernel adaptive filters;
graph diffusion learning; wave operators; random geometric operators; Euclidean random matrices). The mapping,
with verdicts:
Table 8-1. Literature map and rediscovery verdicts.
FIN object
 Closest known object
 Distance
 Verdict
18
8 Part 7 — Literature map: has FIN rediscovered a known operator?
Chung–Yau discrete Green
nRMSE 0.025
 +
 normal
 REDISCOVERY (with a harmless
Kstrict ​
 on Z12
 ​
 function,
 screened
 variant
 ordering
 diagonal twist)
(L + m2 )−1 [1]
spectral-shell band-pass graph
 REDISCOVERY (of a concept,
Klegacy ​
 on Z12
 ​
 filter[7][8]
 corr 0.95 with shell mask
 not a named formula)
REDISCOVERY; repo's "spectral
Kernel on random
 Euclidean random matrices[2]
exact class membership
 tripartitions" = ERM bulk/outlier
clouds
 [3][4]
anatomy
one
 Graph-Kernel-Network
 GNO learns κ ϕ ; ​
 FIN freezes
 known architecture,
 untrained
ψ ↦ Kψ as layer
layer, frozen kernel
[5][6]
 it
 instance
continuous-time
 quantum
˙ ψ​
 = iKψ
 walk[37]; linear reservoir[20][2
 exact
 known
1]
GSP graph filters; ChebNet
Kernel as filter λ(μ)
 exact concept
 known
with fixed coefficients[7][41]
Oja/GHA[9][10];
 kernel
Hebbian Kt
 ​
 adaptive filtering[28]; graph
 exact rules
 known
structure learning
[43]
Kernel
 as
 energy
 Hopfield (linear regime)[13];
missing nonlinearity/readout
 known, truncated
− 1 ​
 ψ ⊤ Kψ
 dense associative memory[15]
2Attention-like
 frozen
 attention[17][18][19];
no data-dependence
 known, static instance
reading
 Hopfield–attention identity[16]
Eigenmodes
 as
 diffusion maps / Laplacian
exact (after PSD-ification)
 known
embeddings
 eigenmaps[22][23]
Kernel on graphs as
 Smola–Kondor graph kernels;
PSD required ⇒ excluded
 adjacent, not a member
regularizer
 Vishwanathan et al.[24][25]
Indefinite
 kernel
Kreĭn-space machines[26][27]
 exact habitat
 known habitat, unexplored by repo
learning
Learning
 Green
 neural Green's functions[44][4
 inverse direction (they learn;
known program, inverse instance
functions from data
 5][46]
; DeepGreen
 FIN posits)
Point-cloud
 kernel
Nyström method[47][51]
 exact
 known (1930)
matrices
Hamiltonian
 Vanchurin[33];
 equilibrium
[29]
 closest: predictive coding
 known
 marriages,
 different
dynamics + learning
 propagation
 ;
 predictive
with learned precision (§4.5)
 ceremonies
marriage
 coding[30][31][32]
Ring
 +
 localized
 continuous
 attractor
 missing
 stabilizing
skeleton of a known object
coupling + bump
 networks
[48][49]
 nonlinearity
19
9 Part 8 — What is genuinely novel (mathematics only)
8.2 Verdict of Part 7
FIN has unintentionally rediscovered known operators — twice.Kstrict_gate ​
 is, to 2.9% Frobenius
residual, a screened Chung–Yau discrete Green function on the cycle with the diagonal removed (Theorem
2.2)[1]; K legacy_ont ​
 is a spectral-shell band-pass graph filter[7]; and the random-cloud experiments are
Euclidean random matrix theory exactly as in Mézard–Parisi–Zee and Goetschy–Skipetrov[2][3]. Nothing in
the kernel's statics is new. The only place the search returns "not found" is the dynamical combination: a
normal-ordered oscillatory Green kernel used simultaneously as a unitary generator and as the slow variable
of a self-referential plasticity loop (Eq. (4-3)). That combination appears to be unstudied as a named object;
its ingredients are all standard. [H for the identifications; M for the negative "unstudied" claim —
absence of evidence is not evidence of absence]
9 Part 8 — What is genuinely novel (mathematics only)
After §8, the genuinely new content is small and can be stated precisely:
1. The normal-ordering twist. Deleting the diagonal of a discrete Green function and using the remainder as
a standalone operator is a legitimate, mildly nonstandard move: it converts a resolvent into a zero-diagonal
Hamiltonian/filter with a built-in negative spectral floor (the contact term), at the price of indefiniteness. I
am not aware of this "normal-ordered resolvent" being studied as an object. It is one paragraph of
mathematics, not a theory. [M]
2. The self-referential bootstrap problem K ⋆ ​
 = Φ(C(K⋆ )) ​
 (Eq. (4-3)): a Hamiltonian whose slow
Hebbian plasticity is driven by the dephased statistics of its own unitary excitations.
Existence/classification of fixed points (they are occupation-weighted projector mixtures, Prop 4.2), their
stability, and the random-dynamical-system behavior under redrawn occupations (§4.5) are, to my
knowledge, unstudied. This is the one place where FIN's actual practice (QW-540's ˙ K = η∣ψi ψj ​
 ∣ ​
 −
decay K ; the hebbian_* audit family) touches an open mathematical question without knowing it. [M]
3. Nothing else. The αgeo ​
 = 4 ln 2 identification is numerology (it enters the operator as a scalar and drops
out of all normalized spectral questions); the "Shannon void asymmetry" is a narrative; the tripartition
gates are ERM genericity. [H]
10 Part 9 — A publishable contribution without physics
Yes, one. Stripped of particles and cosmology, the following is a complete, honest, referee-survivable paper:
20
11 Part 10 — Blind classification: top-10 fields, ranked
Working title:Normal-ordered discrete Green functions as adaptive graph operators: exact identification,
spectral anatomy, and the self-referential Hebbian bootstrap.
Contributions: (1) Thm 2.2-style identification of oscillatory damped kernels with normal-ordered screened
resolvents, with the Möbius equivalence of cutoff-operator and resolvent-plus-contact families (a small but
clean lemma people do get wrong — §3.3 shows two independent audits got it wrong in opposite directions);
(2) the full spectral anatomy (symbols, bands, ERM statistics, spectral dimension) of one parametric family
across its cycle and random-cloud instances; (3) Prop 4.2's self-reference obstruction and the bootstrap fixed-
point problem with first measurements; (4) the sloppy-parameter/non-identifiability analysis (Prop 7.1).
What referees will demand: proofs of (1) as asymptotic statements in N (currently exact at N = 12 —
must be redone as N → ∞ circulant asymptotics, which is standard); one non-toy learning task
demonstrating the adaptive operator doing something (classification-by-diffusion suffices); removal of all
ToE residue.
Realistic targets (ranked):Linear Algebra and its Applications (best fit: exact operator identification +
spectral theory); Applied and Computational Harmonic Analysis (if framed as GSP/filter design); Journal of
Complex Networks or Physical Review E (ERM/random-geometric framing, statistical framing); Neural
Computation (if the learning-loop half leads).
Probability of acceptance at the best-fit venue, honestly: 40–60%, contingent on the N → ∞ upgrade.
[M 85% that the contribution as scoped is real; the acceptance estimate is judgment, not data]
A second, smaller note exists: the counterexample pair "two audits of the same kernel reached opposite operator-
vs-Green assignments, both wrong" — as a methods comment on fitting operator families to spectral symbols.
Not worth a separate paper; fold into the above as a cautionary section.
11 Part 10 — Blind classification: top-10 fields, ranked
Protocol: I ignore all repository text and answer only from the object's properties as established above: a
parametric radial function K(d) = cos(ωd + ϕ)/(1 + βdη ) evaluated (i) as a circulant on Z12 ​
 with zero
diagonal, used as a unitary generator, and (ii) as a Euclidean random matrix on point clouds; symbol a low-pass-
with-notch multiplier; elementwise non-negative; indefinite; normal-ordered-resolvent structure at 0.029 residual.
Ranking by posterior probability that a specialist handing me this object would name that field as its home, with
quantitative justification from the measured discriminants:
Table 11-1. Blind classification ranking.
#
 Field
 Prob.
 Quantitative reasoning
Spectral graph theory /
 Exact circulant; resolvent structure (nRMSE 0.025) = Chung–Yau object[1];
1
 discrete
 harmonic
 0.28
 Laplacian eigenbasis is the object's own diagonalization. No other field explains
analysis
 the 0.9993 Yukawa correlation.
The symbol λ(μ) is literally a graph-filter frequency response[7]; notch, bands,
2
 Graph signal processing
 0.18
 tripartition are GSP primitives; object is used as a multiplier on the graph-Fourier
basis.
21
12 Falsification log, counterarguments, and limitations
3
4
5
6
7
8
9
10
On
 clouds:
 exact
 ERM
 class[2][3];
 Wigner–Dyson
 ⟨r⟩ = 0.51 ;
Random matrix theory
 /
0.14
 bulk+outlier DOS;
 localization
 gradient
 —
 all
 textbook
 ERM
Euclidean random matrices
phenomenology.
Radial
 translation-invariant
 kernel;
 resolvent/contact
 decomposition;
Integral operators / Green's-
fractional-kernel envelope (d−1.8 , Riesz s = 1.2 reading)[38]. Penalized
function
 methods
 (applied
 0.10
because the strict identification is discrete-native (continuum FT fails, corr
analysis)
0.07).
Quantum
 walks
 /
 quantum
 ψ ˙ ​
 = iKψ is CTQW verbatim[37]; PF mode as ground state; dephased
0.08
information on graphs
 covariance is the standard CTQW time-average.
Machine
 learning:
 graph
 One frozen GKN layer[5]; indefinite kernel (Kreĭn habitat)[26]; but no
0.07
kernels & neural operators
 training, no RKHS, no nonlinearity.
Hopfield energy legal; Oja dynamics natural; ring-with-local-coupling =
Computational
 neuroscience
0.06
 attractor skeleton[48]; but MC=0.10 and linearity bar most of the field's
(Hebbian/attractor networks)
content.
Reservoir
 /
 echo-state
 Fixed linear recurrent operator[20]; but ρ = 1.66 unnormalized, N = 12,
0.04
computing
 no readout, MC 0.10.
Koopman/transfer-operator
 eiKt is a unitary Koopman semigroup on its own state space[35][36]; generic
0.03
analysis
 to all linear flows, hence low specificity.
Stochastic processes on graphs /
 One row-normalization away (§7.2)[22]; but as given, row sums 1.66 ≠ 1 —
0.02
diffusion geometry
 not a transition operator.
Residual probability mass 0.00–0.05 distributed over: fractional calculus (Riesz reading), information geometry
(no Fisher metric anywhere in the object — the 4 ln 2 is a scalar prefactor, not a metric), signal processing on
directed graphs (complex variant). The modal answer: spectral graph theory, with GSP the engineering twin.
[H on the discriminants; the probabilities themselves are calibrated judgment]
12 Falsification log, counterarguments, and limitations
12.1 What I tried to break, and what survived
Table 12-1. Falsification log.
Hypothesis attacked
 Attack
 Outcome
Fitted 9 rival families incl. both previous claims; checked
Strict
 =
 Yukawa
 survived
 all;
 continuum
 reading
continuum reading; checked uniqueness over parameter
resolvent (Thm 2.2)
 correctly excluded
box
Legacy
 =
 on-shell
 Tested off-shell families;
 idempotency;
 continuum
 survived,
 with
 projector→filter
object
 decomposition Eq. (3-1)
 correction
Previous
 analysis's
 Reproduced it exactly, then extended search by one
 falsified as "best" (5× worse); its
strict fit
 rational family
 operator-vs-Green inversion falsified
Repo's
 Hebbian
Exact replication, both distance conventions
 falsified (drive-covariance artifact)
emergence
22
Appendix A — Reproducibility of all numerical claims
My own Prop 4.2 protocol-
 Two independent protocols (stroboscopic,
 residual
 corr
 varies
 (−0.07…0.60);
dependence
 continuous-time)
 obstruction invariant
Both kernels; 120 clouds; bulk-trimmed
 ⟨r⟩ 0.51/0.50 — survives; note N = 96
ERM Wigner–Dyson claim
spacings
 is small
log-derivative of return probability, 30 time
 plateau visible but size-limited; downgrade
Spectral dimension ≈ 3
points
 to M
12.2 Honest limitations of this report
1. Finite-size cage. The strict identification is exact-grade at N = 12 ; the N → ∞ asymptotic statement
(circulant Yukawa limits) is standard but not proven here. Everything ERM-side is N ≤ 96 .
2. Negative literature claims. "Unstudied combination" (§9) is a search verdict, not a theorem; the self-
referential bootstrap may live in a literature I did not reach (e.g. self-consistent reservoir adaptation).
3. The bootstrap map's fixed-point classification (§4.5) is empirical; a proof of "fixed points = projector
mixtures" in full generality is missing (Prop 4.2 proves it only for the dephased-covariance limit).
4. Convention sensitivity. The repository mixes distance conventions (Table 1-1); I analyzed the gate-
certified one. Under the absolute-distance convention, strict-kernel statements are numerically similar
(Toeplitz vs circulant at N = 12 differs at the boundary) but not identical; I did not redo the full battery
for the path graph.
5. Confidence calibration. H labels attach to exact or nine-digit-reproducible facts; probability assignments
in §10–11 are structured judgment and should be read as such.
12.3 Final answer to the task's framing question
Does the mathematical core of FIN have a more natural interpretation as an adaptive information-propagation
operator than as a physical field theory? Yes — with a correction of terms. The natural interpretation of the
strict kernel is not yet "adaptive": it is a propagation medium (a normal-ordered discrete Green function), and
"adaptive" is a well-posed upgrade with exactly one obstruction (Prop 4.2) and one trilemma (§7.2). The legacy
kernel is not a propagation operator at all but a shell selector. The field-theoretic reading adds nothing to either
object — every physical claim in the repository reduces to properties these operators already have, or fail to have,
as graph-theoretic objects. The information-dynamics reading is therefore not merely more charitable; it is the
only one that makes the object's actual mathematical content visible. Whether that content is worth developing is
answered in §10: one honest paper, modest venue, no paradigm. [H]
Appendix A — Reproducibility of all numerical claims
Kernels. Eq. (1-1): α = 4 ln 2 , ω = π/4 , ϕ = π/6 , β = 0.01 ; Eq. (1-2): ω = 0.18575, ϕ = 0.1625
, β = 1, η = 1.8 . Domain Z12, ​
 cyclic min-distance, zero diagonal (artifact report_qw2118…json flags:
matrix_symmetric=True , eigenvalues reproduced to 9 digits).
Symbol formula. Eq. (2-1). Fits: least squares over m = 0..6 (curve_fit, maxfev 5 × 10 4 ). nRMSE as defined
in §2.2. Thm 2.2 residual: ∥R∥ F ​
/ ∥K∥ F ​
 = 0.0286; AG 00 ​
 = 1.08836 vs ∣C∣ = 1.08960 .
Continuum transforms. 3D Fourier–Hankel by adaptive quadrature to r = 2000 (legacy: shell-peaked FT,
sign changes; strict: Yukawa fit corr 0.070, power-law ∣FT ∣ ∼ q −0.63 mid-range).
23
References
Learning experiments. Repo replication: Oja, lr=decay=0.01, 30k iters, drive cos(πi/4 + t) + N (0, 0.5) ;
results §4.4. Self-driven stroboscopic: lr=decay=0.02, 30k iters, flight time 0.7, random phase per step; corr
(K, C) = 1.0000, corr(K, Kstrict ​
) = −0.0697 ; rank 8.21→4.90. Continuous-time: eiKΔt , Δt = 0.4,
lr 0.02, 4k iters; residual corr 0.60; rank 8.21→3.57. Pure Hebb: linear ∥K∥ growth. BCM: lr 0.05, 8k iters; corr
0.156, rank 4.68. Bootstrap: K ↦ Φ(C(K)) , 60 steps, occupation redrawn; rank oscillates 4.1–9.4.
ERM. Gaussian clouds N (0, 22I3 ​
); N ∈ {24, 48, 96} ; 60–120 realizations. Spacings: unfolded per
realization, 6 modes trimmed at each edge, pooled; ⟨r⟩ as in Table 6-1. Spectral dimension: ds ​
 ( t) =
−2 ∂ln t​
 ln P (t) , P (t) = mean return probability of ∣K∣ -weighted combinatorial Laplacian heat kernel, 30
log-spaced t ∈ [10 −2 , 102 ] .
Fisher. Eq. (7-1) by forward differences (ε = 10−4 , verified stable at 10 −5 ); null direction residual ∣λmin ∣ ​
 <
10 −6 . Alternate tuple residual 6.1 × 10−5 .
Filter gains. White-noise power gain ∑ λ2 m /12 ​
 = 0.538 (strict); Perron gain 1.6603 2 = 2.757.
Reservoir MC: tanh nonlinearity, input gain 0.1, 11 delays, 4000 steps: MC = 0.10.
All random seeds fixed (7, 11, 42, 1–5 depending on experiment); NumPy/SciPy implementations; no repository
code was executed for any reported number.
References
1. Chung, F., & Yau, S.-T. (2000). Discrete Green's functions. Journal of Combinatorial Theory, Series A, 91(1–2), 191–
214.
2. Mézard, M., Parisi, G., & Zee, A. (1999). Spectra of Euclidean random matrices. Nuclear Physics B, 559(3), 689–701.
3. Goetschy, A., & Skipetrov, S. E. (2013). Euclidean random matrices and their applications in physics.
arXiv:1303.2880.
4. Bogomolny, E., Bohigas, O., & Schmit, C. (2003). Spectral statistics of Euclidean random matrices. Journal of
Physics A, 36(12), 3595–3616.
5. Li, Z., Kovachki, N., Azizzadenesheli, K., Liu, B., Bhattacharya, K., Stuart, A., & Anandkumar, A. (2020). Neural
operator: Graph kernel network for partial differential equations. arXiv:2003.03485.
6. Li, Z., Kovachki, N., Azizzadenesheli, K., Liu, B., Bhattacharya, K., Stuart, A., & Anandkumar, A. (2021). Fourier
neural operator for parametric partial differential equations. ICLR 2021.
7. Shuman, D. I., Narang, S. K., Frossard, P., Ortega, A., & Vandergheynst, P. (2013). The emerging field of signal
processing on graphs. IEEE Signal Processing Magazine, 30(3), 83–98.
8. Hammond, D. K., Vandergheynst, P., & Gribonval, R. (2011). Wavelets on graphs via spectral graph theory. Applied
and Computational Harmonic Analysis, 30(2), 129–150.
9. Oja, E. (1982). A simplified neuron model as a principal component analyzer. Journal of Mathematical Biology, 15,
267–273.
10. Sanger, T. D. (1989). Optimal unsupervised learning in a single-layer linear feedforward neural network. Neural
Networks, 2(6), 459–473.
11. Bienenstock, E. L., Cooper, L. N., & Munro, P. W. (1982). Theory for the development of neuron selectivity. Journal
of Neuroscience, 2(1), 32–48.
12. Bi, G.-Q., & Poo, M.-M. (1998). Synaptic modifications in cultured hippocampal neurons. Journal of Neuroscience,
18(24), 10464–10472.
13. Hopfield, J. J. (1982). Neural networks and physical systems with emergent collective computational abilities. PNAS,
79(8), 2554–2558.
24
References
14. Amit, D. J., Gutfreund, H., & Sompolinsky, H. (1985). Storing infinite numbers of patterns in a spin-glass model of
neural networks. Physical Review Letters, 55(14), 1530–1533.
15. Krotov, D., & Hopfield, J. J. (2016). Dense associative memory for pattern recognition. NeurIPS 29.
16. Ramsauer, H., Schäfl, B., Lehner, J., et al. (2021). Hopfield networks is all you need. ICLR 2021.
17. Vaswani, A., Shazeer, N., Parmar, N., et al. (2017). Attention is all you need. NeurIPS 30.
18. Katharopoulos, A., Vyas, A., Pappas, N., & Fleuret, F. (2020). Transformers are RNNs: Fast autoregressive
transformers with linear attention. ICML 2020, 5156–5165.
19. Tsai, Y.-H. H., Bai, S., Yamada, M., Morency, L.-P., & Salakhutdinov, R. (2019). Transformer dissection: A unified
understanding of transformer's attention via the lens of kernel. EMNLP 2019.
20. Jaeger, H. (2001). The "echo state" approach to analysing and training recurrent neural networks. GMD Technical
Report 148.
21. Maass, W., Natschläger, T., & Markram, H. (2002). Real-time computing without stable states. Neural Computation,
14(11), 2531–2560.
22. Coifman, R. R., & Lafon, S. (2006). Diffusion maps. Applied and Computational Harmonic Analysis, 21(1), 5–30.
23. Belkin, M., & Niyogi, P. (2003). Laplacian eigenmaps for dimensionality reduction and data representation. Neural
Computation, 15(6), 1373–1396.
24. Vishwanathan, S. V. N., Schraudolph, N. N., Kondor, R., & Borgwardt, K. M. (2010). Graph kernels. JMLR, 11, 1201–
1242.
25. Smola, A. J., & Kondor, R. (2003). Kernels and regularization on graphs. COLT 2003, 144–158.
26. Ong, C. S., Mary, X., Canu, S., & Smola, A. J. (2004). Learning with non-positive kernels. ICML 2004.
27. Oglic, D., & Gärtner, T. (2019). Scalable learning in reproducing kernel Kreĭn spaces. ICML 2019 (PMLR 97).
28. Liu, W., Principe, J. C., & Haykin, S. (2010). Kernel Adaptive Filtering: A Comprehensive Introduction. Wiley.
29. Scellier, B., & Bengio, Y. (2017). Equilibrium propagation. Frontiers in Computational Neuroscience, 11, 24.
30. Friston, K. (2010). The free-energy principle: a unified brain theory? Nature Reviews Neuroscience, 11, 127–138.
31. Rao, R. P. N., & Ballard, D. H. (1999). Predictive coding in the visual cortex. Nature Neuroscience, 2, 79–87.
32. Bogacz, R. (2017). A tutorial on the free-energy framework for modelling perception and learning. Journal of
Mathematical Psychology, 76, 198–211.
33. Vanchurin, V. (2021). The world as a neural network. Entropy, 23(9), 1210.
34. Jacot, A., Gabriel, F., & Hongler, C. (2018). Neural tangent kernel. NeurIPS 31.
35. Koopman, B. O. (1931). Hamiltonian systems and transformation in Hilbert space. PNAS, 17(5), 315–318.
36. Mezić, I. (2005). Spectral properties of dynamical systems, model reduction and decompositions. Nonlinear
Dynamics, 41, 309–325.
37. Farhi, E., & Gutmann, S. (1998). Quantum computation and decision trees. Physical Review A, 58(2), 915–928.
38. Kwaśnicki, M. (2017). Ten equivalent definitions of the fractional Laplace operator. Fractional Calculus and Applied
Analysis, 20(1), 7–51.
39. Rasmussen, C. E., & Williams, C. K. I. (2006). Gaussian Processes for Machine Learning. MIT Press.
40. Rahimi, A., & Recht, B. (2007). Random features for large-scale kernel machines. NeurIPS 20.
41. Defferrard, M., Bresson, X., & Vandergheynst, P. (2016). Convolutional neural networks on graphs with fast localized
spectral filtering. NeurIPS 29.
42. Kipf, T. N., & Welling, M. (2017). Semi-supervised classification with graph convolutional networks. ICLR 2017.
43. Franceschi, L., Niepert, M., Pontil, M., & He, X. (2019). Learning discrete structures for graph neural networks. ICML
2019; and Zhu, Y., et al. (2021). Deep graph structure learning for robust representations: A survey. arXiv
preprint.
25
References
44. Boullé, N., Earls, C. J., & Townsend, A. (2022). Data-driven discovery of Green's functions with human-
understandable deep learning. Scientific Reports, 12, 4824.
45. Teng, Y., Zhang, X., Wang, Z., & Ju, L. (2022). Learning Green's functions of linear reaction-diffusion equations.
Proc. Mathematical and Scientific Machine Learning.
46. Negi, P., Cheng, M., Krishnamurthy, K., Ying, W., & Li, S. (2024). Learning domain-independent Green's function for
elliptic partial differential equations. Computer Methods in Applied Mechanics and Engineering.
47. Williams, C. K. I., & Seeger, M. (2001). Using the Nyström method to speed up kernel machines. NeurIPS 13.
48. Amari, S. (1977). Dynamics of pattern formation in lateral-inhibition type neural fields. Biological Cybernetics, 27(2),
77–87.
49. Wu, S., Hamaguchi, K., & Amari, S. (2008). Dynamics and computation of continuous attractors. Neural
Computation, 20(4), 994–1025.
50. Bohigas, O., Giannoni, M.-J., & Schmit, C. (1984). Characterization of chaotic quantum spectra. Physical Review
Letters, 52(1), 1–4.
51. Nyström, E. J. (1930). Über die praktische Auflösung von Integralgleichungen. Acta Mathematica, 54, 185–204.
Primary sources (repository artifacts, cited inline): hyconiek/Fractal-Nadsoliton-Theory (GitHub, state 2026-07-17):
RELEASE_6_2_SHANNON_VOID_ASYMMETRY_EN_PL.md ;
TOE_FINAL_DOCUMENTATION_RELEASE_4_4_LEGACY_FULL.tex ;
report_qw2118_ktotal_spectral_tripartition_gate.json ;
RAPORT_QW2117_KTOTAL_LOCALITY_OPERATOR_AUDIT.md ;
RAPORT_QW2139_KERNEL_GREEN_STATUS_3D_ENERGY_GATE.md ;
RAPORT_QW2134_INTERACTING_MICROCAUSALITY_PERTURBATIVE_GATE.md ; nadsoliton_neural_analysis.py ;
QW-593_Information_Unity.py ;
 AGENTS.md ;
 Recenzja_adwersarialna_i_roadmapa_FIN.md ;
Analiza_operatora_jadra_FIN.md (the previous analysis under test).
26
