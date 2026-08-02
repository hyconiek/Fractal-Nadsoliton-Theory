# FIN Local Research — Release 10.46

## Checkpoint P475–P480–P484

### Dependency-Minimal Exact Replay, Algebraic Elimination Burden, and an Adversarial Phase-Face Audit

**Creator:** Żuchowski, Krzysztof  
**Affiliation:** Independent Researcher — Fractal Information Theory Project  
**ORCID:** 0009-0002-0909-3613  
**Resource type:** Publication — Preprint  
**Version:** 1.0.0  
**Publication date:** 2026-08-01  
**Language:** English  
**License:** CC BY 4.0

## Executive summary

This checkpoint reports three local programs. No network, Firecrawl, external
AI, external audit, laboratory record, remote computation, or empirical data
was used.

1. **P480 — [Computer-assisted proof replayed].** The exact rational O180/O181
   Krawczyk and positivity certificate can now be checked using only the Python
   standard library. The replay recomputes every polynomial and Jacobian
   interval from serialized monomials, applies the rational preconditioner,
   verifies strict inclusion, and pays exact Sylvester and perturbation bounds.
   An adversarial certificate with a zeroed preconditioner is rejected.

2. **P484 — [Strong evidence; exact proof obligation open].** Nine exact
   representation identities reduce the candidate complex tangent to two
   three-dimensional parity sectors. The numerical candidate is extraordinarily
   stable: six independently extracted values of its parameter agree through
   the displayed 80 digits, the Riccati tangent residual is
   (6.65\times10^{-16}), and the fixed-axis balance residual is
   (2.22\times10^{-15}). Nevertheless, an exact counterexample proves that
   active Riccati congruence, positive metric orthogonality, and odd dimension
   alone do **not** force the fixed axis into the causal equal-endpoint plane.
   The repository still lacks an exact identity or interval existence
   certificate showing that its additional shared-entry and standard-sector
   equations pay this final condition. Therefore full-cone optimizer
   nonuniqueness is not marked proved or refuted.

3. **P475 — [Open; resource-bounded no-go].** The trigonometric coefficients
   reduce exactly to the degree-four field
   \(\mathbb Q(\alpha)\), \(\alpha=\sin(\pi/8)\). The 13 Riccati equations are
   cubic over this field and affine in the three normalizer coordinates
   \(A,B,u\). A 45-second, 3-GiB exact lexicographic SymPy elimination did not
   return a univariate relation for the optimum \(L\). This is a reproducible
   limitation of the declared method and hardware envelope, not an
   impossibility theorem and not evidence that \(L\) is transcendental.

The strongest durable positive result of this checkpoint is P480: exact O181
attainment can be independently replayed without a scientific Python stack.
The strongest adversarial result is P484: a previously tempting exact
phase-face conclusion is downgraded because one indispensable algebraic
implication is still absent.

## Inherited exact setting

For the supplied reduced three-slot channel difference \(\Delta\), P473
certified a unique positive root in a rational box of the structured Riccati
system

\[
X_3 N X_3=\frac14\Delta N\Delta .
\]

The corresponding structured normalizer and dual support attain the exact
global value of the declared reduced ordered-comb semidefinite program. P483
extended this existence-and-attainment certificate uniformly to

\[
|q-4/5|\le3\times10^{-9},\qquad
|\theta-\pi/8|\le3\times10^{-9}.
\]

Root uniqueness is local to the certified real polynomial chart. It never
implied uniqueness of the optimizer in the complete complex causal cone.

## P480 — standard-library-only exact replay

### Question

Can a third party replay the finite exact part of O180/O181 without trusting
NumPy, SciPy, SymPy, a nonlinear solver, or precomputed Boolean verdicts?

### Serialized certificate

The builder exports:

- 13 exact polynomial equations as rational monomial lists;
- all 169 exact Jacobian polynomials;
- the rational center and radius \(r=3\times10^{-14}\);
- rational enclosures for the three algebraic sine coefficients;
- a rational approximate inverse \(R\);
- rational center matrices \(N_c\) and \(X_{3,c}\);
- explicit Lipschitz and Sylvester margins.

The checker imports only `fractions`, `json`, `pathlib`, `sys`, `typing`, and
`__future__`. It does not import repository mathematics or numerical packages.

### Exact Krawczyk replay

For \(F:\mathbb R^{13}\to\mathbb R^{13}\), rational center \(x_c\),
box \(X=x_c+[-r,r]^{13}\), interval Jacobian \(J(X)\), and rational
preconditioner \(R\), the checker recomputes

\[
K(X)=x_c-RF(x_c)+(I-RJ(X))(X-x_c).
\]

Every operation is implemented with exact `Fraction` endpoints. The replay
obtains

\[
K(X)\subset\operatorname{int}X,
\]

with minimum inward margin approximately
\(1.8189522605\times10^{-14}\), maximum contraction row sum approximately
\(6.8847446957\times10^{-10}\), and maximum point-residual radius
approximately \(6.6004906310\times10^{-17}\).

The exact certified optimum interval is

\[
\frac{5233281002710117}{10^{16}}
\le L\le
\frac{5233281002710717}{10^{16}},
\qquad |I_L|=\frac{3}{5\times10^{13}}.
\]

### Positivity replay

Exact Gaussian elimination verifies every leading principal minor of

\[
N_c-\frac1{100}I,
\qquad
X_{3,c}-\frac1{100}I
\]

is positive. The declared entrywise variation bounds then give

\[
N\succeq \frac{24999999997}{2500000000000}I,
\qquad
X_3\succeq \frac{124999999997}{12500000000000}I
\]

throughout the root box.

### Adversarial test

The test suite replaces every entry of the preconditioner by zero. The
resulting map can no longer be strictly interior, and the checker rejects the
certificate. Thus acceptance is recomputed rather than copied from a stored
flag.

### Object O186

**Name:** Dependency-Minimal Exact O181 Replay Certificate.  
**Domain:** a finite JSON document of rational monomials, intervals, matrices,
and declared bounds.  
**Codomain:** acceptance or rejection plus exact diagnostic margins.  
**Transformation law:** invariant under loss of the builder stack; its meaning
depends only on rational arithmetic and the ordering of the serialized
coordinates.  
**Necessity:** removing the polynomial/Jacobian data makes independent
recomputation impossible; removing the preconditioner or strict inclusion
test destroys the root certificate; removing positivity payments leaves
attainment inadmissible.  
**Boundary:** the checker validates the finite exact root and positivity
certificate. Riccati-to-support and trace-telescope steps remain analytic
theorems. It does not validate a physical apparatus or observation.

## P484 — adversarial phase-face audit

### Exact block reduction

Let \(\Delta=iK\), with \(K\) real skew, and let the real P473 normalizer be
\(N\). An imaginary Hermitian displacement has the form

\[
N(t)=N+i tQ,\qquad Q^T=-Q.
\]

Since the Riccati equation is linear in \(N\), an exact tangent must obey

\[
X_3QX_3+\frac14KQK=0. \tag{1}
\]

An exact swap/complement basis decomposes the active six-dimensional sector
as

\[
N=N_0\oplus N_0,qquad
X_3=X_+\oplus X_-,qquad
K=\begin{pmatrix}0&C\\-C^T&0\end{pmatrix}.
\]

All nine representation assertions in the executable audit are exact
symbolic identities. Positivity and Riccati congruence imply

\[
X_-=\frac14C^T X_+^{-1}C,
\qquad
B=\frac12X_-^{-1}C^T,
\]

and consequently

\[
BN_0B^T=N_0.
\]

If \(\det B=1\), odd dimensionality supplies a fixed covector and an invariant
skew form. If that covector also lies in the equal-endpoint causal plane, its
lift gives a nonzero exact \(Q\) satisfying (1). In that conditional case the
exact P473 positivity margin yields a nonzero positive phase segment with
constant support and objective.

### Exact falsification of the over-broad lemma

The implication “odd-dimensional metric rotation implies a causal axis” is
false. Set

\[
N_0=X_+=X_-=I_3,qquad
v=(1,0,2)^T,
\]

and take the rational half-turn

\[
B=\frac{2vv^T}{v^Tv}-I_3
=\begin{pmatrix}
-3/5&0&4/5\\
0&-1&0\\
4/5&0&3/5
\end{pmatrix},qquad C=2B^T.
\]

Then \(B\in SO(3)\), \(B^Tv=v\), and both active Riccati blocks hold exactly:

\[
I_3-CC^T/4=0,
\qquad
I_3-C^TC/4=0.
\]

But \(v_1\ne v_3\). Therefore the fixed axis does not satisfy the required
equal-endpoint causality condition. This counterexample does not refute the
specific P473 phase candidate, because it does not impose the additional FIN
shared-entry and standard-sector equations. It proves those equations must do
real work; their effect cannot be replaced by the odd-dimensional argument.

### Evidence retained

At the certified P473 branch the finite calculation gives

\[
\det B\approx1,quad
\sigma(B)=\{1,0.2704387824\pm0.9627371734i\}.
\]

The fixed covector normalized by its first entry is

\[
(1,0.02106502319,0.9999999999999978)^T.
\]

Six tangent equations independently return

\[
\rho\approx0.02979044149309996478959538442548707608969921737694\ldots
\]

through all displayed 80 digits. A 180-digit rerun resolves discrepancies of
order \(10^{-102}\), consistent with the nonlinear-solver tolerance but not a
proof of equality. The candidate tangent residual is
\(6.65\times10^{-16}\).

### Correct status

The exact block representation and exact counterexample are **[Proven]**.
The existence of an exact nonzero causal tangent at the exact P473 root is
**[Strong evidence; open]**. Full-complex-cone optimizer uniqueness is neither
proved nor refuted. O185 is therefore the **Conditional Odd-Dimensional
Complex Comb Phase Face**, not an established physical phase or selector.

## P475 — resource-bounded algebraic elimination

### Exact coefficient field

Let

\[
\alpha=\sin(\pi/8).
\]

Then

\[
8\alpha^4-8\alpha^2+1=0,
\]

and the other coefficients are

\[
\sin(\pi/4)=1-2\alpha^2,
\qquad
\sin(3\pi/8)=3\alpha-4\alpha^3.
\]

The field has degree four: for \(y=2\alpha\), the polynomial
\(y^4-4y^2+2\) is Eisenstein at 2.

After this substitution the exact problem consists of 13 polynomial
equations in 13 operator coordinates over \(\mathbb Q(\alpha)\). Their
operator-coordinate degree is at most three. When \(\alpha\) is represented
as an additional rational-elimination variable, the maximum total degree is
seven, the largest equation contains 29 monomials, and the complete list
contains 262 monomial occurrences.

Every equation is affine in \(A,B,u\). A numerical pivot scout selects rows
5, 10, and 11. Its exact determinant is a nonzero polynomial of total degree
10 with 409 terms, and its value at the floating locator has magnitude
approximately 0.0223899. This identifies a possible structure-aware route but
does not certify that the pivot is nonzero on the entire exact root box.

### Bounded lexicographic attempt

The worker adds the primitive polynomial for \(\alpha\) and asks SymPy 1.12
for an exact lexicographic Gröbner basis in 14 variables, with \(L\) last. Its
hard limits are:

- 40 seconds CPU;
- 45 seconds wall clock;
- 3 GiB address space;
- local Intel i3-10110U / 16 GB host;
- no network and no external CAS.

The parent terminated the attempt after 45.08 seconds. No univariate relation
in \(L\) was returned.

The number \(3^{13}=1,594,323\) is reported only as a crude cubic Bézout
bound over \(\mathbb Q(\alpha)\); it is not a solution count and says nothing
about positivity or the admitted branch.

### Object O187

**Name:** Resource-Bounded Algebraic Elimination Burden Certificate.  
**What it proves:** exact field reduction, polynomial inventory, affine
normalizer structure, a nonzero symbolic pivot polynomial, and faithful
termination at a declared resource limit.  
**What it does not prove:** transcendence, absence of a minimal polynomial,
ideal dimension, algebraic degree, irreducibility of a candidate relation, or
global solution count.  
**Falsification:** a future modular, resultant, homotopy-assisted, or
structure-aware exact workflow can refute the burden boundary by returning and
verifying a univariate relation inside a comparable resource envelope.

## State map after P475/P480/P484

| Lane | Status | Exact boundary |
|---|---|---|
| O181 root and O167 attainment at the nominal point | Computer-assisted proof | Declared reduced three-slot ordered comb |
| Uniform P483 parameter tube | Computer-assisted proof | Tiny dimensionless rectangle and one common local box |
| Dependency-minimal replay | Computer-assisted proof replayed | Standard-library checker validates finite root/positivity data |
| Full complex optimizer uniqueness | Open; strong evidence against | Exact causal-axis identity missing |
| Algebraic minimal polynomial of \(L\) | Open | Direct exact lexicographic method exceeded the declared envelope |
| Lean formalization | Source only | No local Lean toolchain; not machine checked |
| Selector / QW-2191 | Open | No non-premise strict source |
| Dimensional scale | Missing | No time, length, action, mass, or energy unit generated |
| Physical realization | Blocked by external evidence | No apparatus, calibration, custody, or laboratory record |
| Legacy-to-strict completion and role transfer | Incomplete | No silent substitution or inherited physical role |

## Kernel, selector, dimensional, and physical boundaries

The kernel split remains

\[
K_{\mathrm{legacy,ont}}(d)=
\frac{\alpha_{\mathrm{geo}}\cos(\omega d+\phi)}{1+\beta_{\mathrm{tors}}d},
\qquad
K_{\mathrm{strict,gate}}(d)=
\frac{\cos(\omega d+\phi)}{1+\beta d^\eta}.
\]

The present programs are downstream finite-operator analyses. They derive
neither kernel, provide no completion map, and transfer no electroweak,
electromagnetic, gravitational, or other legacy physical role.

QW-2191 remains open. The complex phase-face candidate is not a selector,
orientation source, clock, preparation, instrument, environment, apparatus,
observer, or measurement record. All numbers remain dimensionless. The dual
functional calculus \(e^{-itA}\) and \(e^{-tA}\) remains mathematical
wave/heat duality, not a laboratory implementation.

The internal ontology remains

\[
\text{nadsoliton}\to\text{light}\to\text{matter}\to
\text{emergent observer},
\]

with no informational substrate below the nadsoliton. Nothing in this
checkpoint upgrades that ontology into experimentally validated physics.

## Recommended next programs

| Rank | Program | Exact question | Estimated chance of a decisive local result |
|---:|---|---|---:|
| 1 | P485 | Prove or refute, by exact ideal identity or certified counterexample, that the P473 shared-entry and standard-sector equations force the causal equal-endpoint axis | 55% |
| 2 | P486 | Certify \(\det C>0\), \(\det B=1\), and every orientation premise uniformly on the P473 box | 80% |
| 3 | P487 | Eliminate \(A,B,u\) through a certified pivot and attempt a ten-variable structure-aware resultant/Gröbner calculation | 50% |
| 4 | P488 | Use finite-field Gröbner scouting and Chinese-remainder reconstruction to seek a candidate polynomial for \(L\), followed by exact substitution | 35% |
| 5 | P489 | If P485 succeeds, certify persistence of the complex phase face throughout a nonzero sub-tube of P483 | 30% |
| 6 | P490 | Construct a second implementation of the P480 checker and differential-test both exact verdicts | 75% |
| 7 | P491 | Bound the local algebraic degree of the isolated O181 branch without computing the full minimal polynomial | 45% |
| 8 | P492 | Certify the exact dimension of the full complex contact face or prove a transverse rank lower bound | 40% |
| 9 | P493 | Derive interval-certified sensitivities \(dL/dq\) and \(dL/d\theta\) from the implicit Riccati system | 75% |
| 10 | P494 | Grow the P483 parameter cover adaptively until an exact positivity, inclusion, or conditioning boundary is reached | 70% |
| 11 | P495 | Test whether the exact two-/three-slot support reductions suggest a bounded four-slot theorem or a counterexample | 30% |
| 12 | P496 | Formalize the minimum operational data needed to turn O167 into an observable and prove the local no-laboratory-data boundary | 85% |
| 13 | P497 | Audit which O167 inputs are kernel-split robust and which require an explicit legacy-to-strict completion atom | 70% |

The selected next bounded batch is

\[
\boxed{P485\longrightarrow P486\longrightarrow P487}.
\]

P485 has priority because it decides whether P474/P484 becomes an exact
nonuniqueness theorem or remains numerical evidence. P486 pays a necessary
orientation premise independently. P487 exploits the exact affine structure
found by P475 instead of repeating the exhausted blind lexicographic attempt.

## Reproducibility

From the repository root:

```bash
MPLCONFIGDIR=/tmp/matplotlib python3 fin_program_484.py
MPLCONFIGDIR=/tmp/matplotlib python3 fin_program_480_build_certificate.py
python3 fin_program_480_standalone_checker.py
MPLCONFIGDIR=/tmp/matplotlib python3 fin_program_475.py
python3 -m unittest -v test_fin_programs_480_484.py
```

P475 deliberately consumes up to 45 seconds and is expected to report the
bounded timeout on the present host. The P480 replay itself does not require
the scientific Python stack.

## Artifact inventory

- `FIN_Local_Research_Checkpoint_P475_P480_P484_EN.md`
- `FIN_Local_Research_Checkpoint_P475_P480_P484_EN.tex`
- `FIN_Local_Research_Checkpoint_P475_P480_P484_EN.pdf`
- `FIN_Local_Research_Checkpoint_P475_P480_P484_State.json`
- `FIN_Program_475_Results.json`
- `FIN_Program_475_Algebraic_Inventory.csv`
- `fin_program_475.py`
- `fin_program_475_elimination_worker.py`
- `FIN_Program_480_Standalone_Certificate.json`
- `FIN_Program_480_Standalone_Check_Result.json`
- `fin_program_480_build_certificate.py`
- `fin_program_480_standalone_checker.py`
- `FIN_Program_484_Results.json`
- `FIN_Program_484_Symmetry_Phase_Face.csv`
- `FIN_Program_484_Phase_Face_Witness.npz`
- `FIN_Program_484_Figures/p484_exact_complex_phase_face.png`
- `fin_program_484.py`
- `test_fin_programs_480_484.py`
- `RELEASE_10_46_PROGRAMS_475_480_484.md`
- `FIN_Checkpoint_P475_P480_P484_AGENTS_Guardrail.txt`
- `FIN_PROGRAMS_475_480_484_RELEASE_MANIFEST.sha256`

## Final conclusion

Release 10.46 improves rigor in two opposite but complementary ways. It makes
the exact O181 certificate much easier to audit, while withdrawing an exact
phase-face claim that the available artifacts did not yet prove. The current
mathematical frontier is sharply localized: either the additional FIN
shared-entry equations force the observed causal axis exactly, or the
near-perfect numerical phase direction is only an approximation. That is now
a finite, falsifiable local problem. No physical interpretation is needed to
state or decide it.
