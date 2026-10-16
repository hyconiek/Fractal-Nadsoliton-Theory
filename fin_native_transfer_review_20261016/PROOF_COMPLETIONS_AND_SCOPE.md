# NT proof audit: actual stationary input, complex reflection and operator scope

Date:2026-10-16. Analytical model audit, not mechanized or peer-reviewed proof.
Source files and all earlier accepted packages are preserved unchanged.

## 1. The missing historical input is available for the current theorem

The source executor could not inspect the real post-review completions.
The current intake can. The exact dependency is:

1. MT supplies the connected D_N, reversible count law pi_D,N with
   multinomial multiplicity, local critical static convergence and moments.
2. OBS/PROOF_COMPLETIONS_AND_CORRECTIONS.md §2 explicitly extends UF to
   the random local-equilibrium initial law using those fast/soft moments,
   localization and nonexplosion. It also pays the reflected/full-generator
   departure difference on each fixed critically scaled horizon.
3. This supplies stationary reflected joint laws at (0,t) and (0,2t), not
   merely the original iid-start UF theorem with Z0=0.
4. The quartic limit is stationary reversible with rho_lambda. The exact
   spectral gain defines g_N; the fixture's binary64 g_c is not used as the
   asymptotic parameter sequence.

The actual files and saved predecessor inputs are hash-bound in
[verification.json](verification.json). This pays the current mathematical
dependency. It does not attest historical remote reading, uploaded-plan byte
identity or a checked-out producer Git environment.

Keep the three processes distinct: full global Gibbs, stationary local
reflection, and unchanged late iid-start parent. Only the second supplies
the critical stationary L2 operator family used by NT-002.

## 2. Fix the reflection/conjugation notation

Use either the linear pullback R F=F composed with time reflection, or the
antilinear operation Theta F=conj(R F). The correct reflection form is

    B(F,G)=E[conj(R F) G]=E[(Theta F)G]
          =<E[F|X0],E[G|X0]>_pi.

The source describes Theta as conjugating and later writes an additional
conjugation in its Gram notation. Read literally that double conjugation is
wrong. For F=i constant it would give -1 instead of1. Its primary displayed
formula E[conj(F composed with reflection)F] and its complex code use the
correct form. This is a notation correction, not a counterexample to RP.

Past/future independence is conditional on the **full Markov count state**,
not the observed parity. The finite stationary Markov property and detailed
balance prove the Gram identity for all bounded future cylinders. Include
time-zero functions, or approximate them by stochastic continuity, to identify
the full quotient with L2(pi_D,N).

Time translation acts through the full native semigroup. Positivity here
means a self-adjoint nonnegative transfer semigroup; it is not an assertion
that every matrix representation has positive entries or that its vectors
are physically preparable quantum states.

The intake independently compares the conditional Gram calculation with a
chronological four-time expectation, including complex functions. This tests
the conjugation and order conventions rather than constructing a Gram matrix
and merely checking its defining positivity.

## 3. Strong core and positive resolvent

The source's three-term norm proof is valid. J_N f=f(z_N) is defined on
C_b(R), not on arbitrary L2 equivalence classes. For h=T(t)f,

    ||T_N(t)J_N f-J_N h||_N^2
      =<J_N f,T_N(2t)J_N f>_N
       -2 Re<J_N h,T_N(t)J_N f>_N+||J_N h||_N^2.

The paid stationary two-time laws and static marginal convergence apply to
all three bounded continuous terms. Reversibility identifies the first
limiting term with ||T(t)f||^2. Their combination vanishes. No exact finite-N
parity lumpability, conditional independence of fast modes or new fast gap
enters the proof.

Bounded continuity of T(t)f has an elementary payment here. Couple two
nonexplosive scalar diffusions with the same Brownian path. Their difference
obeys the one-sided Lipschitz estimate

    |Z_t^x-Z_t^y|<=exp(a t)|x-y|, a=chi lambda,

because b(x)-b(y)=(x-y)[a-(x^2+xy+y^2)/3]. Nonexplosion follows from the
quadratic Lyapunov estimate. Dominated convergence for bounded f then gives
T(t)f in C_b. No false sup-norm strong continuity on all C_b is needed.

For each fixed r>0, the source's contraction/Laplace argument pays the
lifted resolvent. This is not a zero-parameter inverse, uniform-in-time
bound, norm-resolvent theorem on one fixed space or a complete spectral limit.

### Dense-section completion, without evaluating null-set representatives

C_b is dense in L2(rho), and ||J_N f||_N->||f|| for f in that core.
Define a strongly convergent section u_N->u by approximation with J_N f_j:
f_j->u in L2(rho) and

    lim_j limsup_N ||u_N-J_N f_j||_N=0.

Native/limit contraction, the core theorem and T(t)f_j in C_b give the
corresponding convergence of T_N(t)u_N by a three-epsilon argument. The same
argument works for fixed positive resolvents with their bound1/r. No choice
of an arbitrary L2 representative at an atom z_N is required.

This supplies a useful strong-section interpretation. The source correctly
does not export a complete independently packaged Kuwae–Shioya/Mosco
application. No standalone full-Mosco or compact-spectral claim is needed
for the admitted core/resolvent theorem or the next physics question.
Do not turn an optional theorem-label upgrade into a default campaign.

## 4. Critical form/domain and time-ordered correlators

The form closure and ground-state transform are sound. With
U=z^4/12-a z^2/2 and rho proportional to exp(-U), use the closure of
integral |f'|^2 rho on C_c^infinity. Smooth positive multiplication by
sqrt(rho) maps its core to the corresponding compactly supported smooth
Schrödinger-form core. The transformed form is ||(partial+U'/2)psi||^2.

The differential expression is

    H=-partial^2+z^6/36-a z^4/6+(a^2/4-1/2)z^2+a/2.

At lambda0 it is -partial^2+z^6/36-z^2/2. To pay the source's compressed
ground-state domain statement, approximate1 in the weighted form by smooth
cutoffs chi_R with derivative at most C/R. Quartic tails give norm convergence
and E(chi_R,chi_R)->0. Thus1 lies in the form kernel and sqrt(rho) in the
operator kernel after transformation. Positivity follows from the closed
Q*Q form even where the displayed potential is negative. A zero-energy
weighted Sobolev function has derivative0 almost everywhere, so the ground
space is one-dimensional. Confinement of this individual limit operator
does not prove compact convergence of the native family.

The source's time-ordered multiplication identities hold for bounded tests
and their complex versions. Products of passive functions can all be realized
on the same path probability space. Noncommutation between multiplication
and transfer does not supply incompatible physical measurements.

## 5. Scope actually admitted

The finite stationary RP/transfer theorem, strong lifted core and positive
resolvent, form-realized critical operator and bounded-cylinder dictionary
are admitted under the inherited model/spectral/local-regime premises.
Many arguments already appear in the preceding research submission; the NT
execution and this intake specialize, bind and review them, not discover
the general Markov/RK construction anew.

The original iid parent inherits bounded fixed-window predictions through
the existing labelled-TV/exit coupling. It does not inherit exact stationarity,
a full-space L2 operator theorem, arbitrary unbounded moments or the global
Gibbs spectrum. Source N2 numerical fixtures are ordinary floating checks
outside the main N in12N asymptotic regime. The intake additionally uses an
N3 parity diagnostic: N2 has exceptional parity closure and must not be used
as a generic finite-N nonlumpability witness.

Source NT-006 is useful scoped reporting, not an automatic architecture
authorization. The next architecture decision belongs to this intake. No
physical QM, kernel source, conserved energy, geometry or empirical status
is promoted by the present success.
