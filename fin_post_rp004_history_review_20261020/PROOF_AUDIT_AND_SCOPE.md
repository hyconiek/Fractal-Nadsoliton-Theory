# Admission of the five studies, with explicit corrections and reuse

Date: 2026-10-20. Analytical intake, not formal verification, a fresh
literature-priority survey or independent external peer review.

## P1 — chosen heat law, exact hidden memory and provider type

Accept the even/odd block decomposition of the supplied symmetric Z12
Laplacian, the exact variation-of-constants memory equation and retained
hidden-initial-data term. D is positive definite without relying on the
printed eigenvalues: for any nonzero even-coordinate vector, extend it by
zero on odd vertices. The full connected positive Laplacian has only constant
vectors in its kernel, and this extension cannot be a nonzero constant.
Equivalently D has a strictly positive diagonal-dominance margin equal to
the sum of the positive odd-shell couplings.

The circulant C has exactly zero alternating Fourier coefficient. Thus its
alternating retained projection is autonomous under this **chosen linear
law**, but the whole six-dimensional retained state is not. The full hidden
initial term cannot be removed from a general preparation. The source's
memory-positivity statement is correct because `exp(-tD)` is positive as an
operator, and congruence by C preserves PSD.

Accept the short-lag diffusion/coherent-endpoint distinctions, but the
`|exp(-itA)|^2` readout assumes a Born-type endpoint rule and is not generally
a Markov semigroup. Accept conditional Gamma_i reconstruction. **Correct the
claim about stationary averaging** using [MEAN_BRACKET_CORRECTION.md](MEAN_BRACKET_CORRECTION.md):
the full mean tensor does identify a reversible one-hot generator with known
pi. A scalar/diagonal-only average, nonreversible current or count-process
average is not the same data.

## P2 — normalized KL functional, not a unique dynamics source

Accept Theorems I--IV in their stated positive C3-curve and finite irreducible
heat-law class. The uniform lower bound on each p_i pays the Taylor remainder;
each normalized summand has O(Delta^2) error, hence total O(Delta).
The `u=2sqrt(p)` map gives a radius-two sphere, so the coefficient
`S_min=2 arccos(sum sqrt(p0_i p1_i))^2/T` is correct. Strictly positive
endpoints have a unique short arc that stays in the positive orthant.

The strict nonminimality argument is sound: a minimizing energy geodesic has
constant speed on the given time interval, while the heat Fisher speed is
analytic for all positive times and tends to zero. Constancy on one open
interval forces it to vanish everywhere, contradicting nonstationary data.
It does not rule out another kinetic geometry or a different variational law.

Accept the entropy-gradient linearization `-I` and the spectral mismatch in
the stated non-flat Laplacian family. Nonuniform weights alone do not imply
all six eigenvalues are distinct in every radial family; the displayed
tuple's non-flat spectrum, or unequal eigenvalues paid for a particular
tuple, is the relevant premise. Roundoff tests do not certify every printed
digit of transcendental kernel values.

The logarithmic-mean representation is **reused**: AGENTS ST8397--ST8547
already gives the identity and its source-law caveat. The proof here is a
correct edgewise algebraic specialization; no fresh external theorem
application or primary-paper reading is asserted by this intake. Likewise
the same-escape activity comparison does not rediscover the older RC/AC
obstructions. The pure-KL scaling/geodesic test is admitted as a scoped
completion/test of GEN-654, not its automatic ontological meaning.

## P3 — six labelled correlations and its precise minimality

Accept H1--H4 for the positive homogeneous D12 radial one-label CTMC with a
supplied comparison lag and labelled Fourier observables. The antipodal
shell has one, not two, neighbors. Derive the inverse independently from
the complete Fourier inversion:

    w_d=-(1/6)sum_(k=1)^5 lambda_k cos(pi*k*d/6)
                          -(1/12)lambda_6(-1)^d.

This matches the displayed inverse. Its maximum absolute row sum is 11/12;
the last row attains it and each other row is smaller. Thus the logarithm
mean-value error bound and positivity-cone test are correct. Rank-six
invertibility leaves a nonzero null direction after any one row is removed;
the positive orthant is open, paying the local counterfamily.

The six-data minimality concerns omission of members of **these specified
modal records**, not arbitrary encodings or every conceivable observable.
Single-lag cone membership supplies a compatible generator in the class;
it does not prove that an independently supplied actual process is Markov
or has that generator at all lags. Semigroup-looking selected correlations
alone do not generally prove full path-law equality outside the class.

## P4 — actual-change indexing and an all-lag return counterpair

Accept E1/E1a, given the declared constant exit rate s. The construction of
the CTMC can be written as iid exponential(s) holds independent of the
stationary jump skeleton; it pays the actual-jump Poissonization without
confusing a null update attempt with a change. One-step modal data identify
`w/s` and the full event-index law is invariant under common rate rescaling.

Accept the exact Q(sqrt(3)) witness pair. Its two spectra swap double modes
1 and 2 and keep the simple mode 6, so all continuous and event-indexed
return functions agree. The directed-edge cubic statistic is invariant
under every vertex permutation; its nonzero per-vertex gap 11/24000 proves
nonisomorphism, stronger than just a non-affine shell rearrangement.

The 120-permutation test at printed decimal weights remains a finite floating
diagnostic. It is not a certified assertion about the exact kernel, a global
classification of every cospectral radial operator, or a neighborhood no-go.

## P5 — four observed states, eleven moments and what remains hidden

Accept E3/E4. Four consecutive observed states `Y0,...,Y3` involve **three
actual changes**, not four jumps. No adjacent equality is allowed, so exactly
five equality partitions remain. Rooted exact enumeration supplies the same
five probabilities without using the source's case classifier. The equal
return moments and different alternating cubic moment give

    Delta Pr(ABAB)=11/34629528,
    TV(full four-state equality laws)=11/17314764.

Minimality is for this exact pair and this consecutive equality observation,
not for all generators or a physical instrument. It does not give iid
repetitions, usable detection performance or native N2 observation costs.

For r transitions, an alternating record of r+1 states has probability
`sum_(j!=0) J_0j^r` by symmetry and translation invariance. Newton identities
through power 11 recover the outgoing eleven-element multiset. The antipodal
value is the sole value of odd multiplicity even when weights collide.
The numerical value does not thereby acquire a spatial shell label. Twelve
states are a sufficient bound, **not a proved optimal record length** for
the six-parameter radial subclass. The rational shell swap keeps every
alternating moment but changes a nonalternating return by -11/7776.

Accept the continuous-time four-point cubic coefficient under its separate
common-clock contract. The d=0 contribution cancels because the pair's heat
traces coincide; finite-dimensional Taylor expansion controls the rest. It
is not an instantaneous native interacting result or a physical measurement.

## Status not changed by this bundle

It contains no independent full-text audit of RP004's new moving-cone proof.
Those actual proof files are now locally available, so do not recommend
acquiring them again. Their status stays SUBMITTED_CONDITIONAL_THEOREM /
INDEPENDENT_AUDIT_PENDING, as in the current execution overlay. This intake
does not execute V-RP4-001. F01/F03/F06, natural source, clock/unit, physical
incidence/composition, selector/legacy roles, B4/FIN-M0 and QM/GR/ToE gates
are unchanged. None of the five comparison carriers is silently promoted
to the primordial informational nadsoliton or to physical spacetime.
