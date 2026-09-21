# R7P-105--112 — weighted Hodge, finite-N walkers, passive memory and firewall

## Weighted Hodge structure (R7P-105)

For a connected weighted graph with oriented incidence `D` and positive diagonal
edge matrix `G`, put `A=D^T G D`.  Then

`P=G^(1/2) D A^+ D^T G^(1/2)`

is the orthogonal projector onto `Range(G^(1/2)D)`.  Therefore

`Sigma_tree=G D A^+ D^T G`,
`Sigma_cycle=G-Sigma_tree`

are PSD, have ranks `n-1` and `m-n+1`, respectively,
`D^T Sigma_cycle=0`, and `D^T Sigma_tree D=A`.  Reorienting edges by a diagonal
sign matrix S replaces D by SD and conjugates both edge covariances by S.  Node
relabelings likewise act by permutations.  Nothing in this algebra refers to
crossings in a planar drawing.  An exact rational weighted-triangle fixture is
included in the checker.

For the 12-label complete weighted graph this gives ranks 11 and 55.

## Cycle spectrum (R7P-106)

Fresh diagonalization reproduces rank 55 and 28 positive numerical clusters at
tolerance 1e-8, with the saved multiplicity pattern.  The edge-space covariance
commutes numerically with the exact D12 edge action to roundoff, explaining the
prevalence of real two-dimensional representation degeneracies.  The value 28
is retained as a numerical spectral count, not promoted from rank 55 to an exact
multiplicity theorem.

## Exact finite-N walker generator (R7P-107)

For occupation counts `n_i` with sum N, independent continuous-time walkers have

`n -> n-e_i+e_j` at rate `n_i W_ij`.

For `p=n/N`, direct generator calculation gives

`E[dp|p]/dt = -A p`,

and predictable quadratic variation

`(1/N) sum_{i<j} W_ij(p_i+p_j)(e_i-e_j)(e_i-e_j)^T dt`.

Only at `u_i=1/12` does this equal `A/(6N) dt`.  The invariant law of counts is
Multinomial(N,u), so

`Cov(p)=(diag(u)-u u^T)/N`,

and every orthonormal mean-zero mode has variance `1/(12N)`.  The N=2 generator
(78 states) is built explicitly and verifies row sums and stationarity.

## OU comparison (R7P-108)

The equilibrium linear-noise OU model gets the exact leading covariance right,
but is not the jump process.  For the normalized alternating mode, a single
uniform walker is a symmetric two-point variable.  Hence the empirical-mode
standardized skewness is exactly zero and its excess kurtosis is `-2/N`, while
the OU Gaussian has zero excess.  At N=12,48,192 these are -1/6, -1/24 and -1/96.
Thus the non-Gaussian diagnostic vanishes at the expected 1/N rate while the
mode variance is exactly `1/(12N)`.

## Even/odd memory and passivity (R7P-109--110)

With even sites observed and odd sites hidden, write the full PSD Laplacian block
matrix as `[[A_EE,B],[B^T,D]]`, with `D>0`.  Since D is real symmetric,

`B(zI+D)^(-1)B^T = sum_alpha R_alpha/(z+d_alpha)`,
`R_alpha=B P_alpha B^T >=0`.

Fresh numerical residue ranks are 1,2,2,0 at four distinct hidden pole groups,
so the visible degree is 5.  The last hidden one-dimensional alternating mode
is exactly uncoupled by cyclic reflection pairing.  Nonzero ranks of the other
residues are presently a numerical realization diagnostic rather than a strict
transcendental interval theorem.

For `Re z>0`, every scalar `1/(z+d_alpha)` has positive real part; hence the
matrix memory is positive-real/Stieltjes.  At z=0, the Schur complement
`A_EE-BD^(-1)B^T` is PSD because it is a Schur complement of PSD A with D>0.
This proves passivity in the supplied symmetric class.  It says nothing about a
changed driven, nonnormal or negative-loading law.

## Conditional finite-copy Gibbs model (R7P-111)

Take N labelled copies with uniform reference measure and tilt by

`exp[(g/(2N)) sum_{a,b} A7[x_a,x_b]]`.

For empirical p the exponent is exactly `N g p^T A7 p/2`.  Sanov/multinomial
counting therefore gives rate function `D(p||u)-g p^T A7 p/2` up to the global
normalization constant.  A pair-only a!=b convention differs by an O(1)
finite-N self term but has the same N-speed rate function.  This is a
conditional realization with supplied `g=beta J`; it does not derive beta, J,
the mediator rank, a pump or a clock.

## Interpretation firewall (R7P-112)

Proved passive objects remain passive.  Finite-N noise does not supply active
gain; an algebraic cycle space is not a newly sourced information substrate;
the conditional Gibbs realization does not source its coupling; and static
Schur passivity is not a universal no-go for every future driven completion.
The gain/source problem remains open unless a separate sourced active law is
introduced.
