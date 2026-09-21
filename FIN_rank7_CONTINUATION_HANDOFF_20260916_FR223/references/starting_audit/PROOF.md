# Audited mathematical results

Date: 2026-09-13. Input: the supplied local pre/post-Discord handoff and its
41 artifacts. This document, not the imported handoff's status labels, states
the accepted scope. No internet conversation was used for this audit.

## 1. Objects and coordinate conventions

Use the frozen strict matrix
`W_ij = cos((743 d+650)/4000)/(1+d^(9/5))`, with zero diagonal and cyclic
distance `d`. Let `A=diag(W 1)-W`, with Fourier eigenvalues `lambda_k`.
The rank-seven mediator retains both real components of sectors 3, 4, 5
and the alternating sector 6. Its factor `X` has columns

`sqrt(lambda_k/6) cos(2 pi k j/12)`,
`sqrt(lambda_k/6) sin(2 pi k j/12)`, for k=3,4,5,
and `sqrt(lambda_6/12)(-1)^j`.

Then `A7=X X^T`. The four-column cosine submatrix is denoted `C` here.
It is a restriction of `X`, not an exhaustive replacement for `X`.
For amplitude coordinates `s`, the supplied field is `h=C s` and
`J_k=sqrt(lambda_k/6)s_k` for k=3,4,5;
`J_6=sqrt(lambda_6/12)s_6`. These definitions remove the handoff's ambiguous
switches between Fourier amplitudes, fields, and dual coordinates.

## 2. Accepted pre-cutoff structural results

### Edge decomposition

Let `D` be the 66-by-12 oriented incidence matrix and `G` the positive
diagonal matrix of edge weights. Then `A=D^T G D`. If `B=G^(1/2)D`,
`B(B^T B)^+B^T` is the orthogonal projection onto `Range(B)`, of rank 11
because the graph is connected. Therefore

`Sigma_cycle = G-G D A^+ D^T G`

is positive semidefinite of rank 55, `D^T Sigma_cycle=0`, and its
complement `Sigma_tree=G D A^+D^T G` satisfies
`D^T Sigma_tree D=A`. Edge orientation changes transform these covariances
by signed permutations; planar crossings do not enter their definition.
The 28 distinct positive eigenvalues and printed residual magnitudes are
numerical observations, not exact multiplicity proofs in this intake.

### Independent walkers and passive memory

For N independent continuous-time walkers with symmetric jump rates W, the
empirical distribution has drift `-A p` and conditional quadratic variation

`(1/N) sum_{i<j} W_ij (p_i+p_j)(e_i-e_j)(e_i-e_j)^T dt`.

Only at uniform p is this `A/(6N) dt`. The Gaussian Ornstein--Uhlenbeck
modal equation is the equilibrium linear-noise description, not the exact
finite-N jump process. Exact equilibrium empirical covariance is
`(diag(u)-u u^T)/N`; every orthonormal mean-zero mode has variance `1/(12N)`.
There is no spontaneous active negative stiffness in this model.

The even/odd Schur complement is PSD since A is PSD and its odd principal
block is positive definite. In particular, eliminating that block cannot
generate a negative static stiffness. The rank-five memory and its trace
are independently reproduced numerically. This is consistent with ST293;
it is not a no-go for unspecified driven or nonnormal systems.

### Spectral budget, distinguishability, and duality

For a D12-invariant mediator B with B1=0,
`V_g(e_j)-V_g(u)=log(12)-g Tr(B)/24`.
The real sector dimensions are 2,2,2,2,2,1. On each irreducible sector a
D12-commuting real symmetric B is scalar; `0<=B<=A` bounds its scalar by
the corresponding eigenvalue. Enumerating the 64 sector subsets shows the
largest budget at rank at most six is `2(lambda_3+lambda_4+lambda_5)`;
at rank seven add `lambda_6`. Exact rational spectral/logarithm intervals
place `6 log(12)` strictly between these values. This proves only the
vertex-versus-uniform energetic threshold, not a global transition or
minimal information dimension. Sector 5 already distinguishes all labels
because 5 is coprime to 12; in the cumulative ladder it appears at rank 3.

For g>0 complete the square in

`F(p,theta)=D(p||u)+||theta||^2/(2g)-theta^T X^T(p-u)`.

Minimizing over theta gives V_g; minimizing over p gives
`Phi_g(theta)=||theta||^2/(2g)-log(sum exp(X theta)/12)`.
These are two orders of joint **minimization**, not a minimax exchange.
For each supplied theta the unique entropy completion is `p=softmax(X theta)`.
At stationarity `theta=g X^T p`. The nonlinearity generates Fourier
frequencies outside the retained field; it does not impose
`p-u in Range(A7)`. Finite addition modulo 12 verifies the reported
frequency-generation orders (a generic algebraic possibility, not nonzero
amplitude for every special input).

## 3. Exact CRT and Ising identities

Set `alpha=pi(j mod 4)/2`, `beta=-2 pi(j mod 3)/3`. Exact trigonometric
identities on the twelve labels give observables
`cos(alpha)`, `cos(beta)`, `cos(alpha-beta)`, `cos(2 alpha)` for modes
3,4,5,6 respectively. The restriction `J_k>=0` is a supplied locked sector.
CRT alone does not prove a global phase-reduction or all mixed-covariance
signs; no universal cooperativity theorem is promoted in this intake.

On the limiting even support (J6 tends to positive infinity), write
`Aspin=+/-1`, `Bspin=1/4+3Y/4`, Y=+/-1. The observables become
`Aspin, Bspin, Aspin Bspin`. Y=-1 has multiplicity two and Y=+1 has
multiplicity one. After summing this degeneracy the four probabilities have
Ising weights with

`H_A=J3+J5/4`, `H_Y=3J4/4-(log 2)/2`, `K=3J5/4`.

For strictly positive probabilities, taking log ratios proves exactly

`K>=0 <=> p1 p4>=p2 p3`,
`H_Y>=-(log 2)/2 <=> 4 p1 p3>=p2 p4`,
`H_A>=K/3 <=> p1 p2^2>=p3 p4^2`.

At zero probabilities these polynomial inequalities are necessary closure
constraints; no converse boundary parametrization is asserted here.
The determinant of the covariance of the three unscaled observables is
`81 p1 p2 p3 p4`, verified symbolically. Scaling gives exactly
`3 lambda3 lambda4 lambda5 p1 p2 p3 p4/8`. These identities do **not**
prove the missing global Ising curvature inequality.

## 4. Parity and the correctly qualified Schur reduction

For the four-coordinate family let q be the even-sector probability.
The law of total covariance gives

`M=W_par+b b^T`,
`W_par=q C_plus+(1-q) C_minus`,
`b=sqrt(q(1-q))(mu_plus-mu_minus)`.

The fourth conditional coordinate is constant within each sector. Thus
W_par has zero fourth row and column, and
`b6^2=lambda6 q(1-q)/3`.
Put `sigma=[2 lambda3(lambda4+lambda5)-lambda4 lambda5]/(24 lambda3)`.
The exact spectral certificate gives `sigma>lambda6/12`, so
`eta=1-b6^2/sigma>0` for all q. Block Gaussian congruence gives

`n_-(sigma I4-M)=n_-(sigma I3-Mtilde)`,
`Mtilde=W_par[345]+b[345] b[345]^T/eta`.

This inertia equality, including equality cases, is the safe general
statement. The scalar resolvent `S=1-b^T(sigma I-W_par)^(-1)b` is defined
only when the inverse exists. If W_par has **exactly one** eigenvalue above
sigma and none equal, S>0 preserves that count and S<0 increases it to two.
S=0 requires a separate equality treatment. Without that premise, the
handoff's unconditional `S>0 iff lambda2(Mtilde)<sigma` is not licensed.

## 5. Corrected one-dimensional resolvent certificate

On s4=s5=s6=0 put `a=lambda3/6`, `J=sqrt(a)s3`, and

**`r=sech(J)`, not `r=exp(-J)`.**

Direct summation gives `q=1/(1+r)`, `t=tanh(J)=sqrt(1-r^2)`,
`W33=a r^2/(1+r)`, `b3^2=a r(1-r)/(1+r)`,
`b6^2=lambda6 r/[3(1+r)^2]`. The other parity-mean differences vanish.
Hence

`S(r)=1-a r(1-r)/(sigma(1+r)-a r^2)
       -[lambda6/(3 sigma)] r/(1+r)^2`.

Let `c=lambda6/(3 sigma)`, `D0=sigma(1+r)-a r^2`. Its endpoints are
positive and D0 is concave, so D0>0 on [0,1]. The numerator

`N=(1+r)^2 D0-a r(1-r)(1+r)^2-c r D0`

is a cubic. Reconstructed Bernstein coefficients on [0,1/2] and [1/2,1]
have strictly positive **rational interval lower endpoints**. This proves
the corrected scalar expression positive on the entire closed interval
for the strict spectrum, not just the supplied decimals. At points where
another W_par eigenvalue equals sigma, the full inverse is not defined;
the scalar expression there is a continuous reduced expression, not an
assertion that a singular matrix has an inverse.

The numerical face minimum is S~0.057549460989 at r~0.6608385948 and
s3~1.703819057. Uniqueness of that minimum has **not** been upgraded here:
the imported claim of a Sturm count has no supplied Sturm certificate.

## 6. Strict-spectrum extreme-face curvature certificate

Now s4=s5=0 and s3,s6>=0. Write `t=tanh(sqrt(a)s3)` and x=q t.
For fixed t, s6>=0 implies
`q>=1/(1+sqrt(1-t^2))`, equivalently `x^2<=2q-1` and q>=1/2.
The Schur-reduced matrix splits into the scalar

`f=a[q+x^2(c(1-q)-1)/eta]`, `eta=1-cq(1-q)`,

and a 2-by-2 block with diagonal entries lambda4/12, lambda5/12 and
off-diagonal entry `sqrt(lambda4 lambda5)x/12`.
Its larger eigenvalue is the handoff's B(x); the smaller is at most
`min(lambda4,lambda5)/12<sigma`.

Let `tstar^2=1-sigma/a`. Symbolic substitution verifies `B(tstar)=sigma`.
If x<=tstar, the two block eigenvalues are at most sigma, so at most the
scalar f can exceed sigma. If x>=tstar, maximize f over physical x^2:

- For `qmin=1-sigma/(2a) <= q <= qcrit=1-1/c`, its coefficient is
  nonnegative: use `x^2=2q-1`.
- For `qcrit<=q<=1`, use the smallest allowed `x^2=tstar^2`.

The polynomial `(sigma-f)eta` has positive interval Bernstein coefficients
on the first interval. On the second its exact factor `(1-q)` is removed
symbolically **before** interval evaluation, and every Bernstein coefficient
of the quotient is strictly positive. Thus f<sigma for q<1 in this case;
the endpoint equality is q=1, x=tstar. Together with Schur inertia this proves

**`lambda2(M)<=sigma` and `lambda2(Mtilde)<=sigma` on this extreme face.**

Both use eigenvalues in decreasing order. Equality occurs only at the
compactified boundary q=1, x=tstar. No Ising-global lemma or numerical
off-face search is used in this face proof. Its exact proof dependencies
are the existing strict spectral enclosures, symbolic identities, rational
interval operations, and Bernstein positivity. This is not a proof-assistant
formalization. Recomputable endpoint fractions are in results.json.

## 7. Counterexample to transfer to the entire seven-coordinate landscape

The handoff's proposed global bound on the four-amplitude covariance must
not be stated for the full seven-dimensional dual. In fact the stronger
everywhere Hessian bound is false.

Take `h_j=2 cos(pi j/2)` and `p=softmax(h)`. This is a legitimate finite
rank-seven field. Let `x=E_p cos(pi j/2)=tanh(1)>3/4`.
In the full seven-coordinate covariance the cosine (4,5) block is

`(1/12) [[lambda4, sqrt(lambda4 lambda5)x],
         [sqrt(lambda4 lambda5)x, lambda5]]`;

the sine (4,5) block has the opposite off-diagonal sign. Cross covariances
between these two blocks vanish by reflection symmetry. The equal-weight
cosine sum and sine difference are orthogonal unit trial directions; their
common covariance Rayleigh quotient is

`(lambda4+lambda5+2 sqrt(lambda4 lambda5)x)/24 > 313/960`.

The strict interval enclosures prove lambda4>2.19, lambda5>2.29 and
sqrt(lambda4 lambda5)>2.23. Also `exp(2)>7` already follows from six
positive Taylor terms, proving x>3/4 without floating arithmetic.
For every supplied g>=3.7, `1/g<=10/37<313/960`. Consequently
`I7/g-Cov_p(X)` is negative definite on this two-dimensional subspace.
It has **at least two negative eigenvalues**. Numerically at the reported
coexistence gain the two eigenvalues are about -0.06125315.

This is an exact counterexample to an **everywhere** full-seven-coordinate
curvature/index bound. The chosen field is not claimed stationary. It does
not refute a separate statement restricted to Morse indices of stationary
points, nor the still-open positive-orthant four-amplitude ceiling. A
stationary-point-only theorem would require an additional proof.

## 8. Spectral/discord comparison without causal promotion

For C_density=I/12+W/20 the exact uniform-to-alternating gap is lambda6/20;
the top neighboring spectral gap is (lambda6-lambda5)/20. Their product
divided by 7392 has the handoff's spectral expression. But the existing
discord theorem uses rational **lower enclosures** Delta_L and delta_L,
not identities equating rounded lower bounds to transcendental values.
Keep its certified d0 unchanged. The growth splitting
`g(lambda6-lambda5)/12` assumes the supplied linear law
`q_dot=(-I+g A/12)q`. It is not derived from the energy alone, and no
causal equivalence between discord and localization follows.
