# MT root-path, nonuniform congestion and lattice-tail audit

Date: 2026-10-15.
Status: analytical intake completion under the admitted exact spectral inputs.
This is not an all-type numerical enumeration or a formal theorem-prover certificate.

## 1. Fixed geometry and the outside-central gradient

The [rational checker](check_geometry.py) pays the displayed epsilon=.01,
delta=.0001 geometry, Hessian floor `8679/98750>0.08`, entropy constant
`C_m=1/(24p_min^2)<6.14` and fast-gradient margin `0.08-C_m delta>0.0793`.
The gain/spectral premises are bound to the previously accepted rational boxes.
The source's unspecified N0 still has to absorb the explicit O(1/N) and
bounded-lambda/N perturbations; no practical N0 is certified here.

Use the source central scales

    a_N=N^(-1/4)(log N)^(1/4),
    d_N=N^(-1/2)(log N)^(1/2).

The exact parity gradient and positive fast Hessian imply

    ||grad_T F|| >= c(|m|^3+||w||)-C Lambda N^(-1/2)|m|.

Outside the central box there are two cases. If |m|>=a_N, the adverse gain
term divided by |m|^3 is at most C Lambda/sqrt(log N), which vanishes. If
|m|<a_N but ||w||>=d_N, that term divided by ||w|| is at most
`C Lambda N^(-1/4)(log N)^(-1/4)`, also vanishing. In both cases the tangent
gradient norm is bounded below by a positive multiple of
`N^(-3/4)(log N)^(3/4)` for sufficiently large N.

Subtracting the mean from an ambient gradient changes no max/min difference.
For a mean-zero q-vector the max/min range controls its Euclidean norm by a
fixed q-dependent factor. Choosing a maximum departure gradient and minimum
arrival gradient therefore gives an actual count-step first-order decrement
at least `c N^(-7/4)(log N)^(3/4)`. The full Hessian norm is uniformly bounded
on a fixed positive-interior enlargement of the cylinder; its quadratic step
error is at most C/N^2. The decrement wins uniformly for all large N.

Every departure count is positive there. The sublevel's fixed buffer from the
cylinder boundary ensures the step cannot leave the cylinder. It decreases F
and remains in the exact sublevel. The finite count graph prohibits infinitely
many strictly decreasing distinct states, so the descent reaches the central
box. This is a constructive all-type argument, not evidence inferred from a
few gradient probes or from the continuous path alone.

## 2. Central count balancing and the energy ceiling

Within each parity class, move one count from a maximal to a minimal entry
until all six entries differ by at most1. Each such move decreases the
within-class squared deviation from its fixed parity mean. It holds m fixed,
so the fast norm never exceeds its initial value apart from O(1/N) rounding.

Next move one count from the overpopulated parity class to the other, choosing
entries/rebalancing so each class remains balanced within O(1) counts. The
absolute parity mass decreases by2/N on every cross-parity move. N divisible
by12 ensures the final m=0 type is exactly `(N/12,...,N/12)`.

Along this route, |m|<=a_N+O(1/N) and ||w||<=C(d_N+1/N). Exact vanishing of
the fast first derivative on the parity line, bounded Hessian and the parity
Taylor series give `|F|<=C(log N+1)/N`. The route is therefore inside the fixed
sublevel and cylinder for all sufficiently large N. A legal path exists for
every central type; every fixed compact box in scaled (Z,Y) is included too.

Concatenating descent with balancing gives the claimed root path. Coercivity
gives `F(start)>=-C_Lambda/N`, so its maximum is at most

    F(start)+C_*(1+log N)/N.

Erase loops to choose simple paths, preserving the energy ceiling. Their
length is bounded by the number of local types. Lexicographic tie-breaking
is a proof device, not a physical orientation selector. Full local sublevel
connectivity and iid-start coverage follow in the stated N-subsequence.

## 3. Nonuniform weighted canonical paths

The cited [Olesker-Taylor–Zanetti Proposition4.7](https://link.springer.com/article/10.1007/s00440-023-01257-x#FPar62)
does address continuous-time paths; its particular weighted-walk convention
has **uniform** stationary measure. It is a related methodological reference,
not a direct license to treat FIN Gibbs weights as uniform.

For FIN use its own identity

    Var_pi f=(1/2) sum_(x,y) pi(x)pi(y)(f(x)-f(y))^2

and the actual conductance `Q(e)=pi(z)q(z,z')`. Telescoping along the
root-concatenated path and Cauchy–Schwarz bounds the variance by weighted
edge congestion times the actual reversible Dirichlet form. No external
uniform-measure formula is substituted.

Type bounds and the root energy ceiling give, at **both endpoints** z of every
edge on P_x,

    pi(x)<=C(N+1)^(12+ceil(C_*)) pi(z).

For a fixed edge, pairs traversing it are a subset of pairs whose x root path
or y root path traverses it. Sum the other endpoint's probability to1; then
bound the remaining sum by `#D_N*C poly(N)*pi(z)`. Path lengths contribute
another `#D_N` factor. Actual legal rates are bounded below by a positive
N-independent constant on the positive-interior domain, so division by Q(e)
does not insert an exponential normalization. With `#D_N<=(N+1)^12`, congestion
is polynomial. The source K=24+12+ceil(C_*) is a valid conservative exponent;
constant-factor conventions in the Dirichlet form can be absorbed in c_*.

Consequently the polynomial lower gap and the source L2/TV burn-in argument
survive the all-type path audit. This is not an optimized N^(-1/2) lower gap.

## 4. Scaled lattice volume and moment domination

The normalized type lattice has tangent covolume `sqrt(12) N^(-11)` in
Euclidean surface measure. In an orthonormal basis `(sigma/sqrt(12),E_fast)`,
the map `(p-u)->(Z=N^(1/4)m,Y=sqrt(N)E_fast^T w)` has determinant
`sqrt(12) N^(21/4)`. Thus a scaled fundamental cell has volume

    v_N=12 N^(-23/4).

Every cell has Z width O(N^(-3/4)) and fast diameter O(N^(-1/2)). The cells
need not be axis-aligned. On bounded scaled boxes, uniform interior Stirling,
the source Taylor expansion and full-sublevel inclusion give ordinary
Riemann convergence with this common volume factor.

For tails, use the source bound

    N F >= Z^4/24+0.04||Y||^2-C_Lambda.

A global bound on the **scaled gradient** is unnecessary and is not generally
uniform on the expanding domain. Instead bound oscillation of this envelope
within each scaled lattice cell. Since |Z|<=epsilon N^(1/4) and
||Y||<=delta sqrt(N), products `|Z|^3*Delta Z` and `||Y||*||Delta Y||` are
bounded by fixed constants. Higher cell increments are bounded as well.
The lattice tail sum is therefore bounded by a fixed multiple of the
neighbouring continuum quartic/Gaussian tail integral. Stirling prefactors
are uniformly bounded in the positive cylinder; normalization has a positive
lower bound from any fixed scaled box. The same argument with polynomial
weights proves every fixed polynomial moment's uniform integrability.

For completeness, the local unnormalized partition satisfies

    W_N(D_N)/N^(1/4)
      -> 12^5 (2pi)^(-11/2)
          integral exp(chi lambda z^2/2-z^4/12-y^T H_fast y/2) dzdy.

Only the convergence/positivity and moment statements, not extra digits of
this prefactor, are needed for MT's mixing argument.

The source MT-003 sentence saying the boundary has weight **"at least"**
exp(-Nb) has the wrong inequality direction. The valid statement is an
**upper** bound, at most a polynomial times exp(-Nb), for the relevant
high-energy boundary layer (and no allowed sublevel types on the excluded
outer cylinder). The proven coercivity/type bounds give that upper bound;
the original wording is preserved, not silently accepted.

## 5. What this completion does not provide

No certified numerical C_*, K, N0 or gap prefactor is exported. The theorem
is asymptotic for fixed bounded lambda and N through sufficiently large
multiples of12. No small-N autonomy, physical seconds, global escape prefactor,
space, particle or fundamental information-conservation claim follows.
