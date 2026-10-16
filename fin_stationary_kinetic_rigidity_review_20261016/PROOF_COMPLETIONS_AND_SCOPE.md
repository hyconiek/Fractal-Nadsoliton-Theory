# KR proof audit: stationary class convergence and its exact open boundary

Date:2026-10-16. Analytical model proof under the inherited exact spectral,
MT local-domain and moment premises. Not mechanized or peer reviewed.

## 1. Pay the dependencies now, not historical access

The actual AC, rate, UF, MT and OBS artifacts are available locally and bound
in verification.json. The source's three extracted UF/MT members and two
archive wrappers match actual bytes, including membership checks inside the
verified ZIPs. Its historical AGENTS snapshot remains unavailable.
This supplies the present proof chain, not a retroactive remote checkout or
evidence that the producer read the current reviews or original research order.

AC already proves the stationary covariance identity and examples. KR's new
contribution is the full stationary class path theorem under NC, plus a
weaker static test identity without NC. Do not count old wheel, Gibbs,
scalar drift or critical coefficients as newly derived laws.

## 2. The class is explicit and nonempty

The generator is a count CTMC on the actual reflected MT domain D_N with
the same Gibbs count law pi_D,N, including multinomial multiplicity. It
has only single replacements, total rate<=CN, exact **reflected parity
drift** Ltilde m=L0 m at every state and exactly stationary initialization.
N runs through multiples of12, g_N=g_c+lambda/sqrt(N), with lambda fixed
in a bounded interval. No pointwise noise, fast gap or Markov scalar quotient
is assumed. The source does not impose a positive lower rate on every edge.

Same-parity thinning is a genuine non-HB example: altered moves have
delta m=0, so both scalar drift and Gamma(m,m) are unchanged at finite N.
It preserves detailed balance and the native proposal capacity. This does
not imply equality of its full finite-N parity histories, since the changed
other coordinates affect later rates. The entire KR class has only its
declared CN cap; do not grant every member unit per-slot attempt capacity.

## 3. The stationary drift error is fully paid

Use the actual UF-002 scalar function, not a conditional expectation or
only its value on the parity line. At p=u+m sigma/12+w the source lists
the correct six terms in

    e_N=N^(3/4)Ltilde m-b(N^(1/4)m),
    b(z)=chi lambda z-z^3/3.

MT supplies all fixed moments of N^(1/4)m and sqrt(N)w. Holder therefore
pays E[(1+|Z_N|)|e_N|]->0 with the listed powers. No stochastic fast-mode
energy inequality for the new generator is needed for this stationary proof.

For completeness, a fixed smaller UF expansion neighborhood can be used
without assuming it contains every MT domain state. Outside that neighborhood,
the global rate/jump cap gives |N^(3/4)Ltilde m|<=C N^(3/4), and |Z_N|<=N^(1/4).
Thus (1+|Z_N|)|e_N|<=C N. Markov's inequality with a sufficiently high MT
moment gives outside probability o(N^-1). This pays the omitted-tail part.
At the native departure layer the scaled weighted correction is at most
C N pi_D(B_D), which vanishes by the actual exponential upper bound.

Exact stationarity converts these one-time payments into integrals over
the entire fixed scaled horizon S, even though model time is sqrt(N)S.
No mixing rate or independence across times is inserted.

## 4. NC is sufficient for the actual bracket, not derived from its mean

Invariance of m^2 gives E Gamma(m,m)=-2E[m Lm]. Hence exact drift matching
and common pi_D give E G_N=E[N Gamma0(m,m)]->2, with G_N=N Gammatilde(m,m).
The cap gives 0<=G_N<=4C. Native interior continuity at u and the paid
boundary mass establish the limit2; it is not inferred just from a rounded
g_c or auxiliary Gaussian covariance.

NC is E|G_N-EG_N|->0. It implies

    E sup_(t<=S)|integral_0^t G_N(X_(sqrt(N)s))ds-2t|
       <=S E|G_N-EG_N|+S|EG_N-2| ->0.

The power is correct: accelerated time sqrt(N) times squared amplitude
N^(1/2) times original Gamma gives N Gamma, not sqrt(N) Gamma.

## 5. Maximum containment, UI and the joint filtration

The true scalar martingale decomposition is

    Z_N=Z_N(0)+integral b(Z_N)ds+integral e_N ds+M_N,
    <M_N>_t=integral_0^t G_N ds.

The source's maximal estimate is valid. Jump Ito for Z_N^2 includes
Gamma as the positive quadratic compensation. The term2z b(z) is bounded
above; E integral |Z_N e_N| vanishes. For its stochastic integral,

    E sup|integral Z_N- dM_N|
       <=C [E integral Z_N^2 G_N ds]^(1/2)
       <=C [4CS E_pi Z_N^2]^(1/2).

Thus stationary MT moments imply a uniform second moment of the maximum,
not merely deterministic-time tightness. Localized drift increments and
the vanishing total residual then pay Aldous tightness.

The compressed fourth-moment step can be written explicitly. The bracket
is bounded by4CS and |Delta M_N|<=2N^(-3/4). Its fourth-jump compensator
is at most4N^(-3/2)<M_N>_S. The jump BDG/fourth-power identity gives
uniform E sup|M_N|^4. This pays UI of the martingale and its square-minus-
bracket; a bounded expected bracket alone would not suffice for every
claimed passage. Vanishing jumps give continuous limits.

Pass both martingale identities jointly with the initial state and the
count-path past. The limit has bracket2t in its joint natural filtration;
Levy characterization gives sqrt(2) Brownian noise with future increments
independent of the initial quartic variable. The inward cubic drift has a
unique nonexplosive law. Consequently admit stationary path convergence

    dZ=(chi lambda Z-Z^3/3)ds+sqrt(2)dB, Z0~rho_lambda,

on every fixed scaled horizon, under NC and the exact frozen class.
No iid preparation, diagnostic local mixing, labelled record law or
unbounded-time/operator statement is supplied by this proof.

## 6. The additional static identity without NC is correct

For h smooth compactly supported, take f''=h. One can choose bounded f',
bounded f''' and f with at most linear growth. Each finite count space
makes the stationary generator test legitimate. Taylor expansion gives

    sqrt(N)Ltilde f(Z_N)
       =f'(Z_N)[b(Z_N)+e_N]+h(Z_N)G_N/2+R_N,
    |R_N|<=C_f N^(-3/4).

The moment/drift payments pass its expectation. Since rho_lambda satisfies
integral[b f'+f'']rho=0, including the zero boundary term from its quartic
tails, one obtains

    E[h(Z_N)(G_N-2)] ->0.

Any joint subsequential law therefore has E[G|Z]=2. It need not have G=2,
and it supplies no trajectory-wise averaging. This is not a posterior or
an instantaneous Markov closure theorem for Z_N.

The proof under NC extends immediately if NC is replaced by the source's
IBA condition, E sup|integral(G_N-EG_N)|->0. All other estimates already
hold independently of NC. This is a sufficient conditional theorem, not
a proof that the class assumptions alone imply IBA.

## 7. Keep the negative controls in their real scope

The exact producer calculation is q4/N2 with a uniform multinomial law.
It is not q12/A7 or a critical-N sequence. Its embedded q12 construction
uses a **six-edge scalar-parity stress**, not the earlier D12-averaged
six-spoke stress preserving every coordinate drift. It need not itself be
D12-equivariant; the KR class does not require that extra symmetry.

The embedding's index3 is exactly the uniform type. Its unnormalized
squared parity-jump coefficient is-8 there and+8 at index4. Positivity,
CN cap and same reflected drift follow from perturbing interior reversible
conductances by a sufficiently small fraction of their minimum. Its support
has ten types whose joint scaled positions approach zero, so its stationary
mass vanishes. The spike satisfies NC while violating pointwise equicontinuity.
It is NOT a different-limit counterexample.

The intake independently checks a nonuniform q3/N3 Gibbs example, not just
the producer's uniform fixture. Finite checks corroborate algebra and
normalizations, not the stochastic theorem. The stronger NC-free IBA claim
remains genuinely OPEN. The new reversible Poisson criterion is a precise
next proof route, not its silent completion.
