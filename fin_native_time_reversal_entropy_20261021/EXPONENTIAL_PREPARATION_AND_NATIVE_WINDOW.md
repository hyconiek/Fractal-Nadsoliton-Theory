# Paying the actual iid exceptions and full-parent likelihood

Date: 2026-10-21. This is the new native application. Inherit the accepted
exact spectral, MT domain and polynomial reflected-gap premises.

## 1. An exponential iid coverage bound, not just the old o(1)

Set q=12 and r=1/10^6. The MT cylinder uses epsilon=1/100, delta=1/10000
and h=1/10^10. If ||p-u||_2<=r, then |m|<=sqrt(12)r<4r<epsilon,
||w||_2<=r<delta and every p_i>1/24. Along the segment from u to p the
entropy Hessian is <=24 I. Since g_N>=0 and A7 is PSD,

    F_g(p)<=D(p||u)<=12||p-u||_2^2<=12r^2<h.               (1)

This ball is inside the ACTUAL component D_N, not only its continuous
sublevel: move a count from an above-uniform label to a below-uniform
label. N in 12N makes the uniform count integer. Each step strictly
decreases ||p-u||_2^2, stays in the same ball/cylinder/sublevel, and reaches
the uniform type. No artificial enlargement or unproved component choice.

Each empirical label frequency has Bernoulli mean 1/q. The elementary
bound log E exp[t(B-EB)]<=t^2/8 follows by integrating the tilted variance
bound <=1/4 twice. Chernoff and a union bound therefore give

    delta_0,N=mu_0(D_N^c)
       <=2q exp[-2N r^2/q] =24 exp[-kappa N],
    kappa=1/(6*10^12)>0.                                 (2)

This strengthens the old central-box Markov estimate. It introduces no
new law or conditioning. The constants are deliberately tiny and provide
no practical finite-N preparation threshold.

## 2. Actual labelled future-path TV at the selected origin

The paid MT gap and initial density estimate give, for the reflected count
chain started in nu_N=mu_0(.|D_N),

    TV(nu_N exp[tL_D],pi_D,N)
      <=C(N+1)^6 exp[-c_* t/(N+1)^K].                     (3)

At tau_N=(N+1)^(K+2), call this bound eta_N. The initial conditioned iid
law, reflection and generator are slot-permutation covariant. Their
labelled laws are uniform on each count orbit. The exact OBS orbit-TV
isometry lifts (3) to the entire labelled state law, without assuming a
separate constant labelled spectral gap.

Couple the full iid parent to that reflected process until first departure
and couple the latter at tau_N to stationary reflected Gibbs. Initial
outside probability is retained. MT's actual weight-domination exit bound
through the COMPLETE interval [0,tau_N+T] gives

    TV(P_N,T,S_D,N,T)<=epsilon_N(T),
    epsilon_N(T)=2 delta_0,N+eta_N
                  +2N(tau_N+T)(N+1)^12 exp[-bN],          (4)

where S_D,N,T is stationary reflected LABELLED path law, b=1/(2*10^10).
The conditional law is only an analysis comparison. P_N,T is always the
unconditioned, unreflected actual parent. Every future record length is
paid through tau_N+T, not just through its observation origin. Markov
contraction carries state-law coupling to full future paths.

## 3. A global finite likelihood bound including rare outside states

The actual labelled chain has total genuine exit rate <=N. Starting iid
uniform, the event of no genuine jump at x before tau_N gives

    mu_N(x)>=q^(-N) exp[-N tau_N]>0                       (5)

for EVERY configuration, including those outside D_N. No attempt marks
are observed or needed for this bound. Let |g_N|<=G and M_A=max|A7_ij|.
The full Gibbs potential has absolute value <=G M_A N/2, so

    Pi_N(x)>=q^(-N) exp[-G M_A N].                        (6)

Using mu_N<=1 and Pi_N<=1, (5)–(6) imply the conservative bound

    osc log(mu_N/Pi_N)
       <=B_N=N[tau_N+2 log q+G M_A].                     (7)

It is polynomial in N, despite very small individual probabilities.
The same conservative bound works on counts, including their multinomial
multiplicity. It is NOT a claim that mu_N is globally Gibbs-close.

Indeed D(P_N,T||S_D,N,T)=infinity: the actual positive mu_N assigns mass
to outside configurations whereas S_D starts only inside D_N. The proof
never invokes such a singular KL comparison. Both P_N,T and its own
reversal have common support, and their endpoint likelihood is bounded
by (7).

There is a further exact simplification for this parent. Both mu_N and
Pi_N are uniform within each labelled count orbit, so mu_N(x)/Pi_N(x)
equals the corresponding count-law ratio. The entire labelled path
reversal likelihood is therefore measurable in the two count endpoints.
Their joint law retains exactly the full reversal KL. This is not a
claim that count endpoints determine general tagged histories or laws.

## 4. Uniform reversal-entropy and model-production bound

Apply the finite path lemma with the FULL reversible parent Pi_N and the
stationary reflected comparison only in TV. Equations (4) and (7) give

    0<=E_N(T)<=A_N(T)<=2 B_N epsilon_N(T).                 (8)

Fix any 0<a<b. Uniformly for 0<=T<=exp(aN), the RHS is at most a fixed
polynomial times exp(-kappa N)+exp[-(b-a)N], plus a polynomial times
exp[-c_*(N+1)^2]. Thus

    sup_(0<=T<=exp(aN)) A_N(T)->0,
    sup_(0<=T<=exp(aN)) E_N(T)->0.                        (9)

This includes the actual critical records T=sqrt(N)S at every fixed S.
Equivariant count, parity and integrated records inherit (9) by data
processing, with their actual reversal maps. No weak-limit entropy
continuity, postselection, global mixing theorem or new noise premise.

## 5. Exact limits of the conclusion

The origin is the paid polynomial tau_N; arbitrary origins covered only
by an unspecified o(1) TV estimate do not receive an entropy upgrade.
Initial preparation can be irreversible. Global Gibbs relative entropy
can remain large while its drop over this locally trapped future is small.
The exponential horizon is an existence statement with unknown K,c_*,N0,
not feasible lifetime guidance or a sharp escape result.

E_N is the Markov path entropy-production functional in nats, not derived
physical heat, temperature, units or a universal cosmological arrow. The
theorem says that these unpostselected near-local-equilibrium passive
histories do not retain a nonvanishing reversal-KL arrow in this regime.
It does not rule out other preparations, active records or new supplied
physical structures, and it does not close F01 or observer ontology.
