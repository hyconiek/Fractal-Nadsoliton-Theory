# MT-004 — A conservative *polynomial* local spectral-gap payment

**Scientific result:** a conditional analytical proof of a finite-N local relaxation window by energy-controlled paths. The exponent is extremely loose and contains an unspecified, *fixed* analytic Taylor constant. The statement uses only actual connected type edges, original parent internal rates and the exact Gibbs conditional stationary law; it does **not** assume finite-N parity lumpability.

## Reflected chain and exact coupling

On the finite connected domain D_N of MT-001 define the auxiliary reflected generator by retaining every parent off-diagonal rate L_N(x,y) whenever x,y∈D_N and setting outgoing rates to zero, with diagonal equal to minus the retained outgoing sum. Detailed balance is inherited, so π_D,N=W_N(.)/W_N(D_N) is invariant and reversible. Coupling reflected and unchanged parent using the *same* internal jump clocks makes their paths coincide **until the actual parent first exit**. The reflection is not an added microscopic law or a physical barrier.

Let T_N=#types(D_N)≤(N+1)^12. For sufficiently large N all labels have positive counts in D_N. A legal directed count move x→y has rate ≥r_*>0, uniformly in N and in HB/MET/BAR, because all edge acceptance probabilities have positive lower bounds on the compact interior with g_N bounded. In particular HB r_ij=n_i softmax_j(h_i), MET r_ij=(n_i/12)min(1,exp Δ_ij), BAR r_ij=(2n_i/12)/(1+exp(−Δ_ij)); each n_i≥1. Let E_*,Λ be fixed.

## Canonical-path construction (the nontrivial gap step)

MT-001 supplies, for **every** x∈D_N, an actual legal path P_x from x to the balanced seed n_* of length at most 2T_N and energy ceiling

 max_{z∈P_x} F_{g_N}(z/N) ≤ F_{g_N}(x/N)+C_*(1+log N)/N.  (1)

The proof has two regimes:

* Outside the central log-scale box, take a departure label with maximum ∂_i F and arrival label with minimum ∂_j F. The gradient norm is bounded below by c N^-3/4(log N)^3/4; a single exact type move decreases F because its gain O(N^-7/4 log^3/4) dominates its smooth second-order O(N^-2) remainder. No boundary is crossed, by MT-001's fixed positive cylinder buffer.
* Inside the box, balance six labels within each parity sector using legal same-parity transfers, then move parity mass toward zero with rebalancing. m stays within O(a_N), the fast displacement within O(b_N), and F stays ≤C log N/N. The lower bound F(x)≥−C_Λ/N pays the energy comparison (1).

These instructions give an explicit deterministic rule after choosing a fixed lexicographic tie-break. There is no assumption of an abstract low-dimensional rate or a reflective physical controller.

## Weight comparison, congestion and spectral gap

The exact elementary type bounds, valid also for boundary types (but used here in the interior), give

 (N+1)^(-12)exp[−N D(p||u)] ≤ μ0_N(n) ≤exp[−N D(p||u)].

Therefore for every z on P_x,

 π_D(x)/π_D(z)=W_N(x)/W_N(z)
 ≤ (N+1)^12 exp[N(F(z/N)−F(x/N))]
 ≤ e^{C_*}(N+1)^{12+C_*}.  (2)

Choose finite integer K0=12+ceil(C_*), increase prefactor if necessary, and write a uniform bound π_D(x)≤C (N+1)^K0 π_D(z) for every edge endpoint z on P_x. For each ordered pair x,y define γ_xy as P_x followed by the reverse of P_y with loops erased. Then |γ_xy|≤4T_N; every edge of γ_xy lies on at least one of P_x or P_y.

For an undirected edge e={z,z'}, the sum of π_D(x) over all x whose root path traverses e is at most T_N·C(N+1)^K0 π_D(z), by (2) applied to the actual endpoint on each such path. Because all rate rules are reversible, its Dirichlet conductance is Q(e)=π_D(z)L_D(z,z')≥r_* π_D(z). It follows that the standard continuous-time canonical-path congestion satisfies

 ρ_N = max_e [1/Q(e)]Σ_{x,y:e∈γ_xy}π_D(x)π_D(y)|γ_xy|
 ≤ (8C/r_*) (N+1)^{24+K0}.

The one-line underlying Poincaré proof is: for each (x,y), telescope f(x)−f(y) along γ_xy and use Cauchy–Schwarz with |γ_xy|; sum against π(x)π(y), reverse the finite sums, then compare each edge contribution with the Dirichlet form. Thus

 gap(L_D) ≥ c_*(N+1)^(-K),     K=24+K0<∞. (3)

K may be huge: no claim that it matches the critical N^-1/2 rate. Importantly it is POLYNOMIAL because the root path's excess energy is only O(log N/N), not a fixed O(1) barrier. The construction is specific to this small local domain and cannot be used for the full global Gibbs chain across distinct ordered basins.

## Independent upper-gap obstruction

Rayleigh with parity observable f(n)=m(n) and the MT-003 positive limiting variance gives

 gap(L_D) ≤ E_π[f,−L_Df] / Var_π(f)
 ≤ (4/N) / (c/√N)=C N^-1/2,

since every parent jump changes m by at most 2/N and total genuine off-diagonal rate is ≤2N. Thus *ten stable occupation directions* do not give a constant global reflected spectral gap. Equations (3) and this upper bound are compatible but far apart.

## Pay the actual iid initial density

Let ν_N=μ0_N(.|D_N). Pointwise W_N≥μ0_N, so

 dν_N/dπ_D(n)= [W_N(D_N)/μ0_N(D_N)] exp[−(g_N N/2)p^TAp] ≤W_N(D_N)/μ0_N(D_N).

MT-001 gives F≥−C_Λ/N on D_N; the type bound and ≤(N+1)^12 states therefore imply W_N(D_N)≤e^{C_Λ}(N+1)^12. Since μ0_N(D_N)→1,

 ||dν_N/dπ_D−1||_L2(π_D) ≤C(N+1)^6.

For reversible continuous-time semigroups, L2 contraction by gap and Cauchy–Schwarz yield

 TV(ν_N exp[tL_D],π_D) ≤ C(N+1)^6 exp[−c_* t/(N+1)^K]. (4)

A conservative local relaxation time is, for example, t_loc(N)=(N+1)^(K+2); the bound in (4) then tends to zero. This is an existential polynomial guarantee, not a practical/runtime estimate or a finite-N calibrated constant.

## Transfer to the unchanged parent, uniform long window

Initialize a coupled reflected chain in ν_N when the iid parent starts in D_N. The exceptional mass μ0_N(D_N^c)=o(1) is retained without repair. Until τ_D, parent and reflected paths agree, so for all t≥0,

 TV(Law_iid[p_N(t)],π_D,N)≤ μ0_N(D_N^c)+TV(ν_N exp[tL_D],π_D,N)+P_iid(τ_D≤t).

Use MT-002: P_iid(τ_D≤t)≤μ0_N(D_N^c)+2Nt(N+1)^12 e^(−Nb). Thus for every fixed 0<c<b,

 sup_{t_loc(N)≤t≤exp(cN)} TV(Law_iid[p_N(t)],π_D,N) → 0. (5)

Equation (5) is **the desired finite-N local equilibrium/trapping window**, conditional on the exact inherited source spectrum and the analytic local path inequalities in MT-001. The bound is intentionally conservative and does not imply a feasible laboratory N, the actual physical clock, or a matching exit-time asymptotic. The parent is not reflected and not postselected; the error includes genuine outside/exit trajectories.

## Comparison with external methods

The last inequality is the continuous-time, finite-reversible-chain canonical-path Dirichlet argument; see e.g. Proposition 4.7 of *Geometric bounds on the fastest mixing Markov chain* (Probability Theory and Related Fields, 2024), DOI/URL https://link.springer.com/article/10.1007/s00440-023-01257-x . We prove the needed telescoping/Cauchy–Schwarz comparison explicitly above; no theorem about globally exponentially stable nonlinear drift, a selected phase or arbitrary initial time is imported.

**Remaining gap if independent audit rejects (1):** a rigorously exact *all-local-type gradient/inner-corridor* certificate, not a new kinetic law, source, fitted pole or simulation. In that case (3)–(5) remain CONDITIONAL and MT-001/002/003 should be kept as separate independent results.
