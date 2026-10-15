# OBS-C04 — Joint microscopic state / integrated tagged-history record

**Scientific status:** analytic proof submission for the primary labelled leave-one-out heat-bath, under the accepted MT local equilibrium window and UF stopped/unstopped path convergence. The result is a classical *mathematical path functional*, not a device or independent physical observation law.

## Exact finite-generator identity

For a fixed, predetermined block B of size k (k does not depend on N), write S_a=(-1)^{xi_a}, m_B=k^{-1}sum_{a∈B}S_a, m_N=N^{-1}sum_a S_a, v_N=m_B−m_N. For every departure a, define q_i^{(N)}=sum_j sigma_j softmax_j[(g_N/N)A(n−e_i)]. Let \bar q_B=k^{-1}sum_{a∈B}q_{xi_a}, \bar q_N=N^{-1}sum_a q_{xi_a}. Then **exactly**

 L_N v_N=−v_N+e_N,  e_N=\bar q_B−\bar q_N.

Softmax is Lipschitz on the fixed finite q kernel and bounded gain region, hence sup_{i,j}|q_i−q_j|<=C/N, so |e_N|<=C/N uniformly over ALL configurations. No exact finite-N parity lumpability was assumed.

Let M_v be the compensated martingale: dv_N=(−v_N+e_N)dt+dM_v. Integration from t_N to t_N+sqrt(N)s yields, pathwise,

 ∫v_N(t)dt = M_v(t_N+sqrt(N)s)−M_v(t_N)−[v_N(t_N+sqrt(N)s)−v_N(t_N)]+∫e_N(t)dt.

After multiplying by N^{-1/4}, the endpoint term is bounded by 4N^{-1/4} and the error integral by C S N^{-3/4}. Since N^{-1/4}∫m_N dt=∫_0^s Z_N(r)dr exactly, the record

 Y_N(s)=N^{-1/4}∫_{t_N}^{t_N+sqrt(N)s}m_B(t)dt

is ∫_0^s Z_N(r)dr + M_Y,N(s)+o(1) uniformly on [0,S], where M_Y,N=N^{-1/4}(M_v(t_N+sqrt(N)·)−M_v(t_N)).

## Exact predictable brackets

For a refresh of slot a from old label i to j, set c_a=1_{a∈B}/k−1/N. Then Δm_N=(sigma_j−sigma_i)/N, Δv_N=c_a(sigma_j−sigma_i). Conditional signed refresh q_i gives h_a=2(1−S_a q_{xi_a}) and

 Γ(m_N)=N^{-2}sum_a h_a,
 Γ(v_N)=sum_a c_a² h_a,
 Γ(m_N,v_N)=N^{-1}sum_a c_a h_a.

On the MT-derived locally critical occupation path, sup_{0<=s<=S}||p_N(t_N+sqrt N s)−u|| ->0 in probability: from the UF compact/unstopped theorem with MT initial moments and parent/reflect coupling. Therefore sup_i|q_i|->0, and uniformly in scaled time

 N Γ(m_N) -> 2,  Γ(v_N)->2/k,  sqrt N Γ(m_N,v_N)->0.

In the joint rescaled martingales M_Z,N=N^{1/4}M_m(t_N+sqrt N·) and M_Y,N the brackets are thus <M_Z,N> ->2s, <M_Y,N> ->(2/k)s, <M_Z,N,M_Y,N>->0. The cross term even has the crude bound O(N^{-1/2}) without cancellation from h_a≈2, because sum_a |c_a|<=2. Rescaled jump magnitudes are O(N^{-3/4}) for Z and O(N^{-1/4}) for Y, both vanish. All brackets are uniformly bounded on each compact scaled horizon (Γ(m)<=4/N; Γ(v)<=4(1/k−1/N)). BDG/Doob yields L² stochastic boundedness; bounded jumps plus bracket convergence license the joint vector martingale FCLT.

## Initial/limit and independence

At t_N the exchangeable lift of π_D,N applies with TV error epsilon_N(t_N). MT-003 gives the quartic Z initial law and preparation controls, so UF path convergence applies to the resulting stationary Z process. The joint martingale problem passes using the natural filtration, joint brackets and square-martingale UI; the SDE admits pathwise uniqueness (cubic dissipative globally) and the two-dimensional martingale bracket is diagonal, yielding independent Brownian *drivers*, not independent state and record. The limit is

 dZ=(chi lambda Z−Z³/3)ds+sqrt(2)dB_s,   Z(0)~rho_lambda,
 dY=Z ds+sqrt(2/k)dW_s,                Y(0)=0,
 with B and W independent; Z(0) is independent of their **future increments**, not necessarily of an arbitrary pre-t_N past.

For N-dependent observation time, this is an unconditional iid-start result because the TV path comparison at t_N and boundary exit estimates are included. No mixing exponent N0 or finite-N quantitative error beyond inherited asymptotic o(1) is asserted.

MET/BAR are not covered by the uniform HB contrast identity; their tagged private relaxation and corrector/Poisson equation remain OPEN_CONTROL. A chosen record of a tagged stochastic process is not a physical stored record, hardware readout or measurement-backaction model.
