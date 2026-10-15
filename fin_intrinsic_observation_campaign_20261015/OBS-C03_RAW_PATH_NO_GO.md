# OBS-C03 — Non-tight critical-time raw snapshot path

**Status:** scoped asymptotic non-tightness proof for k_N~alpha sqrt N, lambda bounded, t_N selected within accepted local window. Not a denial of finite-dimensional limits or a statement about fixed-k integrated records.

Let v_N=m_B−m_N, k_N~alpha sqrt N. The exact labelled generator yields

 L v_N=−v_N+e_N, |e_N|<=C/N,
 Γ(v_N)=sum_a c_a² 2(1−S_a q_{xi_a}),
 c_a=1_{a∈B}/k_N−1/N,
 q_i=sum_j sigma_j softmax_j[(g_N/N)A(n−e_i)].

Under local-critical stationary preparation, p_N->u on every fixed O(1) model-time interval, hence q_i->0 uniformly in i. The private fluctuation U_N(t)=sqrt(k_N)*v_N(t_N+t) obeys the martingale integral representation

 U_N(t)=U_N(0)−∫_0^t U_N(r)dr+M_U,N(t)+o_P(1)

uniformly on each fixed interval. Jumps are ≤2/sqrt(k_N); predictable bracket of M_U tends to 2t; its cross variation with the global collective mode over this model-time interval vanishes. The hypergeometric CLT at t_N gives U_N(0)=>G independent of Z(t_N) with G~N(0,1). A stopped martingale FCLT therefore yields the **stationary OU** dU=−U dt+sqrt(2)dW in model-time t, conditionally asymptotically independent of the slowly changing Z. This uses an actual labelled refresh martingale, not imposed private white noise.

The unfiltered critical-time snapshot readout R_N(s)=N^{1/4}m_B(t_N+sqrt N s) satisfies

 R_N(s)=Z_N(s)+(N^{1/4}/sqrt(k_N))U_N(sqrt N s).

Choose h_N=1/sqrt N ->0. At any fixed s in (0,S), over [s,s+h_N] the slowly varying Z changes by o_P(1), whereas the OU private component changes by a nondegenerate law with variance

 lim Var(R_N(s+h_N)−R_N(s))=2(1−e^{-1})/alpha.

The difference is Gaussian in the limiting bivariate OU law, so its tail probability beyond a fixed epsilon stays strictly positive. This violates the deterministic-time stochastic continuity necessary at almost every fixed time for tightness in D([0,S],J1). More explicitly: if a subsequence were J1-tight, any cadlag limit has at most countably many fixed discontinuity times of positive probability; choose s outside that set, where J1 convergence entails vanishing differences over deterministic s+h_N. Contradiction. This is an **ordinary cadlag J1 no-go**, not a statement that the finite microscopic path has no definition.

An instantaneous Gaussian channel is valid at each time but cannot be promoted to a continuous/readout SDE for its raw critical-speed time series. Applying a filter or integral changes the observation ID and cannot be represented as the same raw output. No physical detector or irreversible measurement law is identified.
