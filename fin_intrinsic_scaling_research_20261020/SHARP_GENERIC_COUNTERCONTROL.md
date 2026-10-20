# IS-003: bounded local ratios cannot force a constant activity bound

Date: 2026-10-20. An exact control family, NOT the actual A7 Gibbs source
and NOT a counterexample to its open constant-a claim or its new window.

Fix q>=3 and define the configuration colour r(x)=sum_i x_i mod q.
For the same unmarked Hamming reference K_N set unnormalized scaling

    w_x=N if r(x)=0, w_x=1 otherwise,
    tilde_pi_x=w_x(K_N w)_x / Z_N,
    z_x=w_x/sqrt(Z_N).

Every colour has q^(N-1) configurations. From each state exactly N legal
neighbours have each different colour. Therefore the stationary weights
per labelled state, before the common normalization, are

    A_high=qN^2,
    A_low=N^2+(q-2)N+1=:L_N.

Thus 1<=A_high/A_low<q for N>1 (equality to 1 at N=1). EVERY adjacent
ratio and its reverse are bounded by the fixed C=q. The source scaling
is exact by construction and is the unique entropy-selected scaling for
this supplied diagnostic pi; no numerical root or extra history state is
used. Slot exchangeability holds. For q=12 and N in12N, global cyclic
label translation preserves r, and reflection sends r to -r, so the
stationary input and source also have the declared D12 label covariance.

For a legal edge between two nonzero colours, the actual normalized edge
factor is exactly

    v_x v_y=d_N/L_N ~ (q-1)/N ->0.                  (1)

The self probability at a low colour is 1/L_N=Theta(N^-2), and the maximum
neighbour ratio z_x/z_y=N. Thus the orders O(N), O(N^2) and Omega(1/N)
in IS-001 cannot generally be replaced by constant scaling or larger
orders using only these bounded-ratio hypotheses. This establishes order
sharpness in THIS BROADER CLASS, not optimality for smooth quadratic A7.

For N in24N one can even place low-low test edges inside the actual MT
geometric cylinder eventually: from the balanced count n_* use
x=n_*-e_0+e_1 and y=n_*-e_0+e_2. Their colours are 1,2, they are adjacent,
and their scaled coordinates tend to zero, hence actual MT compact-box
inclusion contains them. The DIAGNOSTIC stationary measure is nevertheless
different from W_A7. No same-Gibbs/KR/IBA counterfamily is claimed.

This distinction can also be proved by an exact ordering, not merely by
naming the inputs differently. At n_*, A7 n_*=0. For delta=e_1-e_0,
the positive parity eigenvalue gives delta^T A7 delta>=4 Lambda6/q>0.
Thus the actual labelled Gibbs weight at x exceeds that at n_* when
g>0; at g_c their exponent difference is at least 2/N. In this control,
the same two types have colours 1 and 0, so the labelled stationary weight
at x is SMALLER than at n_*. These opposite orderings rule out equality
with the actual A7 stationary input on this subsequence.

The q-colour projected source is also explicit. Its high row is uniform
1/q. For a low colour a,

    P_aa=1/L_N, P_a0=N^2/L_N,
    P_ab=N/L_N for low b!=a.

These rows and detailed balance provide small exact verification without
enumerating q^N states. The colour is a deterministic function of the
existing configuration, not an added physical phase/controller/hidden bit.

The actual A7 source has a smooth mean-field quadratic exponent, which
this microscopic-colour input does not. Using the control to refute that
actual source or silently replacing its Gibbs law would be invalid.
The original stronger constant-a target for A7 remains explicitly OPEN.
