# LI-001 — Exact one-tag feedback-removal inequality (finite N)

**Disposition:** EXACT FINITE PATH-RELATIVE-ENTROPY AND TV THEOREM. This is a comparison of two laws, not a modification of FIN. The result is already supplied in the user's bounded lemma; the derivation here independently checks its constants, direction and scope.

Let the tagged full q12 label be `b`, background type vector `n` with `|n|=N-1`, and original gain `g`. The original P labelled rate-one heat-bath attempts have tag refresh `p_b(n)=softmax((g/N)An)` and background slot i refresh `p_i^P(n,b)=softmax((g/N)A(n-e_i+e_b))`. The auxiliary R has the SAME initial joint measure, the SAME tagged refresh, the SAME N denominator and rate-one clocks, but background refresh `p_i^R(n)=softmax((g/N)A(n-e_i))`. The state is the *fully labelled* chain, with unobserved, rate-one independent Poisson attempts including null refreshes. This construction is valid for any initial joint law, correlated or not.

For `F(h)=log ∑_j exp(h_j)`, `P_h=softmax(h)`, and `δ`,

`D(P_(h+δ)||P_h)=F(h)-F(h+δ)+∇F(h+δ)·δ=∫_0^1 t δ^T Hess F(h+tδ)δ dt`

`=∫_0^1 t Var_(P_(h+tδ))(δ_j) dt <= (osc δ)^2 /8`, since `Var(X)<=range(X)^2/4`. At a changed background attempt, `δ=(g/N)Ae_b`; `osc δ<=2 |g| M_A/N`, where `M_A=max|A_ij|`. Hence `D(p_i^P||p_i^R)<=g² M_A²/(2N²)` uniformly in the complete labelled configuration. Conditional on the common attempt-time and attempted-slot schedule, the chain rule for categorical refresh marks yields

`D(P_marked[0,T] || R_marked[0,T]) <= (N-1)T g² M_A²/(2N²)`.

The initial entropy term is zero because the initial laws are *identical*; the tag's mark-conditional law is identical in P and R. Forgetting the unobserved attempt marks contracts KL; Pinsker's inequality with `TV=sup_A|P(A)-R(A)|` yields

`TV(P_path[0,T],R_path[0,T]) <= |g| M_A sqrt((N-1)T)/(2N)`.

For `T=sqrt(N) S` and fixed S, this is at most `|g|M_A sqrt(S)/(2N^(1/4))`. This comparison covers the **complete** labelled path through the future target time `S=1`, and hence also its background/sign/path sigma-field projections. It does not let the observer use marked attempts. For any finite target and observation experiment with bounded loss in [0,1], the two optimal Bayes risks differ by at most this TV quantity: under any decision rule its loss expectation changes by at most TV; take the infimum in both directions.

**No-go boundaries:** the result does not decouple correlated initial variables; does not show likelihood UI, exact parity-only lumpability for P, full-record sufficiency, source of A7, physical locality, physical clock or observation device. For N=1 the bound is 0 and the background is empty.
