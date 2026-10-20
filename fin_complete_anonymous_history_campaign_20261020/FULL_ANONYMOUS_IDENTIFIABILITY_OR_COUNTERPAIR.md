# AH-001 — complete anonymous equality histories determine the jump law up to relabelling

Date: 2026-10-20. Mathematical claim **conditional on the frozen one-label Markov chain**. No new FIN law is proposed.

## Exact domain and observable

Let `V=Z/12Z`. For `1<=d<=6`, let `w_d>0`, `s=2 sum_{d=1}^5 w_d+w_6`, and define the stochastic matrix

`J_ii=0`, `J_ij=w_{min(|i-j|,12-|i-j|)}/s` for `i!=j`. Start `Y_0` uniformly and subsequently use `J`. For each `L>=1`, the only record is the restricted-growth (first-appearance) sequence `Z_0,...,Z_{L-1}`: `Z_0=0`, a previously visited true vertex receives its previous integer code, and any new vertex receives the next unused integer code. Distances, shell indices, initial real names, transition times, clock, other copies and any direct observational instrument are unavailable.

All `J_ij` off diagonal are strictly positive. In particular `J^2(i,j)>=10 a^2` for every `i,j`, where `a=min_{i!=j}J_ij>0`: for `i!=j`, there are 10 intermediate vertices different from both; for `i=j` there are 11. Thus `J^2(i,j)>= b:=10 a^2`, and with `U` the all-1/12 matrix,

`J^2=beta U+(1-beta)R`, `beta=12b=120a^2 in (0,1)`, where `R` is a stochastic matrix. The inequalities use no representation as a Markov process on FIN's interacting count space.

## Theorem AH1 (isomorphism recovery)

**Statement.** For any two frozen `J,K`, if all finite-prefix probability distributions of the corresponding anonymous codes `Z` coincide, there exists a 12x12 permutation matrix `P` with `K=P J P^T`. Conversely such conjugacy yields identical anonymous histories (including their finite-prefix laws). In fact, from one infinite anonymous sample path, a measurable estimator recovers the conjugacy class almost surely.

**Proof.**

1. **Coverage.** For each vertex `v`, as long as `v` has not appeared, the next transition goes to it with conditional probability at least `a`. Hence `Pr(v is unvisited by time L-1) <= (1-a)^(L-1)`. By a union bound, `Pr(all 12 seen by L-1) >=1-12(1-a)^(L-1)`, and almost surely there is a finite cover time `tau`.

2. **State and edge frequency convergence.** The `J^2` minorization above implies uniform geometric convergence towards the unique stationary uniform measure `pi`. For a bounded mean-zero state function `f`, `||J^(2k) f||_infty <= 2(1-beta)^k ||f||_infty`; consequently `g=sum_{r>=0} J^r f` converges uniformly and `g-Jg=f`. The decomposition

   `sum_{t=0}^{n-1} f(Y_t)=g(Y_0)-g(Y_n)+sum_{t=0}^{n-1}(g(Y_{t+1})-Jg(Y_t))`

   has bounded boundary terms and bounded martingale differences. Azuma-Hoeffding plus Borel-Cantelli gives `n^(-1)sum_{t<n}f(Y_t)->0` almost surely. Take `f_i=1_{Y=i}-1/12` to obtain `N_i(n)/n->1/12` for every `i`. Separately, for each ordered `i,j`, the edge-martingale differences

   `1_{Y_t=i,Y_{t+1}=j} - 1_{Y_t=i}J_ij`

   are bounded by 1 and centered conditionally on `Y_0,...,Y_t`, so their empirical mean tends almost surely to zero. Therefore `N_ij(n)/N_i(n)->J_ij` for each `i,j`.

3. **Random first-visit naming.** After cover time `tau`, the names assigned by the observable `Z` define a random bijection `sigma:{0,...,11}->V`. It is fixed forever after its creation. Let `N_h(n),N_hk(n)` be visits and consecutive-edge counts computed **solely from the code** `Z_0,...,Z_n`. Then, simultaneously for all `h,k`,

   `\hat J_n(h,k):=N_hk(n)/N_h(n) -> J_{sigma(h),sigma(k)}` almost surely.

   Define `\hat J_n(h,k)=0` when a row has zero visits, and pad missing names with zero rows; this makes it a measurable function of a prefix in every case. Because coverage is almost sure, the padding is eventually irrelevant. Thus `\hat J_infty` is a measurable function of the entire anonymous sequence and lies in the finite permutation orbit `O(J)={PJP^T}` almost surely.

4. **Passing from cylinder laws to the orbit.** Equality of all finite-prefix laws on the finite code alphabet gives equality of the induced probability measures on the infinite sequence space by uniqueness on cylinder sets. The distribution of the measurable limit `\hat J_infty` is therefore identical under both models. Under `J` this distribution has probability 1 on `O(J)`; under `K`, probability 1 on `O(K)`. Thus the common distribution gives probability 1 to `O(J) intersect O(K)`; the intersection cannot be empty. Two conjugacy orbits either coincide or are disjoint, so `K=PJP^T`. Conversely, relabel every hidden vertex by `P`; all equality partitions remain unchanged, and uniform initial laws are transported correctly. QED.

**Additional point:** every first-visit bijection actually has positive probability, since the chain can visit 12 distinct vertices in any chosen order with positive probability. This makes the limiting estimator's distribution support equal to `O(J)`; the main proof needed only one shared orbit member.

**Scope:** this establishes completeness of the **entire collection** of finite anonymous equality-history laws, not of any fixed low-order moment, return trace, four-event motif or alternating-moment family. It is an inverse result for an assumed model; it does not produce the path law or choose the radial kernel, activity or a physical clock.
