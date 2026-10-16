# PC-001/PC-002 — quantitative *restricted* weighted-Hall theorem

Date 2026-10-17. **Status: a new conditional analytical lemma, NOT the requested full-Hall theorem.** The compact phase-prefix/current asymptotic is an inherited premise described in the supplied order; its original source proof (PC03/PC04) was not available as bytes here. All finite algebra and geometric facts explicitly paid below are independent of that unavailable source.

## 0. Model and exact target

Fix the actual 12-label A7 count model at `g_c=12/Lambda6`, `N in 12ℕ`, MT reflected D_N and actual departure layer E_N. `m=N^{-1}Σ(-1)^j n_j`, `z=N^{1/4}m`, `C=Σ n_j cos(πj/2)`, `S=Σ n_j sin(πj/2)`, and `R=2(C²+S²)/N`. The frozen tie rule is A if `R≤T(z)` and B otherwise. The permitted upward graph H includes all within-phase edges and cross-phase edges incident to E_N. The first band has successive level parameters z_l,z_{l+1}∈[-1,1].

At a level l, `W_{P,l}=Σ_{x∈P∩X_l}W_x`, `s_x=(N/2)W_x (L0^D m)(x)`, `J_{P,l}=Σ_{x∈P,m(x)≤m_l}s_x`, `I_{P,l}(x)=J_{P,l-1}W_x/W_{P,l}`. Let `a_x=[-s_x-I_{P,l}(x)]_+` for x∈E_N∩P∩X_l; `δ_{P,l}=Σ_E a_x`. The frozen corrections are `I'_E=I_E+a`, `I'_bulk=I_bulk−δ_{P,l}W_x/W_{P,l}^{bulk}`, `O'_x=I'_x+s_x`. By premise PC04 they are nonnegative and sum to the exact phase current.

Full PC-002 asks for all `T⊆X_l`:

`O'_l(T)≤I'_{l+1}(N_H(T))`.

## 1. Uniform exact *single-edge* count-weight ratio (new finite-N proof)

On the real MT cylinder `|m|<1/100`, `||w||_2<1/10000`, where `p=n/N=p_m+w`, `p_{m,j}=(1+m(-1)^j)/12`, the bounds are `p_min=103/1250` and `p_max=158/1875`. For an upward edge `x=n -> y=n-e_a+e_b`, a odd and b even, both in D_N, write `v=e_b−e_a`. Since A7 is positive semidefinite and `g_c A7 p_m=m σ` modulo a constant vector,

`W_y/W_x = [n_a/(n_b+1)] exp{g_c v·A7 p+(g_c/2N) v·A7 v}`.

The second term in the exponent is nonnegative. `|g_c v·A7 w|≤sqrt(2)g_c Λ6||w||≤12 sqrt(2)/10000`, and `g_c v·A7 p_m=2m`. Consequently for N≥10000,

`W_y/W_x ≥ [p_min/(p_max+1/10000)] exp(−2/100−12 sqrt(2)/10000) > 9/10`.

The final strict bound is certified by `sqrt(2)<3/2` and `e^(−t)≥1−t`, wholly with rational arithmetic; the verifier records the explicit lower fraction 1511319/1581875 >.9. The reverse-edge lower bound also follows after swapping endpoints (both are in the same MT cylinder), so `W_x/W_y>.9`. **This is a weight ratio, not an assumption about a new generator's rate.**

## 2. A uniform current window, conditional on the PC03/PC04 local-limit input

Denote the critical density `ρ_0(z)=const exp(−z^4/12)`, `b_0(z)=−z³/3=ρ'_0(z)/ρ_0(z)`, and the phase functions already supplied:

`p_A(z)=1/2+e^(−z²)/8`, `p_B(z)=1−p_A(z)`, `q_P(z)=p_P(z)ρ_0(z)`,

`j_P(z)=∫_{−∞}^z b_0(u)q_P(u)du`.

By oddness of `b_0 q_P`, `j_P(z)=−∫_z^∞b_0(u)q_P(u)du` for z≥0. Since `ρ'_0=b_0ρ_0`, and `3/8≤p_P≤5/8`, both representations give `j_P(z)>0` and

`p_P,min ρ_0(z)≤j_P(z)≤p_P,max ρ_0(z)` for every real z.

For phase A, `p_min,p_max=(1/2,5/8)`, while for B these equal `(3/8,1/2)`. Hence **analytically, with no numerical integration**,

`1/5≤j_A/(4q_A)≤5/16`, `3/16≤j_B/(4q_B)≤1/3`.

The exact stationary prefix-to-local ratio asserted by the source compact local-limit input is

`J_{P,l-1}/(N W_{P,l}) -> j_P(z_l)/(4q_P(z_l))`, uniformly for z_l in a fixed compact interval. (At l+1 use its corresponding current/weight ratio.) This factor 1/4 follows by parity spacing `Δz=2N^{-3/4}`, exact source `s=(N/2)W d`, and `N^{-3/4}b_0(z)` critical drift, but **the uniform Stirling/local-lattice estimates and lattice-to-integral passage are inherited PC03, not freshly audited from unavailable source bytes**.

The MT exponential departure-layer and compact local mass bounds also imply `δ_{P,l}/(N W^{bulk}_{P,l})=o(1)` uniformly in the fixed band: `0≤a_x≤NW_x`, the numerator E_N weight has `F≥h/2`, and the fibre's phase-bulk mass has polynomially weighted `F=O(1/N)` configurations. This is an inherited MT/PC03 dependency, not a new uniform all-D estimate.

Thus there is an N0 (not numerically licensed) for which every level in the band satisfies:

`I'_l(y)≥(9/50)N W_y =.18 N W_y` for **every** y (E increases I, bulk loses o(NW));

`I'_l(x)≤(7/20)N W_x =.35 N W_x` for x outside E_N.

For nondeparture x the original nonreflected LOO drift is used. The logits differ from common parity-linear logits `m σ_j` by ≤`12δ+12/N` using `g_c||A7||=12`. The even-minus-odd choice expectation therefore has magnitude ≤`tanh(|m|+12δ+12/N)`. For N large enough that |z|≤1 gives |m|≤1/200 and 12/N≤.0012, the exact scalar drift obeys

`|L0^D m(x)|≤|m|+tanh(|m|+12δ+12/N)≤.0124`.

Therefore `|s_x|/(NW_x)≤.0062<.007` and, for x not in E_N,

`0≤O'_l(x)≤(357/1000)NW_x`.

The common rate bound alone gives `|s_x|≤NW_x` at **all** states including E_N. At E, the frozen correction gives `O'_x=max(I_x+s_x,0)≤I_x+NW_x`, so if an E vertex happens to have nine legal successors the singleton proof below also works with upper `(27/20)NW_x`. It is NOT asserted that every true departure vertex has nine successors.

## 3. Nine phase-safe upward successors at every x outside E_N

If `|R−T(z)|` exceeds the inherited one-step shell width, all 36 native upward edges remain in the same phase. Otherwise, the SC exact Fourier-step argument applies: upward odd-to-even moves change `(C,S)` in one of four vectors `(±1,±1)`, with 9 distinct label pairs for each vector. Near the phase threshold R is bounded away from zero, and one radially inward (phase A) or outward (phase B) vector remains same-phase despite the `T(z+Δz)` change. Since x∉E_N all original proposals stay within D_N. Thus **at least nine distinct H-upward successors** can be selected for every x outside E_N and adjacent levels in the prescribed band. Their weights satisfy `W_y>(9/10)W_x`. Same-phase successors are sufficient for this argument even though full H also allows E-cross-phase edges.

## 4. New *restricted* weighted-Hall result, without changing the carry

For any nondeparture singleton `T={x}`, using the nine distinct selected successors,

`I'_{l+1}(N_H({x}))≥9*(9/50)*(.9)NW_x=(729/500)NW_x`;

`O'_l(x)≤(357/1000)NW_x`.

Thus the singleton Hall inequality holds with ample positive slack at every nondeparture left vertex, uniformly on the fixed critical band for sufficiently large N.

More generally, for a set T of nondeparture left vertices, choose **the prescribed nine safe successors for each x**, and suppose the number of chosen left predecessors of each right vertex is at most four. Weight and positive-column estimates yield

`I'_{l+1}(N_H(T))≥[(9/50)*(9/10)*9/4]N W(T)=(729/2000)NW(T)`.

But `O'_l(T)≤(357/1000)NW(T)`, and `729/2000−357/1000=3/400>0`. Hence **weighted Hall is proved for every such bounded-overlap subset**. The overlap bound is a *verifiable property of the particular T and its selected nine edges*, not an axiom imposed on all T.

The same local proof can be used for other fixed overlap bounds if sharper current/drift controls are independently paid. **There is no proof that all subsets of X_l admit overlap ≤4.** Indeed dense sets can have many shared successors.

## 5. Exact nonconclusion and terminal open atom

The existence of nine legal successors, strict Gibbs comparability, and Hall for singletons/bounded-overlap T does **NOT** establish Hall for arbitrary T⊆X_l. A separate exact toy fixture shows two individually feasible singletons but a failing pair; this is a logical control, not an actual-FIN adverse family. No actual growing-N violating T was produced. The general PC-002 result stays `OPEN_ADJACENT_FIBRE_COUPLING`, and PC-003 is `NOT_ACTIVATED`. The missing lemma is a genuine **overlap/congestion bound for dense subsets in the actual H graph with the fixed I',O'**, or an explicit actual-N family violating Hall. No new kinetics or physical source is inferred.
