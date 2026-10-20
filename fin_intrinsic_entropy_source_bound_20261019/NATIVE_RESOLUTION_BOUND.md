# ID-001 — native FIN N2 reset law and controlled entropy expansion

## Exact inherited target (not a selected physical law)

For `q=12`, symmetric circulant A7 with first row `a_j=A_{0j}`, and a fixed finite real gain g, let

`p_j=exp[(g/2)a_j] / Σ_{v=0}^{11} exp[(g/2)a_v]`.

Then p_j>0, Σp_j=1. This p is *not* set equal to uniform for positive g. The actual N2 difference quotient of the supplied heat-bath parent has generator (acting on row-transition probabilities)

`B_{ij}=2(p_j−1{i=j})`, equivalently `B=2(1 p^T−I)`.

Its exact semigroup is `P_δ=e^{δB}=e^(−2δ) I+(1−e^(−2δ)) 1 p^T`, and p is stationary: pPδ=p. The q12 label differences are internal coordinates, not measured positions in emergent space. These facts are inherited from the source-bound N2 process, not new derivations of physical kinetics.

## Exact one-step entropy, including diagonal terms

Put `t=1−e^(−2δ)` and `a_i=1−p_i`. Then `Pδ(i,i)=1−t a_i`, `Pδ(i,j)=t p_j` for j≠i. Thus the stationary entropy rate per sampled symbol equals

`hδ =− Σ_i p_i [(1−t a_i)log(1−t a_i)+Σ_{j≠i} t p_j log(t p_j)]`.

Equivalently define `A_p=1−Σ_i p_i²`, `C_p=Σ_i p_i(1−p_i) log p_i (<0)` and obtain the useful **exact identity**

`hδ=−Σ_i p_i(1−t a_i) log(1−t a_i) − t A_p log t − t C_p`.

This formula is a separate check: it keeps the `no state change` probability intact, including refresh-to-same-state events.

## Theorem: small-δ asymptotics

For every fixed strictly positive p (equivalently, fixed finite g), let

`r̄=2 A_p=2(1−Σ_i p_i²)>0`,

`c(p)=r̄(1−log2)−2 Σ_i p_i(1−p_i)log p_i`.

Then as δ↓0,

`hδ = r̄ δ log(1/δ) + c(p)δ + O_p(δ² |log δ|)`.

Proof with explicit bounded remainder route: for `0≤u≤1/2`, `|(−(1−u)log(1−u))−u|≤u²` (integrate `−log(1−v)` between 0 and u and use `−log(1−v)≤v/(1−v)`). Since `t≤2δ≤1/2` for δ≤1/4, the diagonal term is `t A_p + O(t²)` with an absolute error at most t². Also `|t−2δ|≤2δ²`, while `log[t/(2δ)]=O(δ)` by the elementary bounds for `1−e^(−2δ)`. Substitute in the exact identity:

`hδ=t[A_p(1−log t)−C_p]+O(t²)`
`=2δ A_p log(1/δ)+2δ[A_p(1−log2)−C_p]+O_p(δ²|logδ|)`.

The coefficient r̄ is independently the expected **actual state-change rate**, `Σ_i p_i Σ_{j≠i} B_ij=2Σ_i p_i(1−p_i)`. This is not the gross attempted refresh rate 2, since refreshes to the same label make no state change.

No uniform p_min is asserted when g runs over an unbounded range; `O_p` may depend on fixed p, though the elementary diagonal part has a uniform bound. This is an expansion of sampled-process Shannon/KS entropy, not of Lecomte et al.'s *regularized jump-sequence entropy*.

## Conditional application to a supplied QW R26 Hamiltonian

The published QW template, if adopted as a model-time Hamiltonian flow, reads `X=(q,p)∈R^26`, `dot q=p`, `dot p=−∇V(q)` with a fixed polynomial V. No particular QW coefficients, compact invariant energy level, or physical clock are certified by ID-001.

**Conditional corollary.** If a compact forward-and-backward invariant K⊂R^26 is actually supplied, and the vector field F is C¹ on a neighborhood of the convex hull of K (true for a fixed polynomial), then for any invariant probability µ supported in K the time-one map has finite metric entropy. Indeed with `L=sup_{conv(K)} ||DF||<∞`, Grönwall gives a time-one Lipschitz bound `Lip(T1|K)≤e^L`. Covering a compact subset of R^26 at radius ε requires `O_K(ε^(-26))` balls. Initial εe^{−L(n−1)} nets give (n,ε) spanning sets for T1, proving `h_top(T1|K)≤26 L`, and `h_µ(T1)≤h_top<∞` by the variational inequality.

Therefore **if** the same K, µ, and polynomial Hamiltonian flow are offered as a deterministic generator of the N2 heat-bath skeleton with full stationary joining and one shared supplied time parameter for every δ=1/m, the entropy-budget theorem applies and forbids exact all-resolution generation. This is a CONDITIONAL class obstruction, *not* evidence that real FIN chooses this flow. An existing coercive conserved Hamiltonian can establish compact K if the hypotheses are explicitly met; the source family does not establish them for every signed coefficient tuple. No inference is made about the infinite-dimensional NS/RP PDE.
