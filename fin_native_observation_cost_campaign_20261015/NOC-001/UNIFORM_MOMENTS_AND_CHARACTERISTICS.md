# NOC-001 — Native normalized block chart and generator payment

## Exact finite-N generator (new NOC chart)
Let `b,o` be disjoint labelled block/outside q12 counts, `n=b+o`, `N=|n|`, `k=|b|`, `alpha=k/N`. Define `Q_j(v)=softmax_j((g_c/N) A v)` and `f=f_3`. The true generator acts on any test `H(b,o)` as

`L H = sum_{i,j} b_i Q_j(n-e_i)[H(b-e_i+e_j,o)-H(b,o)] + sum_{i,j} o_i Q_j(n-e_i)[H(b,o-e_i+e_j)-H(b,o)]`.

Null (`j=i`) refreshes do not change the state. Reflection blocks only increments for which `n-e_i+e_j` exits `D_N`; no renormalization of allowed refreshes.

Set `F = N^-1/2 sum_a f(x_a)`, `R = k^-1/2 sum_a (1_{a in B}-alpha) f(x_a)`, `W=sqrt(alpha)F+R`. Every actual jump `i->j` of slot `a` changes `(F,R)` by `d_a(f_j-f_i)` with `d_a=(N^-1/2, (1_{a in B}-alpha)k^-1/2)`.

Let `a_f(v)=sum_j f(j)Q_j(v)`. Directly,

`L F = sqrt(N)[a_f(n)-sum_i n_i f_i/N]+N^-1/2 sum_a[a_f(n-e_{x_a})-a_f(n)]`,

`L R = -R + k^-1/2 sum_a (1_{a in B}-alpha)[a_f(n-e_{x_a})-a_f(n)]`.

Because softmax is globally Lipschitz over this bounded fixed-operator set, `|a_f(n-e_i)-a_f(n)| <= C/N` and `sum_a |1_B-alpha|=2k(1-alpha)`. Thus `|L R+R| <= C sqrt(k)/N` **uniformly in configuration**, stronger than a spurious inverse-alpha estimate.

## Stationary moments without a shrinking-block CLT
Under reflected `pi_D`, the total count has the inherited MT local-Gibbs weight and, conditional on `n`, block assignment is exactly multivariate hypergeometric. Thus `E(R|n)=0` and

`E(R^2|n) = (N-k)/(N-1) [sum_j p_j f_j^2 - (sum_j p_j f_j)^2] <= C`.

The inherited MT fast-moment bound gives `E_pi F^2<=C` for r3 (with the proper 12-to-Fourier normalization). Therefore `E_pi(F^2+R^2)<=C` uniformly `1<=k<=N/2`. No transport of this unbounded moment through a TV comparison is performed.

## Bracket density from the native jump rule
Define `J_i(n)=sum_j Q_j(n-e_i)[f_j-f_i]^2`. For all configurations,

`Gamma_FF = N^-1 sum_a J_{x_a}`;
`Gamma_RR = k^-1 sum_a (1_B-alpha)^2 J_{x_a}`;
`Gamma_FR = (Nk)^-1/2 sum_a(1_B-alpha) J_{x_a}`.

As `p -> u`, `J_i(n)=1+f_i^2+O(||p-u||+1/N)`; **not** identically 2 for a given old label. The uniform mean of `1+f_i^2` is 2. Under `pi_D`, `E|p-u|<=C N^-1/4`, while conditional hypergeometric sampling gives `E|p_B-p|<=C k^-1/2` (coordinatewise). It follows that

`E_pi|Gamma_FF-2| <= C (N^-1/4+N^-1)`;
`E_pi|Gamma_RR-2(1-alpha)| <= C (N^-1/4+k^-1/2)`;
`E_pi|Gamma_FR| <= C (N^-1/4+k^-1/2)`.

The last is deliberately conservative (the conditional exchangeability calculation often yields O(N^-1/2)). Mixed brackets must be paid, not assumed independent merely because the stationary covariance vanishes.

For any vector jump `Delta=(Delta F,Delta R)`, `|Delta|<=C(1/sqrt(N)+1/sqrt(k))`. The absolute third-jump compensator obeys `sum_rates E|Delta|^3 <= C(k^-1/2+N^-1/2)`: for inside slots the contribution is `k*k^-3/2`, for outside slots `N*(sqrt(k)/N)^3`; the corresponding F terms are O(N^-1/2).

## Global drift and conditional scope
On the exact parity line `p_m=u+m sigma/12`, `f·softmax(gAp_m)=f·p_m=0`. Uniform transverse derivative is `-gamma F`, `gamma=1-Lambda3/Lambda6`. On the MT cylinder,

`|L F+gamma F| <= C (sqrt(N)|m| ||w|| + sqrt(N)||w||^2 + N^-1/2)`.

By Cauchy–Schwarz with inherited `E m^2<=C N^-1/2` and `E||w||^2<=C/N`, the stationary L1 bound is `C N^-1/4` (plus lower orders). The DKO fixed-alpha formula is **reused algebra**, while its proof of uniformity in shrinking k is the new NOC step.

## Reflected boundary and smooth-generator class
For `v in C^3(R^2)` with bounded derivatives orders 1–3, jumps permitted by parent but suppressed by reflection differ by `<=C ||Dv||/sqrt(k)`. Their total attempt rate is at most `N`. Hence, with departure layer `B_D`,

`E_pi |(L_D-G)(v o Phi)| <= C ||v||_{C^3_1} [N^-1/4+k^-1/2+(N/sqrt(k)) pi_D(B_D)]`,

where `G=-gamma f d_f+d_ff-r d_r+(1-alpha)d_rr` and `||v||_{C^3_1}=sum_{1<=|beta|<=3} ||D^beta v||_infty`.

The boundary term uses genuine `pi_D(B_D)` rather than a freely postselected/truncated ensemble. MT's fixed local positive-energy departure boundary yields `pi_D(B_D)<=poly(N) exp(-bN)`; thus its contribution is negligible for the frozen polynomial sequence, **conditional on the accepted MT domain and bounds for each comparator**.

This is a stationary L1 generator comparison, **not** an assertion of uniform statewise error or growing-time path weak convergence.
