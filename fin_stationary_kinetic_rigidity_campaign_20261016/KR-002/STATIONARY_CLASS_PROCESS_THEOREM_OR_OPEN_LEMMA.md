# KR-002 — Stationary critical parity diffusion under NC

## Conditional theorem (proved analytically from the frozen assumptions)

Assume the KR-001 class, plus the admitted original MT stationary moment/geometry theorems and UF-002 **native exact count drift** expansion. Let `X~_N` start in exactly `pi_D,N`. For each fixed finite `S`,

`Z_N(s)=N^{1/4}m(X~_N(sqrt(N)s)), 0≤s≤S`

converges weakly in `D([0,S],R)` to the unique stationary solution

`dZ_s=(chi lambda Z_s−Z_s³/3) ds+sqrt(2) dB_s`, `Z_0∼rho_lambda(dz) ∝ exp(chi lambda z²/2−z⁴/12) dz`.

Only the bracket concentration `NC` is required beyond stationary Gibbs, reversibility, exact statewise scalar drift matching, single-copy support and uniform `CN` rate cap. The theorem is *not* an all-kinetics statement without NC.

### 1. Native drift transfer without replacing a conditional mean

Write `p=u+m sigma/12+w`, `w` in the ten-dimensional fast tangent. At a native unreflected interior state UF-002 yields for fixed cylinder and critical gain, with `b(z)=chi lambda z−z³/3`,

`|N^{3/4} L^0_N m−b(N^{1/4}m)| ≤ C N^{3/4} [ ||w||² + (|g_N−g_c||m|+|m|³)||w|| + |m|^4 + |g_N−g_c|m² + |g_N−g_c|² |m| + N^{-1} ]`.

Exact KR class equality `L~_N m=L^0_N m` transfers this function of the **entire count state** without averaging and without estimating the diagnostic's fast dynamics. Reflecting departures changes the native drift only on the boundary layer; the absolute difference is `≤ 2 C_0 1_{B_D}` from the native `C_0 N` jump-rate bound and `|Delta m|≤2/N` (at most the full rejected rate). Its integrated stationary contribution is `O(N^{3/4} pi_D(B_D))`. Weighted by `|Z_N|≤N^{1/4}` it is `O(N pi_D(B_D))`. MT's exponential boundary bound makes both vanish.

Define the scaled drift defect `e_N(n)=N^{3/4}L~_N m(n)−b(N^{1/4}m(n))`. By MT's joint stationary polynomial moments `E|m|^p=O(N^{-p/4})`, `E||w||^p=O(N^{-p/2})`, the six interior terms have expected orders respectively

`N^{-1/4}, N^{-1/2}, N^{-1/4}, N^{-1/4}, N^{-1/2}, N^{-1/4}`.

The analogous estimates with an additional `1+|Z_N|` still tend to zero by higher MT moments (Hölder). Thus

`E_pi[(1+|Z_N|)|e_N|]→0`, `E∫_0^S(1+|Z_N(s)|)|e_N(X_{sqrtN s})|ds→0`

by **exact stationarity**. This is the actual paid error, not a supposed deterministic fast tube for the new generator. Rare tails outside the UF fixed cylinder are controlled by MT polynomial/exponential moments together with the global bounded microscopic drift and the polynomial scaling factor.

### 2. Bracket identity and NC

For any finite reversible stationary count CTMC (indeed stationarity suffices), invariance of `m²` gives the exact identity

`E_pi Gamma~_N(m,m)=−2 E_pi[m L~_N m]`.

Exact drift matching and the same invariant law imply

`E_pi G_N = E_pi[N Gamma^0_N(m,m)]`.

The right side tends to **2**: native interior local q12 bracket is continuous at `(u,g_c)` with value 2, `p→u` in local Gibbs probability, bounded rate cap gives UI, and reflected boundary mass vanishes. Hence `E_pi G_N→2`. The class cap and jump bound yield pointwise `0≤G_N≤4C`.

NC and stationarity yield for every fixed S:

`E sup_{0≤s≤S}|∫_0^s G_N(X~_N(sqrtN v))dv−2s|`
`≤ S E_pi|G_N−E_pi G_N|+S|E_pi G_N−2|→0`.

This is the actual predictable quadratic-variation convergence. No claim that stationarity of the **mean alone** entails timewise concentration is made.

### 3. Native martingale, vanishing jumps and maximum containment

Let `M_N(s)=Z_N(s)−Z_N(0)−∫_0^s[b(Z_N(v))+e_N(X~_N(sqrtN v))]dv`.

It is a martingale in the **full natural count filtration**. Its bracket equals the preceding integral. Each jump is at most `2N^{-3/4}`. Its scaled total jump rate is bounded by `C N^{3/2}` and hence the absolute third jump characteristic is `≤8 C S N^{-3/4}`. Its square/bracket martingale is uniformly integrable after localization and on bounded intervals because `〈M_N〉_S≤4CS` deterministically and jumps vanish.

For tight maximum containment, use the exact jump Itô formula for `Z_N²`, stationarity and `2z b(z)=2 chi lambda z²−(2/3)z⁴≤C_lambda`. The bound on `E|Z_N e_N|`, `E Z_N(0)²` and the bracket cap yield

`E sup_{s≤S}|Z_N(s)|² ≤ C_{S,lambda} + E sup_{s≤S}|2∫_0^s Z_N(v-) dM_N(v)|`.

Doob/BDG and stationarity give `E sup|∫ Z_- dM| ≤C (E∫_0^S Z_N(v)² G_N(v)dv)^{1/2} ≤ C (4CS E_pi Z_N²)^{1/2}`, uniformly bounded. Thus `sup_N E sup_{s≤S}|Z_N(s)|²<∞` up to finitely many N (use a harmless bound for small N). This avoids postulating a diagnostic fast gap. The compensator includes the tiny jump correction already in Γ.

For Aldous increments, stop at `|Z_N|≥R`; `b` is bounded on `[-R,R]`, integrated `|e_N|` tends to zero in L1 on the whole horizon, and martingale increments have variance at most `4C δ`. The marginal tightness, vanishing jumps and maximum containment give C-tightness of `Z_N` and joint tightness of `(Z_N,M_N,〈M_N〉)`.

### 4. Joint-filtration identification, initial law, uniqueness

MT-003 yields `Z_N(0)⇒rho_lambda` and bounded polynomial moments. For any bounded continuous functional of a finite set of **past times** of the joint process, martingale increments and the square-martingale increments pass to subsequential limits by the deterministic bracket caps, uniform square/fourth-moment controls and vanishing jump bound. The limit continuous `M` is a martingale with `〈M〉_t=2t` in its *joint* natural filtration and

`Z_t=Z_0+∫_0^t b(Z_s)ds+M_t`.

Lévy characterization identifies `M=√2 B` in that filtration. The cubic-inward drift is locally Lipschitz and nonexplosive (`2z b(z)+2≤C_lambda`), so the martingale problem has a unique law. The limiting quartic `rho_lambda` is its invariant law by the zero-flux/Fokker–Planck identity. Every subsequence has the same limit, completing weak convergence. No coupling of new diagnostic iid initial data is asserted.

## Exact scope

This is a conditional analytic proof, not formal verification. There is no explicit N0, rate of weak convergence, cross-preparation statement, labelled block path law, complete operator limit, uniform growing-horizon result, physical time, or independent kinetic source.
