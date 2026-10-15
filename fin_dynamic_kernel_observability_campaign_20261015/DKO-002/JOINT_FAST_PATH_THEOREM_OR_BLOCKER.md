# DKO-002 — Conditional joint fixed-fast-horizon theorem

**Claim (conditional analytical theorem).** Fix ε∈{0,+.01}, T<∞, and a deterministic sequence k_N/N→α∈(0,1). Under the original iid-uniform labelled q12 LOO-HB parent at exact g_c, observed from a common deterministic post-MT t_N within the accepted local sector and before its paid exit window, the process

(Z_N(t_N+·),F_{1,3,4;N}(t_N+·),U_{1,3,4;N}(t_N+·))

converges weakly on D([0,T],R^7) to (Z_0,F_{1,3,4},U_{1,3,4}), where Z_0 has the inherited quartic distribution and is constant, and (F,U) is the unique stationary, mutually independent Gaussian OU solution

dF_r=−γ_r F_r dt+√2 dB_r,  dU_r=−U_r dt+√(2c_α)dW_r,

F_r(0)~N(0,1/γ_r), U_r(0)~N(0,c_α), c_α=(1−α)/α.

The six independent Brownian drivers are independent of the joint initial law. They are **emergent limiting martingales**, not inserted microscopic noise. The result holds *pointwise in each fixed α*, not uniformly for α→0.

## Proof argument and dependencies

**1. Initial law:** KO-002 gives the joint central CLT for all three F/U coordinates and the critical Z with limiting factorization under conditional local Gibbs, for any fixed α. Conditional label exchangeability supplies the block sampling kernel. MT gives labelled iid-start transfer after local relaxation. No finite-N independence is required.

**2. Characteristic control:** DKO-001 derives exact generator drift/brackets and gives integrable remainder estimates under stationarity of the reflected labelled lift. On every fixed [0,T] the mean absolute accumulated drift errors and bracket errors go to zero; hence convergence in probability. All jump sizes ≤C_α/√N. Local Gibbs polynomial moments control stopped drifts, and quadratic variations are bounded by constants depending only on α, A and T (since f and per-slot rate are bounded), yielding martingale square-integrability and bounded expected quadratic variations. Apply Doob's L2 bound and Aldous' increment criterion to stopped semimartingales; remove localization using MT's uniform moments and the exit coupling.

**3. Identification:** Every subsequential limit solves the martingale problem with generator

G=Σ_r[−γ_r f_r ∂_{f_r} +∂²_{f_r} −u_r ∂_{u_r}+c_α∂²_{u_r}],

with Z held constant and the joint product-Gaussian initial law. The linear globally Lipschitz drift and constant positive diagonal diffusion have a unique path law (explicit variation-of-constants/Gaussian construction); the martingale problem thus determines the full sequence's limit. A zero limiting cross bracket *plus joint continuous Gaussian martingale identification* pays independence; zero covariance alone would not.

**4. Parent transfer:** under the accepted MT stationary reflected coupling, all relevant path laws are TV-close until a fixed T added to t_N. The actual parent is never reflected or conditioned; exceptional paths are retained as errors.

**5. Memory:** V_r=F_r+U_r has covariance C_r(t)=γ_r^-1 e^-γ_r|t|+c_α e^-|t|, since the two stationary OU components are independent. Unless γ_r=1, α=1, or a degeneracy kills a component, V_r alone is generally NOT Markov: the sum has two distinct exponential relaxation rates. It is a Gaussian process with a latent two-dimensional Markov representation, not a newly postulated scalar Langevin process.

## Proof-status ceiling

This is a self-contained **conditional analytical argument** after importing the admitted MT/KO critical local Gibbs, labelled coupling, and comparator extension. No full source reviews or independent expert referee verified the analytic inequalities in this current package. It is not an unconditional mechanically certified finite-N theorem, a growing-T limit, a path-space total variation equivalence, or a finite-N calibrated rate.
