# NT-002 — Strong lifted transfer from the true reversible process

**Theorem (bounded-continuous lifted core; conditional on S1–S4).** Put `H_N=L²(D_N,pi_D,N)`, `H=L²(R,rho_lambda)`, `J_N f(n)=f(z_N(n))`, `S_N(t)=exp(sqrt(N)t L_D,N)` and `S(t)=exp(t G_lambda)`. Then for each fixed `t>=0` and each `f∈C_b(R)` (complex permitted),

`||S_N(t)J_N f-J_N S(t)f||_{H_N} -> 0`.

Also `||J_N f||_{H_N}->||f||_H`. The theorem **does not require** `z_N` to be an exact finite-N Markov quotient.

**Proof.** Norm identification is `int |f|² dnu_N ->int |f|² d rho` by `S1`. Set `h=S(t)f`; by the limiting nonexplosive Feller diffusion, `h∈C_b`. The square is exactly

`||S_N(t)J_N f-J_N h||²
 = <J_N f,S_N(2t)J_N f> -2 Re <J_N h,S_N(t)J_N f> +||J_N h||²`.

This uses selfadjointness of `S_N(t)`, the full native semigroup and no lumpability. The first term converges, via the **stationary (0,2t) pair**, to `<f,S(2t)f>_rho=||S(t)f||²_rho`. The second converges, via the **stationary (0,t) pair**, to `<h,S(t)f>_rho=||h||²_rho`. The third converges by `S1` to `||h||²_rho`. All three bounded continuous functions are admissible in the claimed weak convergence. Their signed sum is zero. QED.

**Corollary (lifted resolvents).** For every fixed `r>0`, let `R_N(r)=(r-sqrt(N)L_D,N)^(-1)`, `R(r)=(r-G_lambda)^(-1)`. For `f∈C_b`,

`||R_N(r)J_N f-J_N R(r)f||_{H_N}->0`.

The resolvents are the Laplace transforms `int_0^∞ e^{-rt}S_N(t)dt` and `int_0^∞e^{-rt}S(t)dt`. The integrand's norm tends to zero at each `t`, and its norm is uniformly bounded by `2||f||∞`; dominate by `e^{-rt}` and use Bochner/Minkowski. `R(r)f` is bounded continuous by the Feller property and dominated time integral. `r=0` is not licensed.

**Precise negative results.** This theorem gives a lifted *core* statement; no assertion is made here of general varying-Hilbert strong/Mosco convergence (which requires a complete identification scheme for all L² sections), compact convergence, uniform-in-time errors, convergence of all finite-N eigenvalues, eigenvectors, spectral counts, gap or trace. Individual compactness of the limiting sextic Schrödinger operator does not promote the family to compact spectral convergence.
