# NT-003 — The actual critical sextic transfer operator

Let `a=chi*lambda` and `U_lambda(z)=z^4/12 -a z²/2`, `rho_lambda(z)=Z_lambda^-1 exp(-U_lambda(z))`. The accepted diffusion generator on the real line is

`G_lambda f=f''-U_lambda' f'`, with `U_lambda'=z³/3-a z`.

Its reversible nonnegative Dirichlet form is `E_lambda(f,g)=int conjugate(f') g' rho_lambda dz` on the closure of `C_c^∞(R)` in the weighted first-order Sobolev/energy norm. Since the polynomial drift is smooth and dissipative, this is the regular conservative diffusion form of `-G_lambda`; its closure specifies the Markov realization (no boundary conditions inserted at infinity).

Define the unitary `U: L²(rho dz)->L²(dz)` by `(Uf)(z)=rho_lambda(z)^(1/2) f(z)`. The positive form in Lebesgue space is

`h_lambda(psi)=int |psi'(z)+(U_lambda'(z)/2)psi(z)|² dz` initially on `C_c^∞`; its closure yields the unique chosen Friedrichs/form realization of `H_lambda=U (-G_lambda) U^(-1)`.

Expanding `Q†Q`, where `Q=∂_z+U_lambda'/2`, gives

`H_lambda= -∂_z² + (U_lambda')²/4-U_lambda''/2
 = -∂_z²+z^6/36 - a z^4/6 +(a²/4-1/2)z² +a/2`.

For `lambda=0`, `H_crit= -∂²+z^6/36-z²/2`. Its normalized ground vector is `psi_0=Z_0^-1/2 exp(-z⁴/24)`, and `Q psi_0=0`, so `H_crit psi_0=0`. The negative quadratic portion of the potential is compatible with `H>=0` because of factorization. The ground state is unique in the supplied 1D irreducible diffusion setting; higher spectrum of the individual limit operator is discrete due to its confining sextic growth, but NO finite-family spectral convergence is claimed.

**Operational meaning:** `z` is the order-parameter amplitude, not physical space; `t` is critical model time, not seconds; `U_lambda` is an invariant-density potential, not conserved energy; `H_lambda` is a selfadjoint transfer operator, not a measured physical Hamiltonian.
