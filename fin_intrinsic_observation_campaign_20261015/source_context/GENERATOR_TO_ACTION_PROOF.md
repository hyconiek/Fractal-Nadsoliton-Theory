# TH-002 — interacting path action from the exact leave-one-out count generator

## 1. Exact finite-N nonlinear generator

Let `p=n/N` and `z_ij=e_j-e_i`. The exact count generator is

`L_N f(p)=N sum_ij p_i s_j^{(-i,N)}(p) [f(p+z_ij/N)-f(p)]`,

where

`s_j^{(-i,N)}(p)=softmax_j[g A(p-e_i/N)]`.

For a linear exponential test `f(p)=xi.p`,

`H_N(p,xi)=(1/N)e^{-Nf} L_N e^{Nf}`

is exactly

`sum_ij p_i s_j^{(-i,N)}(p) [exp(xi_j-xi_i)-1]`.

For bounded `g` and the finite matrix A, `s^{(-i,N)} -> s(p)=softmax(gAp)` uniformly. Thus

`H(p,xi)=sum_ij p_i s_j(p)[exp(xi_j-xi_i)-1]`.

It factorizes as

`H=(sum_i p_i e^{-xi_i})(sum_j s_j e^{xi_j})-1`,

and is invariant under `xi -> xi+c 1`.

## 2. Kraaij Theorem 2.3 applies

Set `v_ij(p)=p_i s_j(p)`. The coordinate softmax derivative satisfies

`sum_k |partial s_j/partial z_k|=2 s_j(1-s_j)<=1/2`.

Hence

`sup_p |p_i s_j[gA(p-e_i/N)]-p_i s_j(gAp)| <= g max_k |A_ki|/(2N)`.

The limiting intensity is continuous; it is strictly positive whenever `p_i>0`; and at a boundary face it factors as `p_i * s_j(p)` with the first factor increasing and the second positive continuous. These are exactly the boundary features required by Kraaij's finite-state Glauber theorem. Consequently, for every initial-law sequence satisfying a good LDP, the empirical path satisfies a path LDP with speed N and Lagrangian rate.

Two initial classes are explicitly paid here:

- deterministic type states `p_N -> p0`, giving `I0(p)=0` at p0 and infinity otherwise;
- the declared Gibbs equilibrium sequence. MP7-031 method-of-types bounds give its LDP with rate `F_g(p)-min F_g`, where `F_g=D(p||u)-(g/2)p^T A p` up to the fixed uniform convention.

## 3. Flux form of the exact limiting Lagrangian

Let `c_ij(p)=p_i s_j(p)` and `ell(q;c)=q log(q/c)-q+c`, with the standard conventions `ell(0;0)=0` and `ell(q;0)=infinity` for q>0. Convex duality for the edge-Poisson Hamiltonian gives

`L(p,v)= inf_{q_ij>=0: sum_ij q_ij z_ij=v} sum_ij ell(q_ij;c_ij(p))`.

This is the non-quadratic path action forced by the generator. At interior p all directed edge rates are positive, so every simplex-tangent velocity (`sum_i v_i=0`) has finite cost. At a boundary face, finite cost forbids outgoing flux from an empty state; equivalently the velocity must lie in the simplex tangent cone.

The zero-cost flux is `q_ij=c_ij`; its velocity is `s-p`, reproducing the mean-field ODE without assuming a quadratic action.

## 4. Exact reversal identity and a basin-local quasipotential

For interior p define, modulo an irrelevant constant,

`(grad F)_i=log p_i-g(Ap)_i`.

Then

`p_i s_j exp[(grad F)_j-(grad F)_i]=p_j s_i`.

Relabelling i,j yields the exact symmetry

`H(p,xi)=H(p,grad F-xi)`.

Legendre duality therefore gives

`L(p,v)=grad F(p).v + L(p,-v)`.

Since `L>=0`, every path from p_a to p_b satisfies

`Action >= F(p_b)-F(p_a)`.

If p_b lies in the deterministic basin of an interior stable minimum m, the time-reversed zero-cost relaxation trajectory from p_b toward m approaches this lower bound. Thus, for targets connected to m by the basin relaxation orbit, the basin-local quasipotential is

`V_m(p)=F(p)-F(m)`.

This is local/basin-scoped. It does not identify the global communication height, `Gamma=B4`, a capacity prefactor, or a physical least-action principle.

## 5. Gaussian truncation is not the path law

At the fold equilibrium, the actual soft theta jump exponent has second cumulant rate `2 g_f` but a nonzero fourth cumulant rate about `120.5038353`; the cubic rate vanishes there by the equilibrium i<->j symmetry. The quartic coefficient in the local Hamiltonian expansion is about `5.02099314` before any critical rescaling. Hence equal drift/covariance does not make the full jump Hamiltonian Gaussian. The diffusion in TH-005, if proved, is a critical scaling limit in which higher jump cumulants vanish—not an exact replacement of the path action.
