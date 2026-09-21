# MP7-011 — global-minimizer reduction to aligned nonnegative C4

## Exact dual/primal bridge

Let `u0=(1/12,...,1/12)`, `A7=X7 X7^T`, and `X7^T u0=0`. Define

`L(p,theta)=D(p||u0)-theta^T X7^T p+||theta||^2/(2g)`.

For fixed `p`, the unique minimum in `theta` is `theta=g X7^T p`, giving

`min_theta L = D(p||u0) - (g/2)||X7^T p||^2 = V_g(p)`.

For fixed `theta`, the finite Gibbs variational identity gives

`min_p [D(p||u0)-(X7 theta).p] = -log(mean_j exp((X7 theta)_j))`,

hence `min_p L = Phi_g(theta)`.

The simplex is compact and the entropy convention `0 log 0=0` makes the first minimization attained.  The quadratic term makes `L`/`Phi` coercive in `theta`, because the finite feature field grows at most linearly in `||theta||`. Therefore

`min_p V_g(p) = min_{p,theta} L(p,theta) = min_theta Phi_g(theta)`

and global minimizers correspond through the two minimizing conditions.

## Phase alignment in mediator coordinates

Write each conjugate Fourier pair `(cos k, sin k)`, `k=3,4,5`, in polar form with radius `rho_k>=0` and phase `phi_k`.  Its field amplitude is `a_k=sqrt(lambda_k/6) rho_k`.  Write the alternating coordinate as signed amplitude `s6`, with field amplitude `b=sqrt(lambda6/12)s6`.

The Euclidean quadratic dual cost

`||theta||^2/(2g) = (rho3^2+rho4^2+rho5^2+s6^2)/(2g)`

is independent of the three phases and invariant under label translations.  If `b<0`, an odd label translation flips its sign while rotating the paired phases and preserving the norm.  After this sign correction, MP7-008 gives

`Z(phi;b)<=Z(0;b)`.

Therefore replacing any mediator point by its aligned representative with the same radii and nonnegative alternating amplitude cannot increase `Phi_g`.

## Global-minimizer statement

For every global minimizer `theta`, the aligned representative is also a global minimizer.  Thus **there exists** an aligned nonnegative C4 global minimizer.

Moreover, because a global minimizer cannot strictly decrease under alignment, equality in MP7-008 must hold.  MP7-009/010 classify equality on every amplitude stratum: after discarding phases of zero-amplitude modes, the field is a label translation of the aligned field.  Hence **every global minimizer is D12-equivalent (indeed translation-equivalent) to an aligned nonnegative C4 representative**.

This is a theorem about global minimizers. It does not align arbitrary stationary points and does not imply a unique minimizing orbit.
