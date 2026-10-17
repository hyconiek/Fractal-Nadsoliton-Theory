# UE-002: exact noncommuting operator bounds, negative direction and controls

## A. Two unconditional implications **inside the declared operator class**

Let `V_inf=A+B`, `B=2v²Y`, `Phi_K=v tanh(kx)`. Neither implication requires `[A,B]=0`.

1. **Pointwise positive-semidefinite sufficiency.** If `A>=0` and `V_inf>=0` in Loewner order, then

   `V(x)=sech²(kx) A + tanh²(kx) V_inf >=0`.

   Therefore `L_perp>=0`. This gives non-negativity, **not** necessarily a spectral gap, full complex-sector stability or nonlinear orbital stability.

2. **Scalar Pöschl--Teller comparison.** Put `mu=lambda_min(V_inf)` and `b=max(0,lambda_max(B))`. Since `V_inf>=mu I` and `B<=b I`,

   `L_perp >= [-d_x²+mu-b sech²(kx)] I` in quadratic-form order.

   Scalar ground energy is `mu-k² s²`, `s=(sqrt(1+4b/k²)-1)/2`. Hence `inf spec(L_perp)>=mu-k²s²`. In particular, `mu>=k² s²` is a *rigorous sufficient* nonnegativity criterion, but its failure proves **nothing** by itself. No commutativity assumption was used.

## B. Exact negative-direction test for arbitrary noncommuting arrays

For any constant unit vector `u`, take `f(x)=sech(kx)` and `w(x)=u f(x)`. Exact one-dimensional integrals give

`||f||²=2/k`, `||f'||²=2k/3`, `int sech²(kx) f² dx=4/(3k)`.

Thus the *exact* Rayleigh quotient is

`R(u) = k²/3 + u^T[V_inf-(2/3)B]u`.

Consequently, if

`lambda_min(V_inf-(2/3)B) < -k²/3`,

then `L_perp` has a strictly negative quadratic-form direction. If `V_inf` is strictly positive definite, the negative spectrum is necessarily discrete, so at least one genuine eigenvalue is negative. This is a sufficient **instability direction**, not a necessary-and-sufficient spectral classification. It bypasses a full spectral-digit campaign.

## C. Exact noncommuting adverse countercontrol

Choose explicitly, **for independent illustration, NOT as FIN arrays**,

```
V_inf = [[1, 1/2], [1/2, 1]],
B     = [[4, 0],   [0,   0]],
A     = [[-3,1/2],[1/2, 1]],    k=v=1.
```

Then `[V_inf,B] = [[0,-2],[2,0]] !=0` and `spec(V_inf)={1/2,3/2}>0`. Nevertheless choose `w=(sech x,0)`. Exact values:

- `||w||²=2`;
- kinetic energy quadratic form `2/3`;
- constant far-vacuum contribution `2`;
- attractive well contribution `-16/3`;
- **total quadratic form `-8/3`, Rayleigh quotient `-4/3`.**

The essential spectrum starts at `1/2`; hence there is a negative discrete eigenvalue. **Positive far-vacuum mass is insufficient.** This is not a demonstration that the unavailable actual FIN arrays yield a negative mode.

## D. Stop condition

No source-bound verdict `STABLE_FIN_WALL` or `UNSTABLE_FIN_WALL` is admissible until `M2`, `K_s`, `Y`, all actual phase blocks, the full wall, sign conventions and operator domain are bound and replayed. The commuting PT seed inherited in the order is treated as a *control*; the new contribution is the noncommuting form theorem and exact negative-direction test, not a relabeling of that seed.
