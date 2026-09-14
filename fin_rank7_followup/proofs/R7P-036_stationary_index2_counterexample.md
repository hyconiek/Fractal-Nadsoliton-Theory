# R7P-036 — certified stationary full-seven-coordinate index-two witness

Status: **INTERVAL_CERTIFIED** for the accepted strict spectral intervals, at the exact supplied gain `g=5`.

## Exact invariant family

Use

\[
h_j=J\cos(\pi j/2)+K(-1)^j.
\]

The softmax probability is period four and reflection even. Summation over the three repeated period-four cells makes every omitted full-`X7` stationarity component vanish identically. Hence a root of the two reduced equations is a genuine stationary point of the full seven-coordinate dual, not merely of `C4`.

Writing

\[
D=\cosh J+e^{-2K},\quad
x={\sinh J\over D},\quad
y={\cosh J-e^{-2K}\over D},
\]

the full stationarity condition on this invariant family is exactly

\[
J=5{\lambda_3\over6}x,\qquad
K=5{\lambda_6\over12}y.
\]

## Root isolation

A parametric interval Krawczyk calculation, using the accepted outward strict intervals for `lambda3` and `lambda6`, proves a unique root in

\[
J\in[1.244948095164464,1.244948295164464],
\]

\[
K\in[0.7795555209881978,0.7795557209881978].
\]

The rationalized preconditioner, interval residual, interval Jacobian, Krawczyk image and strict inclusion margins are serialized in `certificates/R7P-036_stationary_index2_witness.json`. The inclusion is uniform over the accepted spectral boxes; in particular it contains the fixed strict tuple.

## Full H7 inertia without floating eigenvalue decisions

At the stationary root,

\[
x={J\over 5(\lambda_3/6)},\qquad
y={K\over5(\lambda_6/12)}.
\]

Period-four symmetry block-diagonalizes the full seven-coordinate Hessian into:

1. the `3sin` singleton;
2. a `(3cos,6alt)` block;
3. a `(4cos,5cos)` block;
4. a `(4sin,5sin)` block.

The last two blocks have the same determinant

\[
\Delta_{45}=
\left({1\over5}-{\lambda_4\over12}\right)
\left({1\over5}-{\lambda_5\over12}\right)
-{\lambda_4\lambda_5x^2\over144}.
\]

The interval certificate gives `sup Delta45 < -0.0202278506`. Thus each disjoint 2-by-2 block has exactly one negative eigenvalue.

The remaining blocks are strictly positive: the `3sin` entry has lower bound `>0.16711362`; the `(3cos,6alt)` first diagonal is `>0.09563222` and its determinant is `>0.01087251`. Therefore there are no additional negative or zero directions.

Hence, at this genuine stationary point,

\[
\boxed{n_-(H_7)=2,\qquad n_0(H_7)=0.}
\]

By the already-derived stationary primal/dual inertia theorem (R7P-014), the 11-dimensional primal tangent Hessian has the same two negative and zero directions, plus four additional positive directions.

## Scientific consequence

This is a certified counterexample to the **unrestricted** conjecture “every stationary point of the full seven-coordinate rank-seven active-gain landscape has index at most one.” The conjecture is false already at exact `g=5`.

It does **not** refute a separately formulated gain-restricted statement near the local equal-energy event `g_eq≈3.718344898...`. That narrower question remains open and must be named with its gain domain.
