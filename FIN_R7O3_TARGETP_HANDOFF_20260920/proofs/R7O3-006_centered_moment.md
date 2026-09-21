# R7O3-006 — fixed-witness centered-moment certificate

Let `M4 = Cov(F)` for the four normalized C4 feature coordinates under the exact shared-field distribution, and let `tau=67/250`. Fix an exact rational `4x3` matrix `B` of column rank three and define `Z=B^T F`. Fix any exact rational vector `c in Q^3`.

The identity

`E[(Z-c)(Z-c)^T] - Cov(Z) = (E[Z]-c)(E[Z]-c)^T`

shows in Loewner order that

`B^T M4 B = Cov(Z) <= E[(Z-c)(Z-c)^T]`.

Hence, if the interval proof certifies on an entire parameter cell

`K = tau B^T B - E[(Z-c)(Z-c)^T] > 0`,

then for every nonzero `u in R^3`, with `v=Bu != 0`,

`v^T (tau I4-M4) v = u^T B^T(tau I4-M4)B u >= u^T K u > 0`.

Thus `M4` has a three-dimensional subspace `range(B)` on which every Rayleigh quotient is strictly below `tau`. Courant–Fischer gives

`lambda_2(M4) = min_{dim S=3} max_{0 != v in S} R_M4(v) < tau`,

and therefore in particular `lambda_2(M4) <= 67/250`.

The rational basis is not required to be orthonormal. Its exact Gram matrix `B^T B` is used in `K`, and exact nonzero `3x3` row minors establish rank three. The proposal routine may use floating eigenvectors only to choose a rational `B`; the proof checker receives fixed rational `B,c` and performs no eigensolver or optimization.

For the production chart `A=sqrt(r), u=1-s, v=1-t, y=exp(-2J6)`, the seven aggregate weights are

`1`, `2 s t^3`, `A^2 t^4`, `2 A^2 s t`, `2 A t^2 y`, `2 A s t^(2-sqrt(3)) y`, `2 A s t^(2+sqrt(3)) y`.

All are nonnegative on the compact hull and the first weight is exactly one, so the common normalization denominator is at least one. The proof keeps the shared variables through interval jets; it does not replace the seven weights by an independent weight box.

This lemma proves only the four-amplitude practical threshold `tau=67/250`. It does not prove the sharper sigma target, full-X7 globality, a sourced physical gain, or any selector/physical-role claim.
