# H4 RESIDUAL O(2) REDUCTION OF LIGHT-FEEDBACK ANSATZ

Status: `PASS_PARTIAL_REDUCED_SELECTOR_TEST_PACKET_READY`
As of: `2026-03-06`

## Goal

Zredukowac ansatz `H3` do residualnego sektora selektora `O(2)` i zapisac najprostszy test:
- czy `K_obs` daje efekt trywialny,
- czy daje rzeczywisty symmetry-breaking,
- czy pozostaje nierozstrzygalny bez jawnych wspolczynnikow bloku zredukowanego.

## Setup

For each current mode pair `i`, let
- `V_i = span{c_i, s_i}`
- `u_i(theta_i) = cos(theta_i) c_i + sin(theta_i) s_i`

Let the residually projected operator induced by `H3` on `V_i` be

`A_i := P_i E_i^* G_light^* R_mat^* O_obs R_mat G_light E_i P_i`

where `P_i` is the projection to `V_i`.

The reduced selector energy is then

`q_i(theta_i) := lambda_obs <u_i(theta_i), A_i u_i(theta_i)>`.

## 2x2 reduction formula

In the orthonormal basis `{c_i, s_i}`, write

`A_i = [[a_i, b_i], [b_i, d_i]]`

with real symmetric coefficients.

Then

`q_i(theta_i) = lambda_obs (a_i cos^2 theta_i + 2 b_i sin theta_i cos theta_i + d_i sin^2 theta_i)`

which is equivalently

`q_i(theta_i) = lambda_obs ( (a_i + d_i)/2 + (a_i - d_i)/2 cos(2 theta_i) + b_i sin(2 theta_i) )`.

## Immediate consequence

This reduces the selector question to one sharp dichotomy:

1. **Isotropic case**
   - if `a_i = d_i` and `b_i = 0`, then `q_i(theta_i)` is constant,
   - so `K_obs` does not break the residual `O(2)` degeneracy on that pair.

2. **Anisotropic case**
   - if `(a_i - d_i, b_i) != (0, 0)`, then `q_i(theta_i)` is nonconstant,
   - so `K_obs` induces an orientation-selecting effective energy on that pair.

Therefore the first real test of `H3` is not a global theorem but a coefficient-level question:
- is the projected `2 x 2` block isotropic or anisotropic?

## What is gained

`H4` shows that the future test of `K_obs` is mathematically concrete.
The operator no longer needs to be discussed as a narrative loop; it reduces to a check of the projected `2 x 2` anisotropy data.

## What remains open

`H4` does **not** compute actual `a_i, b_i, d_i`.
So it does not yet show whether the present ansatz is:
- symmetry-breaking,
- symmetry-preserving,
- or zero on the residual selector sector.

## Best current conclusion

The repository now has a packet-ready residual-sector reduction of the light-feedback ansatz. The next decisive question is whether the reduced `2 x 2` block is isotropic or anisotropic.

## Frontier after H4

- `H4_B1 := no explicit computed projected 2x2 coefficient block A_i = [[a_i,b_i],[b_i,d_i]] has yet been extracted for any actual mode pair`
- `H3_B1 := no tested reduction of the H3 ansatz onto the residual O(2) selector sector has yet been constructed` is reduced to coefficient extraction level
- `T12_B1 := strict-core typing judgment with totality and uniqueness remains undischarged`
- `C32_B2 := raw cross-pair overlap route remains degenerate`

## Anti-overclaim boundary

This reduction does **not** show that:
- the ansatz is physically correct,
- the ansatz breaks `O(2)` in the actual theory,
- actual `theta_1`, `theta_2` are exported,
- strict core is closed,
- theorem-level or full-closure PASS has been achieved.
