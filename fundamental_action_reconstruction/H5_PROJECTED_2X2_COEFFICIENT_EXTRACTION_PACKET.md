# H5 PROJECTED 2X2 COEFFICIENT EXTRACTION PACKET

Status: `PASS_PARTIAL_PROJECTED_2X2_EXTRACTION_PACKET_READY`
As of: `2026-03-06`

## Goal

Zapisac jawny packet ekstrakcji wspolczynnikow `a_i, b_i, d_i` dla zredukowanego bloku `2x2`
pojawiajacego sie w `H4`, tak aby nastepny krok byl juz tylko obliczeniem konkretnego projected bloku dla wybranej pary modow.

## Starting point

From `H4`, for each current mode pair `i`:

- `V_i = span{c_i, s_i}`
- `A_i = P_i E_i^* G_light^* R_mat^* O_obs R_mat G_light E_i P_i`
- `A_i = [[a_i, b_i], [b_i, d_i]]` in the basis `{c_i, s_i}`

## Coefficient extraction formulas

The coefficients are fully determined by basis-projected matrix elements:

- `a_i = lambda_obs <c_i, A_i c_i>`
- `b_i = lambda_obs <c_i, A_i s_i>`
- `d_i = lambda_obs <s_i, A_i s_i>`

Equivalently, because `A_i` is intended as a real symmetric projected operator on `V_i`,

- `b_i = lambda_obs <s_i, A_i c_i>` as a consistency check.

## Minimal extraction packet

A future actual extraction requires exactly the following ingredients:

1. an actual mode pair label `i`,
2. explicit basis vectors `c_i, s_i`,
3. explicit projected carrier map `P_i`,
4. explicit action of `E_i`,
5. explicit action of `G_light`,
6. explicit action of `R_mat`,
7. explicit action of `O_obs`,
8. the three scalar contractions listed above.

## Immediate reduction of the problem

`H5` reduces the next decisive task to one finite target:

> compute or formally export the three scalars `a_i, b_i, d_i` for at least one actual current mode pair.

No additional selector narrative is needed at this stage.

## What remains open

`H5` does **not** yet provide:
- any actual evaluated `a_i, b_i, d_i`,
- any nonzero off-diagonal entry,
- any proof of anisotropy,
- any proof of isotropy,
- any selector export.

## Best current conclusion

The repository now has a packet-ready finite extraction target for the residual `2x2` test of the light-feedback ansatz. The next strict question is whether any actual current mode pair yields a computable anisotropic block.

## Frontier after H5

- `H5_B1 := no explicit extracted triple (a_i, b_i, d_i) has yet been computed or exported for any actual current mode pair`
- `H4_B1 := no explicit computed projected 2x2 coefficient block has yet been extracted` is reduced to scalar extraction level
- `T12_B1 := strict-core typing judgment with totality and uniqueness remains undischarged`
- `C32_B2 := raw cross-pair overlap route remains degenerate`

## Anti-overclaim boundary

This packet does **not** show that:
- any actual pair has anisotropy,
- the ansatz works,
- the ansatz fails,
- `QW-2191` is discharged,
- theorem-level or full-closure PASS has been achieved.
