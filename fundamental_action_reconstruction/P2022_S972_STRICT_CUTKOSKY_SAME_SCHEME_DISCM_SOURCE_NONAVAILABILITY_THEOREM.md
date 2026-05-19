# P2022 S972 Strict Cutkosky Same-Scheme DiscM Source Audit Theorem

Status: `P2022_SAME_SCHEME_DISCM_SOURCE_PARTIAL_TREE_VERTEX_AVAILABLE_WITH_TRACE_NO_FALSE_PASS`
As of: `2026-05-18`

## Goal

P2020 exported the exact tree-level real-linear-polarization CutSum matrix for
`graviton -> gauge_gauge`.  P2021 then corrected the route by rejecting the P1994
`Z(s)` proxy as a P1953 dressed-residue input.

P2022 asks whether the current strict-side repository already exports enough
source data to compute `DiscM_common_basis` for the exact P2020 `{plus,cross}`
matrix normalization.

The corrected answer is more nuanced than the first P2022 draft:

```text
the minimal tree-level hAA source chain is available,
but the loop-derived, same-scheme DiscM source data are not yet available.
```

## Source-status correction

The previous P2022 wording was too strong when it marked D1 as a hard missing
source.  The repository already contains:

1. P1955: minimal tree-level `hAA` vertex from the gauge-metric density,
2. P2019: the corresponding local transverse tree component,
3. P2020: the angularly transported, phase-space integrated `{plus,cross}`
   CutSum matrix.

Therefore D1 is now classified as

```text
PARTIAL_TREE_LEVEL_SOURCE_AVAILABLE_NOT_LOOP_DISCM_READY.
```

This partial discharge does **not** make `DiscM_common_basis` available.  It only
means the obstruction is not the total absence of a tree `hAA` vertex; the real
blockers are loop-level and projector/scheme-level.

## Remaining blocking source objects

P2022 still records four hard failures:

1. gauge fixing, ghost sector, and BRST physical-state projector for the same
   channel,
2. loop-derived dressed pole/residue data, not seed inheritance and not proxy
   `Z(s)`,
3. evaluated `DiscM_common_basis` in the exact P2020 `{plus,cross}` basis,
4. a single renormalization/gauge scheme lock across the vertex, propagator,
   DiscM backend, and CutSum.

This is not a theory no-go result; it is a repository-state theorem saying the
data needed for theorem-grade Cutkosky comparison are still incomplete.

## Symbolic underdetermination witness

P2022 records the current algebraic obstruction in matrix form.  If the still
unexported loop side is represented as

```text
DiscM = (alpha_gr*J_R2 + beta_gr*J_Rmunu2 + gamma_gr*J_EH_mix) I_2,
```

and the P2020 no-identical-symmetry CutSum is

```text
CutSum = (1/pi) I_2,
```

then `DiscM - CutSum` cannot be decided from P2020/P2021 alone.  One assignment
of the unexported loop coefficients gives trace defect `-2/pi`; another gives
trace defect `0`.  Thus the equality is underdetermined until the actual loop
coefficients and master integrals are exported.

## No-false-pass boundary

P2022 is not:

1. a computation of `DiscM_common_basis`,
2. a proof of `DiscM = CutSum`,
3. a BRST cohomology theorem,
4. all-state unitarity,
5. `QW-2191` discharge,
6. ToE closure.

It is a precise blocker map for the next strict-side unitarity move, now with
the existing tree-level `hAA` source chain credited correctly.

## Progress toward ToE

P2022 moves the unitarity route away from proxy-factor manipulation and onto the
actual source-data frontier.  The exact P2020 CutSum matrix and minimal tree
`hAA` source chain are ready.  The left side of the optical theorem is still not
computable until the loop-level vertex/self-energy data, BRST projector, ghost
trace, and scheme lock are exported.

## Next honest step

Build P2023 by deriving or explicitly tabulating the BRST physical-state
projector and ghost-sector exclusion trace for `graviton -> gauge_gauge` in one
fixed scheme.  In parallel, specify the loop-level self-energy/vertex insertion
rules that upgrade the tree `hAA` chain into a same-scheme `DiscM_common_basis`
calculation.
