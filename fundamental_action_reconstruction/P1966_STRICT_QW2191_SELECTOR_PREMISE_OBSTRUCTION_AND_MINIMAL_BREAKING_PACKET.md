# P1966 Strict QW-2191 Selector Premise Obstruction and Minimal Breaking Packet

Status: `QW2191_MINIMAL_SELECTOR_PREMISE_IDENTIFIED__STRICT_AXIS_ONLY_SOURCE_FOUND__GLOBAL_SELECTOR_CLOSURE_NOT_PASSED`  
As of: `2026-05-17`

## Goal

Execute the next honest strict-lane step for the `QW-2191` selector block after
`P1965`, without importing the retired legacy kernel and without promoting a
local or scoped selector result into global strict-core closure.

The target is not to claim ToE closure.  The target is to identify, with a
machine-checkable symbolic calculation, the exact mathematical shape of the
selector premise/source needed to bypass the `O(2)` uniqueness obstruction,
and to distinguish the already-exported lane-scoped axis sources from global
sign-sensitive selector closure.

## Inputs respected

This packet respects the current guardrails:

1. `K_strict_gate(d) = cos(omega*d + phi)/(1 + beta*d^eta)` remains the only
   strict working kernel lane used for the surrounding chain.
2. No legacy physical-role claim is imported.
3. The canonical ontology remains
   `nadsoliton -> light -> matter -> emergent observer`.
4. `QW-2191` is treated as a real strict-core obstruction unless a genuine
   selector source is exported.
5. `N487` is not promoted beyond its scoped diagonal/local lane.
6. `AX16` remains `strict_extension_only`, not a strict-core theorem.

## Symbolic obstruction theorem

Let `V_m = span(c_m, s_m)` be one Fourier-degenerate pair plane.  Let
`R(t)` be the standard `O(2)` rotation acting on this plane.  If a symmetric kernel-only / pure `O(2)`-equivariant operator `A` is equivariant under every rotation, then

```text
A R(t) = R(t) A  for every t.
```

The script computes this condition for a generic real matrix

```text
A = [[x, y], [z, w]].
```

Differentiating the commutator at the identity gives the linear constraints
forcing the equivariant symmetric normal form

```text
A = [[x, 0], [0, x]].
```

Therefore the pair-plane eigenvalue has multiplicity two.  A scalar operator
on `V_m` cannot select a unique axis.  This is the precise kernel-only `QW-2191` obstruction restated as an executable symbolic result.  It does not deny the existence of additional strict lane-scoped symmetry-breaking sources; it says they cannot be obtained from an `O(2)`-equivariant scalar pair-plane operator alone.

## Minimal symmetry-breaking premise

The minimal local datum that can cut the continuous `O(2)` family is a nonzero
traceless symmetric selector tensor on each affected pair plane:

```text
Delta_sel_pair_m = [[sigma_c_m, sigma_s_m],
                    [sigma_s_m, -sigma_c_m]].
```

Its characteristic polynomial is

```text
lambda^2 - sigma_c_m^2 - sigma_s_m^2.
```

Thus the spectral gap is nonzero exactly when

```text
sigma_c_m^2 + sigma_s_m^2 > 0.
```

Under that condition the continuous axis-choice family is reduced to the
residual `Z2` sign ambiguity of the selected eigenvector.  This is the
mathematical mechanism required to bypass `QW-2191`.

## Current-state audit verdict after repo grep

The repository does export strict, lane-scoped **axis-only** sources that were
too coarsely classified in the first P1966 audit:

1. the Shannon element-order mode-index assignment (`F454`, with theorem
   support in `N480/N488/N496/N514`),
2. the diagonal/local mode-index assignment (`F453/N487/N492`), and
3. the `P455/N497` alignment audit showing that the two axes agree up to
   residual `Z2`.

However, the global/sign-sensitive selector closure is still not exported:

1. `P1499` still records missing global closure requirements,
2. `P1552` keeps theorem-level uniqueness unproved for its candidate,
3. `N487` remains scoped and does not claim global selector closure, and
4. `AX16` remains `strict_extension_only`, not admissible strict-core
   `S_sel_int`.

Therefore corrected `P1966` identifies the minimal tensor shape, recognizes the
existing strict axis-only source, and still refuses global `QW-2191` closure.

## Machine-checkable outputs

- Script:
  `p1966_s916_strict_qw2191_selector_premise_obstruction_and_minimal_breaking.py`
- Witness JSON:
  `generated/p1966_s916_strict_qw2191_selector_premise_obstruction_and_minimal_breaking.json`

The JSON contains:

1. the symbolic `O(2)` commutant calculation,
2. the minimal traceless symmetric breaking tensor,
3. numerical eigen-gap probes using `numpy`,
4. the current-source audit against `P1499`, `P1552`, `N487`, `AX16`, and
   `P1965`,
5. an explicit false-pass guard distinguishing axis-only `O(2)->Z2` from global selector closure.

## Theorem export

### Negative strict-core theorem

On a Fourier-degenerate pair plane, every symmetric `O(2)`-equivariant kernel-only operator is scalar.  Hence the pair-plane `O(2)`-symmetric layer alone cannot choose a unique selector axis on that plane.

### Conditional extension theorem

If each affected degenerate pair plane is equipped with a nonzero traceless
symmetric tensor `Delta_sel_pair_m` satisfying
`sigma_c_m^2 + sigma_s_m^2 > 0`, then the pair-plane operator has two distinct
eigenvalues and the continuous `O(2)` selector ambiguity is reduced to residual
`Z2` sign choice.

### Axis-source refinement

Repo grep shows that the Shannon element-order strict lane already supplies an
axis-only internal source.  `P1967` maps that source into this `Delta_sel`
tensor form for `pair1..pair5` and verifies nonzero gaps.

### Non-claim

This is not a global `QW-2191` closure theorem.  Axis-only `O(2)->Z2` is not
the same as sign-sensitive physical orientation, admissible `S_sel_int`, or
full ToE closure.

## Next honest step

Use `P1967` as the constructive axis source, then test whether each downstream
observable needed for strict ToE closure is insensitive to the residual `Z2`
sign.  If not, construct or refute a strict sign-sensitive datum.

## Lay explanation

If the mathematics is perfectly round inside a two-dimensional plane, it cannot itself say which direction is special.  Corrected `P1966` proves that statement and also notes that the repo already has a strict lane-scoped compass for the axis.  What remains open is stronger: choosing a signed arrow/global selector, not merely an unsigned axis.
