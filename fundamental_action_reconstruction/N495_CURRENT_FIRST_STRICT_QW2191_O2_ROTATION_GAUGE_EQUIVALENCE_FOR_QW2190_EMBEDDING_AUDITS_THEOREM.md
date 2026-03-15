# N495 Current First Strict `QW-2191` O(2) Rotation Gauge‑Equivalence for `QW-2190` Embedding Audits Theorem

Status: `N495_DISCHARGED_CURRENT_FIRST_STRICT_QW2191_O2_ROTATION_GAUGE_EQUIVALENCE_FOR_QW2190_EMBEDDING_AUDITS_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

`QW-2191` states a strict uniqueness obstruction:

```text
kernel-alone translation-invariant data does not canonically select an axis inside each degenerate Fourier pair plane;
there is a continuous O(2) family of mode-index assignments.
```

Separately, `QW-2190` defines SU(3)/SU(2) embedding audits by:

1. selecting an orthonormal basis `B` of a declared subspace, and then
2. embedding fixed internal generators `g_a` by:

```text
G_a(B) := B g_a B^T.
```

This theorem records one strict clarification needed to avoid a false next blocker:

```text
for the QW-2190 embedding audits, continuous O(2) basis rotations inside the chosen subspace act only as conjugations
and do not change invariance or Lie-closure audits; therefore the O(2) family is gauge for those audits.
```

This does **not** discharge `QW-2191` as a canonical-representative obstruction and does not claim ToE closure.

## Strict-admissible evidence reused

1. `QW-2190`
   - embedded generator definition and invariance / Lie‑closure audit definitions,
2. `QW-2191`
   - kernel-alone continuous `O(2)` mode-index assignment family (kept true in its scope),
3. `N493`
   - residual `Z2` sign-flip conjugation argument (used here only as a special case template),
4. `P454`
   - executable audit artifact confirming the continuous-rotation statement numerically on a tested grid (no promotion beyond this theorem’s scope).

## Theorem (orthogonal basis changes are conjugations for the QW‑2190 embedding audits)

Let `B ∈ R^{n×k}` have orthonormal columns (`B^T B = I_k`) and let `P := B B^T` be the orthogonal projector onto `span(B)`.

Let `R ∈ O(k)` be any orthogonal matrix (`R^T R = I_k`) and define a new orthonormal basis:

```text
B' := B R.
```

### Claim 1. Projectors (and thus kernel-subspace invariance audits) are invariant under orthogonal basis changes.

Since `R R^T = I_k`:

```text
P' := B' (B')^T = B R R^T B^T = B B^T = P.
```

Therefore any invariance audit depending only on `P` (e.g. `||(I-P) K P||`) is unchanged by the full `O(k)` basis freedom. ∎

### Claim 2. Embedded generators are conjugate under a lifted subspace conjugator.

Fix internal generators `g_a ∈ C^{k×k}` and define:

```text
G_a(B)  := B g_a B^T,
G_a(B') := B' g_a (B')^T.
```

Define the lifted conjugator on the full space:

```text
U := B R B^T + (I - P).
```

Then `U` acts as `R` on `span(B)` and as identity on its orthogonal complement, and:

```text
U G_a(B) U^T
  = (B R B^T + (I-P)) (B g_a B^T) (B R B^T + (I-P))^T
  = B R g_a R^T B^T
  = G_a(B').
```

Therefore the embedded generator families obtained from `B` and `B'` are conjugate.
Conjugation preserves commutators, so Lie‑closure audits are unchanged under the full `O(k)` basis freedom. ∎

### Claim 3. Application to the `QW-2191` `O(2)` family on QW‑2190 degenerate pair planes.

On the `QW-2190` carrier, the `QW-2191` obstruction is precisely an `O(2)` freedom inside degenerate pair planes.
For the `QW-2190` embedding audits, those planes enter only through chosen orthonormal bases.

By Claims 1–2, replacing any such basis by a rotated basis (`O(2)` in a pair plane, or a block-diagonal `O(k)` in
multi-column subspaces) does not change:

1. kernel-subspace invariance audits (projector invariance), nor
2. Lie‑closure audits (conjugation equivalence).

`P454` confirms this numerically for a grid of continuous rotations applied to the current exported basis objects.
Therefore, for the **QW‑2190 embedding audits**, the continuous `O(2)` family from `QW-2191` is a gauge/basis convention. ∎

## What `N495` does not claim

`N495` does not claim:

1. axiom‑free global physical uniqueness of a canonical representative basis in ambient coordinates,
2. global discharge of `QW-2191` as a kernel-alone canonical-representative obstruction,
3. strict-core selector closure or admissible `S_sel_int`,
4. ToE closure.

