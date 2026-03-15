# N493 Current First Strict `QW-2191` Residual `Z2` Sign‑Flip Gauge‑Equivalence (Conjugation) Theorem

Status: `N493_DISCHARGED_CURRENT_FIRST_STRICT_QW2191_RESIDUAL_Z2_SIGN_FLIP_GAUGE_EQUIVALENCE_CONJUGATION_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

After the canonical local‑diagonal `O(2) -> Z2` cut (`N484/N485/N487`) and the exported lane‑scoped internal orientation datum
(`F453/N492`), the remaining ambiguity on each degenerate pair plane is residual `Z2`:

```text
u  ->  -u
```

This theorem records the narrow strict statement needed to avoid a false next blocker:

```text
for the QW‑2190 SU(3)/SU(2) embedding audits, residual Z2 sign flips are pure basis conjugations and do not change
kernel‑subspace invariance nor Lie‑closure audits.
```

This does **not** claim global physical uniqueness, does **not** discharge `QW-2191` in kernel‑alone scope, and does **not** claim ToE closure.

## Strict-admissible evidence reused

1. `QW-2190`
   - mode scaffold and the embedded SU(3)/SU(2) audit definition,
2. `QW-2191`
   - kernel‑alone continuous `O(2)` obstruction (kept true in its scope),
3. `N484/N485/N487`
   - diagonal/local lane cuts continuous `O(2)` down to residual `Z2` on the degenerate pair planes,
4. `F453/N492`
   - exports an explicit strict‑derived eigenbasis representative (up to sign) on each pair plane,
5. `P452`
   - executable audit artifact confirming the conjugation equivalence numerically on the exported `F453` basis (no promotion beyond this theorem’s scope).

## Theorem (residual sign flips are conjugation equivalences for QW‑2190 audits)

Let `B ∈ R^{n×k}` have orthonormal columns (`B^T B = I_k`) and let `P := B B^T` be the orthogonal projector onto `span(B)`.

Let `S ∈ R^{k×k}` be a diagonal sign matrix with entries in `{+1,-1}` and define the sign‑flipped basis:

```text
B' := B S.
```

### Claim 1. Projectors (and thus kernel‑subspace invariance audits) are invariant under residual sign flips.

Since `S^T S = I_k`:

```text
P' := B' (B')^T = B S S^T B^T = B B^T = P.
```

Therefore any invariance audit depending only on the projector (e.g. `||(I-P) K P||`) is unchanged under the residual `Z2` sign flips. ∎

### Claim 2. Embedded generators are conjugate under a lifted subspace conjugator.

Define the embedded generators (as in `QW-2190`) from fixed internal generators `g_a` by:

```text
G_a(B)  := B g_a B^T,
G_a(B') := B' g_a (B')^T.
```

Define the lifted conjugator on the full space:

```text
U := B S B^T + (I - P).
```

Then `U` acts as `S` on `span(B)` and as identity on its orthogonal complement, and:

```text
U G_a(B) U^T
  = (B S B^T + (I-P)) (B g_a B^T) (B S B^T + (I-P))^T
  = B S g_a S B^T
  = G_a(B').
```

Therefore the sign‑flipped embedded generators are conjugate to the baseline embedded generators.
Conjugation preserves commutators, so Lie‑closure relations (and their residual audits) are unchanged under residual `Z2` sign flips. ∎

### Claim 3. Application to the current diagonal/local lane exports.

On the current repo state, `F453` exports orthonormal bases on the `QW-2190` carrier, in particular:

- `B3 := [e0, u_{1,+}, u_{1,-}]` for the SU(3) subspace, and
- `B2 := [u_{2,+}, u_{2,-}]` for the SU(2) subspace,

with residual `Z2` sign freedom on the `u_{m,±}` vectors.

By Claims 1–2, all such residual sign flips yield conjugate embeddings and unchanged invariance/closure audits in the `QW-2190` sense.

`P452` confirms this numerically for the current strict‑derived exported `F453` basis (max conjugation residuals `~1e-15`). ∎

## What `N493` does not claim

`N493` does not claim:

1. axiom‑free **global** physical uniqueness beyond the QW‑2190 embedding audits (e.g. beyond the diagonal/local lane and `n=12` scope),
2. global `QW-2191` discharge (kernel‑alone obstruction remains true),
3. strict‑core selector closure or admissible `S_sel_int`,
4. ToE closure.

