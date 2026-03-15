# N494 Current First Strict `QW-2190` Diagonal/Local Mode‑Index Canonicalization — Uniqueness Up To Conjugation Theorem

Status: `N494_DISCHARGED_CURRENT_FIRST_STRICT_QW2190_DIAGONAL_LOCAL_MODE_INDEX_CANONICALIZATION_UNIQUENESS_UP_TO_CONJUGATION_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

`QW-2191` proves that **kernel alone** leaves a continuous `O(2)` family inside each degenerate Fourier pair plane, so full
mode‑index assignment is not axiom‑free unique from the translation‑invariant host.

The current repo now exports a strict diagonal/local lane which cuts this freedom:

- continuous `O(2)` is reduced to residual `Z2` on all degenerate pairs (`N487`, supported by `N484/N485`),
- and an explicit strict‑derived mode‑index assignment basis object is exported (`F453`, packaged by `N492`).

This theorem records the next necessary clarification:

```text
for the QW‑2190 SU(3)/SU(2) embedding audits, the residual Z2 sign freedom is only a conjugation/basis convention,
so the diagonal/local lane canonicalizes the embedding uniquely up to conjugation in the declared scope.
```

No global selector closure and no ToE closure are implied.

## Strict-admissible evidence reused

1. `QW-2190`
   - mode scaffold + embedded SU(3)/SU(2) audits,
2. `QW-2191`
   - kernel‑alone continuous `O(2)` obstruction (kept true in its scope),
3. `N484/N485/N487`
   - diagonal/local `O(2) -> Z2` cut criterion + current strict‑derived all‑pairs cut decision + scoped discharge of the continuous family,
4. `F453/N492`
   - exported strict‑derived mode‑index assignment basis object and its packaging as a lane‑scoped internal orientation datum (axis‑only),
5. `N493` (+ `P452`)
   - residual `Z2` sign flips are conjugation equivalences for the QW‑2190 SU(3)/SU(2) embedding audits.

## Theorem (canonicalization up to conjugation on the diagonal/local lane)

### Claim 1. Continuous `O(2)` degeneracy is eliminated on the diagonal/local lane (residual `Z2` remains).

By `N487`, on the exported canonical local‑diagonal residual sector for `n=12`, every degenerate pair plane `pair_m (m=1..5)`
has nonzero diagonal/local Fourier defect `|F_{2m}(d)| ≠ 0`, hence the diagonal/local sector cuts the continuous `O(2)`
family down to residual `Z2` (sign) on each pair plane.

Therefore the `QW-2191` continuous family obstruction is discharged **in the declared diagonal/local lane scope**, but only
up to residual sign. ∎

### Claim 2. The repo exports an explicit strict‑derived representative basis on that lane.

By `F453` (packaged in `N492`), the repo exports explicit vectors:

- `B3 := [e0, u_{1,+}, u_{1,-}]` (SU(3) subspace basis),
- `B2 := [u_{2,+}, u_{2,-}]` (SU(2) subspace basis),

constructed from the diagonal/local defect angles `θ_*(m)` on each `pair_m`.

This provides an explicit canonical representative basis point (axis‑fixed) on the declared diagonal/local lane, with only
residual sign freedom remaining. ∎

### Claim 3. The residual `Z2` sign freedom does not change the QW‑2190 embedding audits (conjugation equivalence).

By `N493`, residual sign flips in the exported bases correspond to conjugations of the embedded generators by a lifted
subspace conjugator, and therefore:

- kernel‑subspace invariance audits are unchanged, and
- Lie‑closure audits are unchanged,

in the QW‑2190 SU(3)/SU(2) embedding sense.

Therefore, for the QW‑2190 embedding audits, the diagonal/local lane fixes the embedding uniquely **up to conjugation** in
the declared scope. ∎

## What `N494` does not claim

`N494` does not claim:

1. axiom‑free **global** physical uniqueness beyond the declared diagonal/local lane and QW‑2190 audit scope,
2. global discharge of `QW-2191` (kernel‑alone scope remains true),
3. strict‑core selector closure nor admissible `S_sel_int`,
4. ToE closure.

## Consequence (next honest frontier)

After `N494`, the remaining honest uniqueness frontier is no longer “residual sign inside the diagonal/local lane”.
It is instead:

1. extension beyond the diagonal/local lane without importing an external selector postulate, and/or
2. continuation of strict‑only ToE closure on the remaining non‑uniqueness / non‑closure frontiers explicitly listed by the current strict audits.

