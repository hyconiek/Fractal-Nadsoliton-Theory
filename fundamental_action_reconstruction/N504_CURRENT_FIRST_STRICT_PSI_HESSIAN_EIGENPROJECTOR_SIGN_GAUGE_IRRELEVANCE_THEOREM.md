# N504 Current First Strict `Psi` Hessian Eigenprojector Sign Gauge‑Irrelevance Theorem

Status: `N504_DISCHARGED_CURRENT_FIRST_STRICT_PSI_HESSIAN_EIGENPROJECTOR_SIGN_GAUGE_IRRELEVANCE_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

`F459` exports a strict-derived **numeric** eigensystem for the full diagonal/local `Psi` Hessian `H_psi`, including a concrete eigenvector
representative matrix.

Any eigenvector representative necessarily carries a residual sign freedom:

```text
v  ->  -v.
```

This theorem closes the narrow “no false‑PASS” hygiene point needed for safe downstream use:

```text
if downstream objects are expressed projectively (rank‑one eigenprojectors) or via squared overlaps,
then residual eigenvector sign is gauge and cannot change the exported downstream object.
```

This theorem does **not** derive a sign‑sensitive physical orientation datum, does **not** discharge global `QW-2191`, and does **not** claim ToE closure.

## Strict-admissible inputs reused

1. `F459`
   - strict-derived diagonal/local `Psi` Hessian eigensystem value instantiation (eigenvalues + eigenvector representatives).
2. `F460`
   - strict-derived export of the rank-one spectral projectors `P_j := |v_j><v_j|` (sign-gauge invariant packaging of `F459`).
3. `P463`
   - the eigenmode projection audit uses squared overlaps `|<b_i, v_j>|^2`, which is sign-invariant.
4. `A10`
   - anti-overclaim boundary.

## Theorem

### Claim 1. Eigenprojectors are invariant under eigenvector sign flips.

Let `v` be any real column vector and define the rank‑one projector:

```text
P := |v><v| = v v^T.
```

Under the residual sign flip `v' := -v`:

```text
P' = v' (v')^T = (-v)(-v)^T = v v^T = P.
```

Therefore the eigenprojector object is sign‑gauge‑invariant. ∎

### Claim 2. Squared overlaps are invariant under sign flips.

For any real vectors `b` and `v` define the scalar overlap `⟨b, v⟩` and its squared magnitude:

```text
q(b,v) := |⟨b, v⟩|^2.
```

Under either sign flip `v -> -v` or `b -> -b`:

```text
|⟨b, -v⟩|^2 = |-⟨b, v⟩|^2 = |⟨b, v⟩|^2,
|⟨-b, v⟩|^2 = |-⟨b, v⟩|^2 = |⟨b, v⟩|^2.
```

Therefore any mixing/intensity matrix built from squared overlaps is sign‑gauge‑invariant. ∎

## Meaning (no false‑PASS)

This theorem means only:

1. for the diagonal/local `Psi` Hessian eigensystem exported by `F459`, residual eigenvector sign conventions cannot affect downstream objects
   expressed projectively as eigenprojectors (as exported by `F460`),
2. the `P463` eigenmode projection audit outputs based on squared overlaps are sign‑gauge‑invariant,
3. therefore residual sign is a tracked gauge/convention layer for these eigenmode/eigenprojector downstream objects; lifting it to a
   sign‑sensitive physical orientation datum remains a separate open problem (cf. `B3` scope) and must not be smuggled in.

## What N504 does not claim

`N504` does not claim:

1. a sign‑sensitive physical orientation datum derived,
2. strict‑core selector closure / admissible `S_sel_int`,
3. global discharge of `QW-2191`,
4. ToE closure.

