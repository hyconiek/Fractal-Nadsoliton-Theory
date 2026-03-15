# N506 Current First Strict Pair1/Pair2 Chart-Transport Operator Projector Sign Gauge‑Irrelevance Theorem

Status: `N506_DISCHARGED_CURRENT_FIRST_STRICT_PAIR1_PAIR2_CHART_TRANSPORT_OPERATOR_PROJECTOR_SIGN_GAUGE_IRRELEVANCE_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

`F461` exports an explicit lane-scoped `pair1↔pair2` chart-transport operator `O_12` on the declared `n=12` Fourier carrier.

This theorem records the narrow strict hygiene point needed for downstream use **without smuggling a sign-sensitive orientation**:

```text
projector-level transport under O_12 is invariant under residual Z2 sign flips u -> -u.
```

Therefore `O_12` can be used as an **axis/projector gluing ingredient** between the two declared local charts
without claiming any physical directed-orientation datum.

## Strict-admissible inputs reused

1. `F451/N489`
   - strict sigma-int slot-free theta-pair supply exporting `u_1,u_2` (declared scope).
2. `F457`
   - strict-derived transition angle `alpha_12 := (theta_2-theta_1) mod 2π`.
3. `F461`
   - explicit lane-scoped orthogonal chart-transport operator `O_12` on the `n=12` Fourier carrier.
4. `A10`
   - anti-overclaim boundary.

## Theorem

### Claim 1. Rank‑one projectors are invariant under sign flips.

For any real unit vector `u`, define the rank‑one projector:

```text
P(u) := u u^T.
```

Under the residual sign flip `u' := -u` we have:

```text
P(u') = (-u)(-u)^T = u u^T = P(u).
```

So the projector-level representative is sign-gauge invariant. ∎

### Claim 2. Conjugation transport of projectors is sign-gauge invariant.

Let `O` be any orthogonal operator (`O^T O = I`). Define the transported projector:

```text
P_O(u) := O P(u) O^T.
```

Then for `u'=-u`:

```text
P_O(u') = O P(u') O^T = O P(u) O^T = P_O(u).
```

Therefore, for the exported lane-scoped chart transport `O_12` of `F461`,
the downstream transported projector object `O_12 P(u) O_12^T` is sign-gauge invariant. ∎

## Meaning (B3 continuation without false‑PASS)

This theorem means only:

1. `F461` supplies a **lane-scoped** (pair1/pair2) chart-transport object that is safe to use at the **axis/projector level**,
2. residual `Z2` sign does not carry physical content for such projector-level transport,
3. lifting residual `Z2` to a **sign-sensitive physical orientation datum** remains open (`B3` scope),
4. none of this upgrades to a global selector atlas (`H41`) nor discharges global `QW-2191` (`H40`).

## What N506 does not claim

`N506` does not claim:

1. a sign-sensitive physical orientation datum derived,
2. a global selector atlas or cocycle/gluing data exported,
3. global discharge of `QW-2191`,
4. strict-core selector closure / admissible `S_sel_int`,
5. ToE closure.

