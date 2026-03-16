# F472 Current Strict Pair1 Sign‑Distinction Obstruction From Aut‑Invariant Reference Weights Packet (No False‑PASS)

Status: `F472_EXECUTED_CURRENT_STRICT_PAIR1_SIGN_DISTINCTION_OBSTRUCTION_FROM_AUT_INVARIANT_REFERENCE_WEIGHTS_PACKET_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

`H37` asks for a strict-core sign-sensitive state object or observable distinguishing `u` from `-u` on `pair1`.

After `F471/N517`, we have a narrow obstruction showing that the **particular** even reference family `ord_{Z12}` / `r_ord`
cannot distinguish the sign on the current exported `pair1` sine axis via a scalar observable of the form:

```text
S(u) := Σ_x w(x) u(x).
```

This packet records a **strictly stronger** obstruction on the same exported strict instance:

1. any `Aut(Z_12)`‑invariant reference weight family `w : Z_12 -> R` is even under reflection `x↦-x` because `-1 ∈ Aut(Z_12)`,
2. the current exported `u_1` is odd under reflection (it is `± s1`),
3. therefore the weighted sum cancels for **every** direction‑free (`Aut(Z_12)`‑invariant) weight:
   `Σ_x w(x) u_1(x) = 0`.

So no such weight family can supply a sign-distinction observable on `pair1` of the form `Σ_x w(x) u_1(x)` on the current exported instance.

This still does **not** prove that `H37` is impossible in strict core: it only blocks the entire class of
`Aut(Z_12)`‑invariant weight‑based scalar attempts.

## Strict‑admissible inputs reused

1. `F456`
   - exports the current strict `pair1` representative vector `u_1` (and its `(c1,s1)` coordinates).
2. `F446/N479` (optional cross-check only)
   - `ord_{Z12}` (and `r_ord`) is `Aut(Z_12)`‑invariant and hence direction‑free.

## Exported artifacts

Executed by:

```text
python3 fundamental_action_reconstruction/f472_current_strict_pair1_sign_distinction_obstruction_from_aut_invariant_reference_weights_packet.py
```

Exports:

- `fundamental_action_reconstruction/generated/sign_distinction_obstruction_pair1_aut_invariant_reference_weights_strict_v1.json`
- `fundamental_action_reconstruction/generated/f472_current_strict_pair1_sign_distinction_obstruction_from_aut_invariant_reference_weights_packet_summary.json`

## Hard limits

This packet does **not** claim:

1. discharge of `H37` (it records an obstruction, not a sign-sensitive physical datum),
2. strict-core selector closure / admissible `S_sel_int`,
3. global discharge of `QW-2191`,
4. ToE closure.

