# F471 Current Strict Pair1 Sign‑Distinction Obstruction From Even Reference Weights Packet (No False‑PASS)

Status: `F471_EXECUTED_CURRENT_STRICT_PAIR1_SIGN_DISTINCTION_OBSTRUCTION_FROM_EVEN_REFERENCE_WEIGHTS_PACKET_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

`H37` asks for a strict-core sign-sensitive state object or observable distinguishing `u` from `-u` on `pair1`.

A tempting strict move is to try to use the already exported strict internal reference shape
`ord_{Z12}` / `r_ord(x) ∝ exp(-alpha_geo*ord_{Z12}(x))` (`F446/N479/N480`) to define a sign-sensitive scalar observable,
e.g. a weighted sum or correlation of the form:

```text
S(u) := Σ_x w(x) u(x)
```

and then claim `S(u)` distinguishes `u` from `-u` by sign.

This packet records the **narrowest honest obstruction** on the current exported strict instance:

- the current exported `pair1` representative is `u_1 = ± s1` (odd under reflection `x↦-x`),
- the `ord_{Z12}` and `r_ord` weights are **even** under reflection,
- therefore the weighted sum cancels exactly: `Σ_x w(x) u_1(x) = 0`.

So the `ord`/`r_ord` reference family cannot by itself supply a sign-distinction observable on `pair1` in this current exported instance.

This does **not** prove that `H37` is impossible in strict core; it only blocks one specific, recurrent route.

## Strict‑admissible inputs reused

1. `F446`
   - exports `ord_z12_by_x` inside `r_ord_z12_v1_reference_distribution`.
2. `F456`
   - exports the current strict `pair1` representative vector `u_1` (and its `(c1,s1)` coordinates).

## Exported artifacts

Executed by:

```text
python3 fundamental_action_reconstruction/f471_current_strict_pair1_sign_distinction_obstruction_from_even_reference_weights_packet.py
```

Exports:

- `fundamental_action_reconstruction/generated/sign_distinction_obstruction_pair1_even_reference_weights_strict_v1.json`
- `fundamental_action_reconstruction/generated/f471_current_strict_pair1_sign_distinction_obstruction_from_even_reference_weights_packet_summary.json`

## Hard limits

This packet does **not** claim:

1. discharge of `H37` (it records an obstruction, not a sign-sensitive physical datum),
2. strict-core selector closure / admissible `S_sel_int`,
3. global discharge of `QW-2191`,
4. ToE closure.
