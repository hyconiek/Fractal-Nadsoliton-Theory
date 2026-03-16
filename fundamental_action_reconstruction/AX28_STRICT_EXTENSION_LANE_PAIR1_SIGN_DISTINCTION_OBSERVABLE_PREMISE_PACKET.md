# AX28 Strict‑Extension Lane Pair1 Sign‑Distinction Observable Premise Packet

Status: `AX28_EXECUTED_STRICT_EXTENSION_LANE_PAIR1_SIGN_DISTINCTION_OBSERVABLE_PREMISE_PACKET_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

`H37` remains open in strict core: the repo exports no **strict** sign‑sensitive physical state object or observable distinguishing `u` from `-u` on `pair1`.

After `N518` and the probe hygiene `P472`, the honest strict conclusion is:

- any direction‑free (`Aut(Z_12)`‑invariant) reference weights are even and cannot distinguish sign via a scalar `Σ_x w(x)u(x)`,
- no strict(-derived) **weight-like** reflection‑breaking per‑site reference is currently exported beyond non‑canonical marked‑site `K_total` rows,
- therefore a strict discharge of `H37` requires an explicit reflection‑breaking / orientation source (e.g. a generator/orientation fixing datum `T164`).

This packet records the minimal forward step **without false pass**:

> In `strict_extension_only` scope, explicitly fix a `Z_12` generator/orientation, and then export a sign‑distinction observable on `pair1` of the form `S(u)=Σ_x w_dir(x)u(x)` using a directed reference weight `w_dir`.

This does **not** change strict core and does **not** claim `T164` is discharged.

## Extension-scope premises (explicit)

Accepted only in:

```text
strict_extension_only
```

1. Fix a preferred `Z_12` generator/orientation (marked direction):
   - successor map `suc_fix(k) := (k+1) mod 12` (generator `g_fix := 1`).
2. Use the directed coordinate `x ∈ {0,...,11}` relative to the identity site `0`.

## Observable definition (sign-sensitive scalar)

Let `alpha_geo := 4 ln 2` (strict-side numeric source exists as `alpha_geo_strict_derived_v1`).

Define a directed exponential weight (normalization irrelevant for sign):

```text
w_dir(x) := exp(-alpha_geo * x).
```

For any `pair1` representative vector `u(x)` in the site basis, define:

```text
S_dir(u) := Σ_{x=0..11} w_dir(x) u(x).
```

Then `S_dir(-u) = -S_dir(u)` and (for typical non-even weights) `S_dir(u) ≠ 0`, so `S_dir` distinguishes the sign of `u`.

## Exported artifacts

Executed by:

```text
python3 fundamental_action_reconstruction/ax28_strict_extension_lane_pair1_sign_distinction_observable_packet.py
```

Exports:

- `fundamental_action_reconstruction/generated/strict_extension_lane_pair1_sign_distinction_observable_packet.json`
- `fundamental_action_reconstruction/generated/ax28_strict_extension_lane_pair1_sign_distinction_observable_packet_summary.json`

## Hard limits (no false pass)

`AX28` does **not** claim:

1. strict-core discharge of `H37` (this is extension scope only),
2. strict-core canonical generator/orientation fixing datum export (`T164` remains open),
3. strict-core selector closure / admissible `S_sel_int`,
4. global discharge of `QW-2191`,
5. ToE closure.

