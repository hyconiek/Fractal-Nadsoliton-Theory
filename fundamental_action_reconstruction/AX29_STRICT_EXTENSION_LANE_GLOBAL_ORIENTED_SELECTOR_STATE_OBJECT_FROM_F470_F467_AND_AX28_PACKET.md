# AX29 Strict‑Extension Lane Global Oriented Selector State Object From F470/F467 + AX28 Packet

Status: `AX29_EXECUTED_STRICT_EXTENSION_LANE_GLOBAL_ORIENTED_SELECTOR_STATE_OBJECT_FROM_F470_F467_AND_AX28_PACKET_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

Strict core now exports:

1. a global **projective/ray-level** selector state object on `C_v1` (`F470/N516`), and
2. a lane-scoped **vector-level oriented transport** package on `{pair1..pair5}` tracked explicitly as a gauge/convention layer (`F467/N511`),

but still does **not** export any strict sign-sensitive physical orientation datum (`H37` remains open).

After `AX28`, in `strict_extension_only` scope we now have an explicit **sign-fixing observable** on `pair1` that selects a preferred representative of the ray by requiring `S_dir(u_1)>0`.

`AX29` packages the minimal global consequence (still in extension scope):

> Export an explicit global oriented vector selector state object on `C_v1` by applying the `AX28` sign fix to the existing `F467` vector section, while keeping the underlying strict projective state (`F470`) unchanged.

## Exported artifacts

Executed by:

```text
python3 fundamental_action_reconstruction/ax29_strict_extension_lane_global_oriented_selector_state_object_from_f470_f467_and_ax28_packet.py
```

Exports:

- `fundamental_action_reconstruction/generated/strict_extension_lane_selector_state_global_c_v1_oriented_vector_v1.json`
- `fundamental_action_reconstruction/generated/ax29_strict_extension_lane_global_oriented_selector_state_object_from_f470_f467_and_ax28_packet_summary.json`

Probe-level hygiene:

- `P473` audits that the exported extension-lane oriented vectors reproduce (up to numeric tolerance) the same strict-core projector operators in each chart `(c_m,s_m)` as the strict global projective selector state (`F470`), so the extension-lane step fixes only a sign-gauge representative and does not change the underlying strict ray/projector state.

## Hard limits (no false pass)

`AX29` does **not** claim:

1. strict-core discharge of `H37` (this is extension scope only),
2. strict-core canonical generator/orientation fixing datum export (`T164` remains open),
3. strict-core selector closure / admissible `S_sel_int`,
4. global discharge of `QW-2191`,
5. ToE closure.
