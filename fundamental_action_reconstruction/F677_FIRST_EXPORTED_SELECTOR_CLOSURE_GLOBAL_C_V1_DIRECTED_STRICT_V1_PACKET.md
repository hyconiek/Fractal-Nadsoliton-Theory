# F677 First Exported `SelectorClosure_global_C_v1_directed_strict_v1` Packet

Status: `F677_FIRST_EXPORTED_SELECTOR_CLOSURE_GLOBAL_C_V1_DIRECTED_STRICT_V1_PACKET_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

Export one explicit **global selector closure object** on the strict configuration space `C_v1` in the **directed (vector-level) scope**:

```text
SelectorClosure_global_C_v1_directed_strict_v1
```

This closure object is intentionally **premise-based**: it is built on top of the exported directed global selector state object (`F474/N524`), which itself depends on the explicit fixing datum `T164` (`F473/N523`) and does not claim `Aut(Z_12)`-invariant sign canonicity (`N462` boundary).

## Construction (explicit, vector-level, sign-lift made explicit)

For each chart `pair_m`, take:

1. the exported directed representative vector `u_m` from the global directed state object (`F474`),
2. the exported chartwise output channel matrix `Y_sel(pair_m)` in the fixed output basis `(o_+, o_-)` from the global output operator object (`F660`),
3. the induced output vector in `(o_+,o_-)` coordinates:

```text
v_out_raw(pair_m) := Y_sel(pair_m) · u_m
```

On the current export set, `v_out_raw(pair_m)` is a chartwise **directed** representative which is only defined up to a per-chart sign unless one exports an explicit sign-lift/section choice (audit `P674`, boundary `N675`).

This packet therefore exports an explicit **chartwise sign-lift** `s(pair_m) ∈ {±1}` and defines the corrected chartwise output vectors:

```text
v_out(pair_m) := s(pair_m) · v_out_raw(pair_m),
```

with the sign-lift rule declared explicitly in the exported JSON object (no implicit sign fixing).

The exported global directed closure object packages:

- the chartwise raw and sign-lifted output vectors `v_out_raw(pair_m)` and `v_out(pair_m)`,
- an explicit **well-definedness/gluing certificate** (max-norm equality check across the sign-lifted vectors),
- explicit references to the global atlas/transition objects on `C_v1` and to the global directed state/output operator objects used.

## Hard limits (no false pass)

This packet does **not** claim:

- strict-core selector closure,
- kernel-alone/global `QW-2191` discharge,
- ToE closure,
- any operator-level transition groupoid promotion beyond section/projector level (`N512` boundary),
- any `Aut(Z_12)`-invariant sign canonicity (depends on the explicit fixing datum path via `F474`).

