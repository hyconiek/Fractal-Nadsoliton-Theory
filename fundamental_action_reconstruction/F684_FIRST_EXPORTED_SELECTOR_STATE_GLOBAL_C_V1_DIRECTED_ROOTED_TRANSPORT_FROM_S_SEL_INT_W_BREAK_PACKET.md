# F684 First Exported `SelectorState_global_C_v1_directed` Rooted Transport from `S_sel_int` `w_break` Packet (No False‑PASS)

Status: `F684_FIRST_EXPORTED_SELECTOR_STATE_GLOBAL_C_V1_DIRECTED_ROOTED_TRANSPORT_FROM_S_SEL_INT_W_BREAK_PACKET_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

Export an explicit **global directed (sign-sensitive) selector state object** on `C_v1` by constructing
a vector representative section on `{pair1..pair5}` via:

1. fixing sign on `pair1` using the exported strict-core reflection-breaking weight payload `w_break_by_x` (from `F647`),
2. propagating that sign to other charts via the exported rooted transports `O_1m`,
3. keeping the result explicitly scoped as a **tracked convention/section choice** (no physical orientation claim).

## Output object

```text
SelectorState_global_C_v1_directed_rooted_transport_from_S_sel_int_w_break_strict_convention_v1
```

## Hard limits

- Does not claim `Aut(Z_12)`-invariant sign canonicity.
- Does not claim strict-core selector closure / admissible `S_sel_int`.
- Does not claim kernel-alone/global `QW-2191` discharge.
- Does not claim ToE closure.
