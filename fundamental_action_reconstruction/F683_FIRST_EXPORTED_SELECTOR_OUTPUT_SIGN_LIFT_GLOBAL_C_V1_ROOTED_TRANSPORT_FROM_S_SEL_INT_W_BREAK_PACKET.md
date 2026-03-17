# F683 First Exported Global Selector Output Sign-Lift (Rooted Transport from `S_sel_int` `w_break`) Packet

Status: `F683_FIRST_EXPORTED_SELECTOR_OUTPUT_SIGN_LIFT_GLOBAL_C_V1_ROOTED_TRANSPORT_FROM_S_SEL_INT_W_BREAK_PACKET_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

`N675` records a strict boundary: raw **directed** closure outputs are obstructed without an explicit sign-lift / section choice.

`P683` audits that, once a root sign is fixed on `pair1` using the exported reflection-breaking seed weight `w_break` (`F647`),
one can propagate a consistent sign choice to `{pair2..pair5}` by aligning the exported rooted transports `O_1m` with the exported representative vectors.
This yields chartwise output vectors with consistent `o_+` sign under the exported promoted output channels `Y_sel(pair_m)`.

`F683` performs the narrowest honest export below any physical orientation claim:

1. export one explicit **global** sign-lift/section-choice object on `C_v1`,
2. with its rooted rule and its dependency on axis-only transport representatives kept explicit,
3. intended to be usable downstream wherever directed sign is only a tracked gauge/convention.

This does **not** claim strict-core selector closure, kernel-alone/global `QW-2191` discharge, or ToE closure.

## Construction (rooted transport sign lift; explicit)

1. Fix root sign on `pair1` by:

```text
s_1 := sign( Σ_x w_break(x) u_1(x) )
```

2. For each `m=2..5`, propagate sign by requiring alignment with the rooted transport:

```text
v_m := O_1m (s_1 u_1)
s_m := sign( <v_m, u_m> )
```

3. Apply `s_m` to the exported chart coordinates and push through `Y_sel(pair_m)` to form the chartwise output vectors in `(o_+,o_-)`.

## Outputs

`F683` exports:

1. global rooted sign-lift object:
   - `fundamental_action_reconstruction/generated/selector_output_sign_lift_global_c_v1_rooted_transport_from_s_sel_int_w_break_strict_convention_v1.json`
2. packet summary:
   - `fundamental_action_reconstruction/generated/f683_first_exported_selector_output_sign_lift_global_c_v1_rooted_transport_from_s_sel_int_w_break_packet_summary.json`

## Hard limits (no false pass)

`F683` does **not** claim:

- a directed/sign-sensitive physical orientation datum in strict core,
- `Aut(Z_12)`-invariant canonicity,
- kernel-alone/global `QW-2191` discharge,
- strict-core selector closure,
- ToE closure.

