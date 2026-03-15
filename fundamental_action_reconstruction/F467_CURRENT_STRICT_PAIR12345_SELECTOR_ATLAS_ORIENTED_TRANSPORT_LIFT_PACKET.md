# F467 Current Strict `pair1..pair5` Selector‑Atlas — Oriented Transport (α mod 2π) Lift Packet (No False‑PASS)

Status: `F467_EXECUTED_CURRENT_STRICT_PAIR12345_SELECTOR_ATLAS_ORIENTED_TRANSPORT_LIFT_PACKET_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

After `F466/P469/N510`, the repo exports a lane‑scoped five‑chart selector‑atlas ingredient on `{pair1..pair5}` at **projector level**,
with full triple cocycle/path‑independence audits. This is sign‑gauge‑safe, but it is **axis-only** (angles only `mod π`) and therefore
does not provide oriented transport (`mod 2π`) at the level of representative vectors.

This packet performs the next honest minimal continuation under strict no‑false‑PASS discipline:

```text
export a tracked gauge/convention lift of the five‑chart projector atlas to oriented transport (α mod 2π),
induced by the currently exported representative vectors u_1..u_5, and audit full triple cocycle at vector level.
```

This remains:

- **lane‑scoped** (exported‑artifact overlap on the declared `n=12` Fourier carrier; not a global open cover on `C_v1`),
- explicitly a **sign‑tracked convention layer** (not a sign‑sensitive physical orientation datum),
- strictly below strict-core selector closure and below any global discharge of `QW-2191`.

## Strict‑admissible inputs reused

1. `F456`, `F462`, `F464`, `F465`
   - exported representative vectors and projector operators `A_1..A_5` on `{pair1..pair5}`.
2. `F461`
   - exported `O_12` and `alpha_12 mod 2π` derived from strict sigma‑int corridor transition data.
3. `F466`
   - projector-level five‑chart ingredient with full triple cocycle data (used only as prior infrastructure; this packet does not promote it to a global atlas on `C_v1`).

## Exported artifacts

Executed by:

```text
python3 fundamental_action_reconstruction/f467_current_strict_pair12345_selector_atlas_oriented_transport_lift_packet.py
```

Exports:

1. oriented theta family (tracked convention):
   - `fundamental_action_reconstruction/generated/theta_family_pair12345_oriented_mod_2pi_strict_convention_v1.json`
2. oriented chart transport operators (α mod 2π; tracked convention):
   - `fundamental_action_reconstruction/generated/o13_pair1_pair3_selector_chart_transport_operator_oriented_alpha13_mod_2pi_strict_convention_v1.json`
   - `fundamental_action_reconstruction/generated/o14_pair1_pair4_selector_chart_transport_operator_oriented_alpha14_mod_2pi_strict_convention_v1.json`
   - `fundamental_action_reconstruction/generated/o15_pair1_pair5_selector_chart_transport_operator_oriented_alpha15_mod_2pi_strict_convention_v1.json`
   - `fundamental_action_reconstruction/generated/o23_pair2_pair3_selector_chart_transport_operator_oriented_alpha23_mod_2pi_strict_convention_v1.json`
   - `fundamental_action_reconstruction/generated/o24_pair2_pair4_selector_chart_transport_operator_oriented_alpha24_mod_2pi_strict_convention_v1.json`
   - `fundamental_action_reconstruction/generated/o25_pair2_pair5_selector_chart_transport_operator_oriented_alpha25_mod_2pi_strict_convention_v1.json`
   - `fundamental_action_reconstruction/generated/o34_pair3_pair4_selector_chart_transport_operator_oriented_alpha34_mod_2pi_strict_convention_v1.json`
   - `fundamental_action_reconstruction/generated/o35_pair3_pair5_selector_chart_transport_operator_oriented_alpha35_mod_2pi_strict_convention_v1.json`
   - `fundamental_action_reconstruction/generated/o45_pair4_pair5_selector_chart_transport_operator_oriented_alpha45_mod_2pi_strict_convention_v1.json`
3. oriented vector section on `{pair1..pair5}` with full triple cocycle audits:
   - `fundamental_action_reconstruction/generated/u_12345_pair12345_chart_glued_orientation_vector_section_oriented_mod_2pi_strict_convention_v1.json`
4. oriented atlas object (lane‑scoped; convention layer):
   - `fundamental_action_reconstruction/generated/selector_atlas_pair12345_oriented_transport_mod_2pi_strict_convention_v1.json`
5. summary:
   - `fundamental_action_reconstruction/generated/f467_current_strict_pair12345_selector_atlas_oriented_transport_lift_packet_summary.json`

## Meaning (no false‑PASS)

This packet means only:

1. the repo exports a **tracked** oriented‑transport lift (angles `mod 2π`) of the already exported five‑chart projector atlas ingredient,
2. the lift is induced by the currently exported representative vectors `u_1..u_5` and is therefore a gauge/convention layer,
3. full triple cocycle/path‑independence is audited at the level of those exported vectors, without claiming a physical sign convention.

## Hard limits

This packet does **not** claim:

1. a sign‑sensitive physical orientation datum (lifting residual `Z2` as physics),
2. strict-core selector closure / admissible `S_sel_int`,
3. a global selector atlas on the full strict domain `C_v1`,
4. global discharge of `QW-2191`,
5. ToE closure.

