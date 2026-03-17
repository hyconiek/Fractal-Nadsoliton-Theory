# P714 Current Strict `T176` `w_break` Parity Root-Support Profile Audit Probe (No False-PASS)

Status: `P714_CURRENT_STRICT_T176_W_BREAK_PARITY_ROOT_SUPPORT_PROFILE_AUDIT_PROBE_SPEC_NO_FALSE_PASS`
As of: `2026-03-17`

## Goal

After `P713`, the next honest question is no longer merely:

```text
which roots survive?
```

but:

```text
why do only those roots survive on the current exported strict-core payload?
```

This probe audits whether the current `w_break` root-support profile is explained by a clean parity mechanism on `Z_12`.

## Inputs

1. `generated/f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json`
   - exported strict-core reflection-breaking payload `w_break_by_x`.
2. `generated/a_m_pairm_orientation_projector_operator_strict_core_v1.json` for `m=1..5`
   - exported representative vectors `u_m` and their chart coordinates.
3. `generated/p713_current_strict_t176_multiroot_rooted_sign_lift_root_independence_audit_probe_summary.json`
   - current supported-root corridor result.

## Audit question

Test whether the current surviving supported-root corridor

```text
{pair1, pair5}
```

is explained by the following strict parity profile:

1. `w_break` is odd under reflection `x -> -x` on `Z_12`,
2. the currently exported representatives on `pair1` and `pair5` are odd-axis states,
3. the currently exported representatives on `pair2`, `pair3`, `pair4` are even-axis states,
4. therefore the linear anchor `Σ_x w_break(x) u_m(x)` vanishes on the even-axis charts and survives only on the odd-axis charts.

## Acceptance

The probe passes only if it can show all of:

1. `w_break` is odd (within tolerance),
2. `u_1`, `u_5` are odd and have nonzero `w_break` anchor,
3. `u_2`, `u_3`, `u_4` are even and have zero `w_break` anchor,
4. the supported roots reported by `P713` match exactly that parity support profile.

## Hard limits

This probe does **not**:

- claim that no future strict-core provider can anchor all charts,
- claim `T176` discharge,
- claim kernel-alone/global `QW-2191` discharge,
- promote any directed/sign-sensitive physical orientation datum into strict core,
- claim ToE closure.
