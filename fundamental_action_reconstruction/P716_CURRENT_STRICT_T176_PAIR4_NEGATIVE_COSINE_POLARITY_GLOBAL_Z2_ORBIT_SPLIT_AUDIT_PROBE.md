# P716 Current Strict `T176` `pair4` Negative Cosine-Polarity Global `Z2` Orbit-Split Audit Probe (No False-PASS)

Status: `P716_CURRENT_STRICT_T176_PAIR4_NEGATIVE_COSINE_POLARITY_GLOBAL_Z2_ORBIT_SPLIT_AUDIT_PROBE_SPEC_NO_FALSE_PASS`
As of: `2026-03-18`

## Goal

After `P715`, the next honest question is narrower than

```text
why is exact all-root agreement still failing?
```

The sharp current continuation is:

```text
is the remaining exact Z2 split explained completely by the current exported
pair4 chart polarity under the even-anchor fallback, or is there still a
deeper multi-chart mismatch?
```

This probe audits that exact issue.

## Inputs

1. `generated/p715_current_strict_t176_parity_completed_dual_anchor_multiroot_audit_probe.json`
   - rooted dual-anchor results.
2. `generated/a_m_pairm_orientation_projector_operator_strict_core_v1.json` for `m=1..5`
   - exported chart coordinates and representatives.
3. `generated/f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json`
   - exported strict-core payloads `w_break_by_x` and `w_ref_unnormalized_by_x`.

## Audit rule

Test whether the exact residual split found by `P715`

```text
same orbit roots     = {pair1, pair2, pair3, pair5}
negated orbit roots  = {pair4}
```

is fully explained by the following joint profile on current exports:

1. `pair4` is the unique currently exported **negative cosine-axis** chart,
2. `pair2` and `pair3` are the positive cosine-axis charts,
3. the dual-anchor rule uses the even fallback `w_ref_unnormalized` on
   `pair2`, `pair3`, `pair4`,
4. `dot(w_ref_unnormalized, u_4) < 0` while the corresponding cosine-anchor
   scalars for `pair2`, `pair3` are positive,
5. the `pair4` rooted section is exactly the **global negation** of the
   reference section recovered from the other roots.

## Result discipline

The probe passes only if it shows that no additional root-specific transport
failure is needed to explain the current exact split:

- all roots remain supported (`P715`),
- the projective orbit remains common,
- and the remaining exact branch split is exhausted by the unique
  `pair4` negative cosine polarity under the current even-anchor fallback.

This is still a boundary result, not a discharge.

## Hard limits

This probe does **not**:

- claim a strict-core directed/sign-sensitive physical orientation datum,
- claim kernel-alone/global `QW-2191` discharge,
- claim `T176` discharge,
- claim that the residual split is physically irrelevant in strict core,
- claim ToE closure.
