# P470 Current Strict `pair1..pair5` Oriented Transport (α mod 2π) — Full Cocycle Audit Probe (No False‑PASS)

Status: `P470_DRAFT_EXPECTED_EXECUTED_BY_p470_SCRIPT_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

Audit the oriented-transport (`α mod 2π`) convention-layer lift exported by `F467` on `{pair1..pair5}`:

1. orthogonality + involution of each exported transport operator,
2. vector-level transport: `O_ij u_i = u_j` on the exported representative vectors,
3. full triple cocycle/path-independence at vector level for all triples on `{pair1..pair5}`.

This is a probe-level numeric audit on the exported current `n=12` artifacts only.

## Inputs

- `F456`, `F462`, `F464`, `F465`: exported representative vectors `u_1..u_5`.
- `F461`: exported `O_12` and `alpha_12 mod 2π`.
- `F467`: exported oriented transport operators (`O_13`, `O_14`, `O_15`, `O_23`, `O_24`, `O_25`, `O_34`, `O_35`, `O_45`) and the oriented vector section.

## Outputs

Executed by:

```text
python3 fundamental_action_reconstruction/p470_current_strict_pair12345_oriented_transport_full_cocycle_audit_probe.py
```

Exports:

- `fundamental_action_reconstruction/generated/p470_current_strict_pair12345_oriented_transport_full_cocycle_audit_probe.json`
- `fundamental_action_reconstruction/generated/p470_current_strict_pair12345_oriented_transport_full_cocycle_audit_probe_summary.json`

## Hard limits

This probe does **not** claim:

- a theorem-level global selector atlas on `C_v1`,
- discharge of `QW-2191`,
- sign-sensitive physical orientation datum (this is a tracked gauge/convention layer),
- ToE closure.

