# P442 Current Strict `T168` Vacuum + Self‑Coupling Provider Frontier Audit Probe

Status: `P442_EXECUTED_CURRENT_STRICT_T168_VACUUM_AND_SELF_COUPLING_PROVIDER_FRONTIER_AUDIT_PROBE_NO_FALSE_PASS`  
As of: `2026-03-14`

## Goal

`P438/P441` currently recommend the next strict move:

```text
T168 (strict-derived provider for vpsi[0..11], g4[0..11], g6[0..11])
```

because `P437` (the `N477` evaluation harness) is not computable on current repo exports.

`P442` audits this frontier at a practical level:

1. confirm which **strict-chain scalar** vacuum/self-coupling inputs already exist (`QW-2122/2123/2124`),
2. confirm that the repo still does **not** export a strict-derived **lift** of those scalars into the canonical
   per-site arrays `vpsi/g4/g6` required by `P437`,
3. keep the “no false pass” boundary explicit (no promotion of mapping premises into strict core).

## Strict scalar sources (already exist)

The strict scalar-side vacuum closure chain exists and is deterministic:

1. `QW-2122` exports the psi-potential scalar parameters (broken branch) including:
   - `lambda_psi_strict`,
   - `rho_star_sq`,
   - `diag_floor`.
2. `QW-2123` exports a strict branch selection rule (broken branch required) tied to `lambda_min(K_total)<0`.
3. `QW-2124` promotes the scalar vacuum closure to a branch-resolved strict PASS.

These are necessary ingredients for vacuum-floor consistency (`R15`), but they do **not** by themselves discharge
`T168`.

## What remains missing (the actual `T168` provider)

`T168` requires a strict-derived value-provider that instantiates:

```text
vpsi[0..11], g4[0..11], g6[0..11]
```

as the **canonical** FIN constant-vacuum and local self-coupling families used by `QW-2165/2166` and consumed by
`N477/P437`.

On current exports:

- the designated harness input `p437_input_vpsi_g4_g6_candidate.json` intentionally keeps these values `null`,
  and the repo contains only toy/probe/extension instantiations elsewhere.
- therefore there is still an explicit missing “lift/mapping” layer from the strict scalar chain into the canonical
  per-site families, with no hidden “choose a representative / choose a branch / choose a site” slot.

`N478` (audited toy-level by `P443`) additionally confirms that fixing the strict scalar outputs alone does not
canonically determine the needed per-site arrays; the missing lift/mapping is structural, not cosmetic.

## Outputs

`P442` writes:

- `fundamental_action_reconstruction/generated/p442_current_strict_t168_vacuum_and_self_coupling_provider_frontier_audit_probe.json`
- `fundamental_action_reconstruction/generated/p442_current_strict_t168_vacuum_and_self_coupling_provider_frontier_audit_probe_summary.json`

## Hard limits (no false pass)

`P442` does not claim:

1. discharge of `T168`,
2. any strict-derived canonical per-site vacuum vector `vpsi[0..11]`,
3. any strict-derived per-site self-coupling families `g4/g6`,
4. discharge of `T167/T166` or `QW-2191`,
5. ToE closure.
