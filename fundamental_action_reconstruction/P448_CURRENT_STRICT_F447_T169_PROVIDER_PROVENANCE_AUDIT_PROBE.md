# P448 Current Strict `F447` `T169` Provider Provenance Audit Probe (No False‑PASS)

Status: `P448_EXECUTED_CURRENT_STRICT_F447_T169_PROVIDER_PROVENANCE_AUDIT_PROBE_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

`F447` exports a per‑site value provider intended to discharge the `T169 → T168 → T167 → T166` diagonal/local lane
computability bottleneck by instantiating:

```text
vpsi[0..11], g4[0..11], g6[0..11]  +  Sigma_psi{k}_psi{k+6} (k=0..5).
```

Downstream harnesses (`P437 → P434`) can then compute `F2(d)` and `theta_*`.

However, `P437` explicitly tags its run status as:

```text
PASS_COMPUTED_FROM_INPUTS_REQUIRES_PROVENANCE_REVIEW
```

meaning the computation is reproducible but the **classification / strict‑derived provenance** of the upstream value
provider must be reviewed before any discharge is verbally promoted.

`P448` is that mechanical provenance audit:

1. it validates schema completeness and strict-input references of the `F447` provider,
2. it lists all nontrivial selection premises used by the lift (as exported by the provider),
3. it compares the provider against the `T169` acceptance criteria and the `P446` “real discharge blueprint” checklists,
4. it outputs a conservative recommendation: keep `strict_derived` **only if** a theorem-level admissibility/uniqueness
   story is exported; otherwise downgrade to `strict_extension_only` (no false‑PASS).

This probe does **not** add a new selector premise and does **not** claim `QW-2191` discharge or ToE closure.

## Execution

Default (audits the designated generated provider):

```bash
python3 fundamental_action_reconstruction/p448_current_strict_f447_t169_provider_provenance_audit_probe.py
```

Audit an explicit provider JSON:

```bash
python3 fundamental_action_reconstruction/p448_current_strict_f447_t169_provider_provenance_audit_probe.py --provider <path>
```

## Output

`P448` persists:

- `fundamental_action_reconstruction/generated/p448_current_strict_f447_t169_provider_provenance_audit_probe.json`
- `fundamental_action_reconstruction/generated/p448_current_strict_f447_t169_provider_provenance_audit_probe_summary.json`

## Hard limits (no false pass)

`P448` does not claim:

1. discharge of `T169` / `T168` / `T167` / `T166`,
2. that the audited provider is strict-derived (it only audits and recommends),
3. global `QW-2191` discharge,
4. ToE closure.

