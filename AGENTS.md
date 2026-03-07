# Project Guardrails

Before continuing work in `fundamental_action_reconstruction`, read:

1. `fundamental_action_reconstruction/K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md`
2. `fundamental_action_reconstruction/K2_STRICT_GATE_KERNEL_DERIVATION_CHAIN_NOTE.md`
3. `fundamental_action_reconstruction/F2_STRICT_GATE_KERNEL_PROVENANCE_AND_FAR_INPUT_CLASSIFICATION_PACKET.md`
4. `fundamental_action_reconstruction/F3_CURRENT_FAR_FRONTIER_KERNEL_ARTIFACT_SENSITIVITY_CLASSIFICATION_PACKET.md`

## Kernel split guardrail

- `K_legacy_ont(d) = alpha_geo*cos(omega*d+phi)/(1+beta_tors*d)` is the legacy ontological/effective kernel.
- `K_strict_gate(d) = cos(omega*d+phi)/(1+beta*d^eta)` is the later operational strict working kernel.
- The repo does **not** currently export a rigorous bridge identifying them.
- Do **not** silently substitute `K_strict_gate` for `K_legacy_ont` or for the canonical parameter layer
  `D_f / alpha_geo / beta_tors`.

## Ontology guardrail

- The nadsoliton itself is the primordial information of the universe in a solitonic state.
- There is no separate informational layer underneath the nadsoliton.
- Preferred internal order remains:
  `nadsoliton -> light -> matter -> emergent observer`.

## FAR continuation rule

- In `A1/A4/A8`, treat `K_strict_gate` only as a later-pipeline operational control/import unless an explicit bridge is added.
- When choosing the next FAR move, prefer routes classified as `kernel-split-robust` by `F3`.
- Do not reopen artifact-sensitive upstream classes unless the task is explicitly to build a `legacy -> strict` bridge or non-bridge theorem.
