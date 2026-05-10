# P786 Current QW2093 Kernel-Invariants-Only Formula Layer Audit Probe

Status: `P786_CURRENT_QW2093_KERNEL_INVARIANTS_ONLY_FORMULA_LAYER_BLOCKED_NONEXPORT`
As of: `2026-03-19`

## Goal

Run one narrow audit on the current `QW-2093` formula layer:

1. test whether `sin2_theta_w_eff` is already kernel-invariants-only,
2. test whether `alpha_s_boundary_mu0_alpha0` is already kernel-invariants-only,
3. block any silent promotion if frozen ansatz, legacy touchpoint, or GeV-mass-chain residue remains.

## Scope

`P786` does not try to improve the formulas.
It only decides whether the present layer may be honestly promoted into the
minimal strict bridge without hidden semantic or calibration leakage.

## Main Checks

1. detect `alpha_geo` residue in the `sin2` formula,
2. detect frozen-source labeling in the EW chain,
3. detect `m_top/m_bottom` dependence in the `alpha_s` boundary,
4. detect explicit hierarchy-log ansatz labeling in the alpha_s gate,
5. check whether the mass-chain foundation is already promoted to strict first principles,
6. keep the `QW-2064` wide-CI warning visible.

## Hard Limit

No object is promoted by `P786`.
If any residue remains, the result stays `blocked/nonexport` even if downstream
phenomenology looks numerically good.
