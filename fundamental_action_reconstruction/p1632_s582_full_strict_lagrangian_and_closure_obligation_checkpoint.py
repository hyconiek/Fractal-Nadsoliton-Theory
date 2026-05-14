#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1631 = GEN / "p1631_s581_cover_wise_jacobian_invertibility_summary.json"


def build_full_lagrangian(coeff: dict[str, float]) -> dict[str, str]:
    lam = coeff["lambda_sm_eff"]
    kap = coeff["kappa_gr_eff"]
    eps = coeff["epsilon_mix_eff"]

    return {
        "L_strict_scalar": (
            f"1/2 g^{{μν}}(∇_μ ψ)(∇_ν ψ) - 1/2*{lam:.12f}*ψ^2 - λ4/4*ψ^4"
        ),
        "L_SM_gauge": "-1/4 G^A_{μν}G_A^{μν} -1/4 W^I_{μν}W_I^{μν} -1/4 B_{μν}B^{μν}",
        "L_SM_fermions": "Σ_f i ψ̄_f γ^μ D_μ ψ_f - Σ_f y_f(ψ̄_{fL} H ψ_{fR} + h.c.)",
        "L_SM_higgs": "(D_μ H)^†(D^μ H) - μ_H^2 H^†H - λ_H(H^†H)^2",
        "L_GR": f"{kap:.12f}*R - Λ + c1*R^2 + c2*R_{{μν}}R^{{μν}} + c3*R_{{μναβ}}R^{{μναβ}}",
        "L_mix": f"{eps:.12f}*ψ*R + ξ(H^†H)R + χ ψ^2(H^†H)",
        "L_total": (
            "L_strict_scalar + L_SM_gauge + L_SM_fermions + L_SM_higgs + L_GR + L_mix"
        ),
    }


def main() -> None:
    s31 = json.loads(IN1631.read_text(encoding="utf-8"))
    coeff = s31["forward_backward_chain"]["forward_chain"]["coeff_mapping"]
    lag = build_full_lagrangian(coeff)

    summary = {
        "checkpoint": "P1632_S582",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1632_FULL_STRICT_LAGRANGIAN_EXPORTED",
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "kernel_to_coeff": s31["forward_backward_chain"]["forward_chain"]["kernel_params"],
        "coefficients": coeff,
        "full_lagrangian_density": lag,
        "action": "S = ∫ d^4x sqrt(-g) * L_total",
        "eom_bundle": {
            "psi": "□ψ + m_ψ^2 ψ + λ4 ψ^3 - ε_mix R - 2χψ(H†H) = 0",
            "higgs": "D_μD^μH + μ_H^2H + 2λ_H(H†H)H + ξRH + χψ^2H - Yukawa_sources = 0",
            "gauge": "D_μG_A^{μν}=J_A^ν, D_μW_I^{μν}=J_I^ν, ∂_μB^{μν}=J_Y^ν",
            "metric": "κG_{μν}+Λg_{μν}+H_{μν}[c1,c2,c3]=T_{μν}^{SM}+T_{μν}^{ψ}+T_{μν}^{mix}",
        },
        "strict_core_closure": s31["strict_core_closure"],
        "closure_readout": {
            "state": "OPEN",
            "requires": [
                "export E_selector_internal_source_full_domain",
                "witness W_noncyclic_provider_for_selector_uniqueness_PROVED_GLOBAL",
                "theorem T_qw2191_selector_uniqueness_strict_GLOBAL",
                "theorem T_global_toe_closure_strict",
            ],
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "next_honest_step": "Dowód theorem-level: overlap-compatibility na C_global_noncyclic_cover => global selector uniqueness => strict-core closure.",
        "lay_summary": "Mamy pełny wzór lagranżianu i równania ruchu z kernela strict. Nadal brakuje globalnego dowodu, że lokalne kawałki składają się w jedną, unikalną teorię.",
    }

    out = GEN / "p1632_s582_full_strict_lagrangian_and_closure_obligation_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
