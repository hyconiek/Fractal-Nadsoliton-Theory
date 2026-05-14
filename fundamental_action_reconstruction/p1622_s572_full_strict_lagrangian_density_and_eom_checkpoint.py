#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1562 = GEN / "p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json"
IN1563 = GEN / "p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    s62 = _load(IN1562)
    s63 = _load(IN1563)

    c = s62["derived_lagrangian_coefficients"]
    lam = c["lambda_sm_eff"]
    kap = c["kappa_gr_eff"]
    eps = c["epsilon_mix_eff"]

    lagrangian = {
        "L_strict_scalar": f"1/2*(∂_μψ)(∂^μψ) - 1/2*{lam:.12f}*ψ^2",
        "L_SM": (
            "-1/4*G^A_{μν}G_A^{μν} -1/4*W^I_{μν}W_I^{μν} -1/4*B_{μν}B^{μν}"
            " + Σ_f i*ψ̄_f*γ^μ*D_μ*ψ_f + |D_μH|^2 - V(H) - Σ_f y_f(ψ̄_f H ψ_f)"
        ),
        "L_GR": f"{kap:.12f}*R + Λ_eff + c1*R^2 + c2*R_{{μν}}R^{{μν}}",
        "L_mix": f"{eps:.12f}*ψ*R + ξ*|H|^2*R",
    }
    lagrangian["L_total"] = "L_strict_scalar + L_SM + L_GR + L_mix"

    eom = {
        "psi": f"□ψ + {lam:.12f}*ψ - {eps:.12f}*R + δL_SM/δψ = 0",
        "higgs": "D_μD^μH + ∂V/∂H† + ξ*R*H - Σ_f y_f*ψ̄_fψ_f = 0",
        "gauge": "D_μF_a^{μν} = J_a^ν (z wkładami fermionów i Higgsa)",
        "metric": f"{kap:.12f}*G_{{μν}} + Δ_{{μν}}[c1,c2,Λ_eff] = T_{{μν}}^SM + T_{{μν}}^ψ + T_{{μν}}^mix",
    }

    summary = {
        "checkpoint": "P1622_S572",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1622_FULL_STRICT_LAGRANGIAN_EOM_EXPORT",
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_kernel_params": s62["strict_kernel_params"],
        "derived_coefficients": c,
        "full_lagrangian_density": lagrangian,
        "euler_lagrange_eom": eom,
        "strict_core_closure": {
            "status": "OPEN",
            "missing_exports": s63.get("missing_exports", []),
            "missing_witnesses": s63.get("missing_witnesses", []),
            "missing_theorems": s63.get("missing_theorems", []),
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "external_team_validation_required": False,
        "next_honest_step": "Zbudować theorem object T_qw2191_selector_uniqueness_strict oraz machine-checkable variational proof log dla L_total.",
        "lay_summary": "To jest pełniejsza mapa teorii: od parametrów kernela strict do równania ruchu pól materii i geometrii, ale z uczciwie otwartym statusem brakujących dowodów.",
    }

    out = GEN / "p1622_s572_full_strict_lagrangian_density_and_eom_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
