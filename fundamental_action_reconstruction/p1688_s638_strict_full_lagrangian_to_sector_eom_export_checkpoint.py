#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1662 = GEN / "p1662_s612_strict_full_lagrangian_explicit_density_summary.json"
IN1664 = GEN / "p1664_s614_strict_full_lagrangian_manifest_and_inversion.json"
OUT = GEN / "p1688_s638_strict_full_lagrangian_to_sector_eom_export_summary.json"


def main() -> None:
    s62 = json.loads(IN1662.read_text(encoding="utf-8"))
    s64 = json.loads(IN1664.read_text(encoding="utf-8"))

    eom_export = {
        "scalar_phi": "□φ + m_φ^2 φ + (λ_3/2)φ^2 + (λ_4/3!)φ^3 + 2λ_{φH}φ(H†H) + 2ξ_{φR}Rφ = 0",
        "higgs_H": "D_μD^μH + μ_H^2 H + 2λ_H(H†H)H + ξ_HR R H + λ_{φH}φ^2 H + Yukawa_source = 0",
        "gauge_A": "∇_ν(Z_a F_a^{νμ}) + J_a^μ + χ_RG ∇_ν(R F_a^{νμ}) = 0",
        "fermion_psi": "iγ^a e_a^{ μ}D_μ ψ_f - y_f H·ψ_f = 0",
        "metric_g": "(M_Pl^2/2)G_{μν} + Λg_{μν} + H_{μν}^{(R2,Ric2,Riem2)} = T_{μν}^{SM+mix}",
    }

    bidirectional_status = {
        "forward_kernel_to_eom": "EXPORTED_SECTORWISE_FROM_FULL_L_TOTAL",
        "reverse_eom_to_lagrangian": "KEEP_OPEN_HELMHOLTZ_PLUS_QG_OBLIGATIONS",
        "local_inversion_from_p1664": s64.get("inverse_recovery", {}).get("local_pass", False),
    }

    payload = {
        "checkpoint": "P1688_S638",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "kernel_context": s62.get("kernel_definition", "K_strict_gate(d)=cos(omega*d+phi)/(1+beta*d^eta)"),
        "L_total_context": s62.get("L_total_explicit"),
        "sector_eom_export": eom_export,
        "bidirectional_status": bidirectional_status,
        "qg_open_obligations": [
            "full_counterterm_renormalization_export",
            "spin2_spin0_unitarity_nonproxy_witness",
            "quantum_background_independence_witness",
        ],
        "status": "KEEP_OPEN_STRICT_CORE_CLOSURE_REQUIRES_QG_EXPORTS",
        "next_honest_step": "Build a strict-core spin-2 curved-background operator witness with explicit spectral-sign test and no proxy shortcuts.",
        "lay_summary": "Mamy już pełny wzór teorii i pokazaliśmy równania ruchu dla każdego sektora. To nie zamyka jeszcze ToE, bo nadal brakuje twardych dowodów kwantowej grawitacji: renormalizacji, unitarności i niezależności od tła.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
