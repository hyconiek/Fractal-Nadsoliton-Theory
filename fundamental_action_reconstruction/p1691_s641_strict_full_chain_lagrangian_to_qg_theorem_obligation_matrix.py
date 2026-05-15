#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1688 = GEN / "p1688_s638_strict_full_lagrangian_to_sector_eom_export_summary.json"
IN1690 = GEN / "p1690_s640_strict_1loop_counterterm_and_brst_cutkosky_joint_witness.json"
OUT = GEN / "p1691_s641_strict_full_chain_lagrangian_to_qg_theorem_obligation_matrix.json"


def main() -> None:
    p1688 = json.loads(IN1688.read_text(encoding="utf-8"))
    p1690 = json.loads(IN1690.read_text(encoding="utf-8"))

    full_lagrangian_density = {
        "L_gravity": "sqrt(-g)*[(M_Pl^2/2)R - Λ + c_R2 R^2 + c_Ric2 R_{μν}R^{μν} + c_Riem2 R_{μνρσ}R^{μνρσ}]",
        "L_gauge": "sqrt(-g)*[-(1/4) Σ_a Z_a F^a_{μν}F_a^{μν}]",
        "L_fermion": "sqrt(-g)*[Σ_f i ψ̄_f γ^a e_a^{ μ}D_μ ψ_f]",
        "L_higgs": "sqrt(-g)*[(D_μH)^†(D^μH) - μ_H^2(H^†H) - λ_H(H^†H)^2 - ξ_HR(H^†H)R]",
        "L_yukawa": "sqrt(-g)*[-(y_u Q̄H~u + y_d Q̄Hd + y_e L̄He + h.c.)]",
        "L_scalar_phi": "sqrt(-g)*[(1/2)(∇φ)^2 - (1/2)m_φ^2φ^2 - (λ_3/3!)φ^3 - (λ_4/4!)φ^4 - λ_{φH}φ^2(H^†H) - ξ_{φR}φ^2R]",
        "L_mix": "sqrt(-g)*[χ_RG Σ_a R F^a_{μν}F_a^{μν} + CT_{1loop}(R^2,Ricci^2,Riemann^2,SM-mix)]",
    }

    gate_table = p1690.get("gates", {})
    theorem_level_obligations = {
        "counterterm_flow_closure": "KEEP_OPEN",
        "brst_global_nilpotency_and_cohomology": "KEEP_OPEN",
        "cutkosky_full_sector_unitarity": "KEEP_OPEN",
        "background_family_independence": "KEEP_OPEN",
    }

    payload = {
        "checkpoint": "P1691_S641",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "kernel": "K_strict(d)=cos(omega*d+phi)/(1+beta*d^eta)",
        "chain_forward": "K_strict -> coefficients -> full L_total -> sector EOM -> spin2 operator -> local 1-loop/BRST/Cutkosky gates",
        "chain_reverse_obligation": "(QG theorem witnesses) -> consistency of EOM family -> integrable variational origin for full L_total -> strict coefficient closure",
        "full_lagrangian_density": full_lagrangian_density,
        "sector_eom_anchor": p1688.get("sector_eom_export", {}),
        "local_gate_status_from_p1690": gate_table,
        "theorem_level_obligation_matrix": theorem_level_obligations,
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "next_honest_step": "Export theorem-grade counterterm-flow closure and BRST/Cutkosky compatibility on a background family for spin-2 + SM mix.",
        "lay_summary": "Mamy pełną listę składników równania teorii oraz lokalne testy bezpieczeństwa. To wciąż nie jest certyfikat końcowy: brakuje globalnych dowodów kwantowej grawitacji.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
