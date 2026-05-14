#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1661 = GEN / "p1661_s611_strict_h2_boundary_flux_lemma_scaffold_summary.json"


def main() -> None:
    s61 = json.loads(IN1661.read_text(encoding="utf-8"))

    explicit_density = {
        "L_scalar": "sqrt(-g)*[1/2 g^{μν}(∇_μφ)(∇_νφ) - (m_φ^2/2)φ^2 - (λ_3/3!)φ^3 - (λ_4/4!)φ^4]",
        "L_gauge": "sqrt(-g)*[-1/4 G^A_{μν}G_A^{μν} - 1/4 W^I_{μν}W_I^{μν} - 1/4 B_{μν}B^{μν}]",
        "L_fermion": "sqrt(-g)*[Σ_f i ψ̄_f γ^a e_a^{ μ}D_μ ψ_f - Σ_f y_f(ψ̄_{fL}Hψ_{fR}+h.c.)]",
        "L_higgs": "sqrt(-g)*[(D_μH)^†(D^μH) - μ_H^2(H^†H) - λ_H(H^†H)^2]",
        "L_gravity": "sqrt(-g)*[(M_Pl^2/2)R - Λ + c1 R^2 + c2 R_{μν}R^{μν} + c3 R_{μναβ}R^{μναβ}]",
        "L_mix": "sqrt(-g)*[ξ_HR(H^†H)R + ξ_φR φ^2R + λ_{φH}φ^2(H^†H)]",
    }

    summary = {
        "checkpoint": "P1662_S612",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "KEEP_OPEN_P1662_EXPLICIT_FULL_LAGRANGIAN_DENSITY_EXPORTED",
        "strict_only": True,
        "legacy_bridge_used": False,
        "kernel_definition": s61["kernel_definition"],
        "forward_chain_context": "K_strict -> coefficients -> full L_total -> EOM",
        "full_lagrangian_density_explicit": explicit_density,
        "L_total_explicit": "L_scalar + L_gauge + L_fermion + L_higgs + L_gravity + L_mix",
        "reverse_chain_status": "OPEN",
        "qg_gates": s61["qg_gates"],
        "next_honest_step": "Przepiąć P1662 jako referencyjny L_total do eksportu jawnych równań EOM sektor-po-sektor i aktualizacji H2/H3 witnessów.",
        "lay_summary": "Zapisaliśmy pełny lagranżian w jednej jawnej formule, żeby kolejne kroki nie opierały się na skrótach. Dzięki temu łatwiej sprawdzać, czy równania ruchu naprawdę odpowiadają całej teorii.",
    }

    out = GEN / "p1662_s612_strict_full_lagrangian_explicit_density_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
