#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1713 = GEN / "p1713_s663_strict_metric_residual_zero_attack_plan_checkpoint.json"
IN1662 = GEN / "p1662_s612_strict_full_lagrangian_explicit_density_summary.json"
OUT = GEN / "p1714_s664_strict_metric_first_explicit_emunu_candidate_checkpoint.json"


def main() -> None:
    p1713 = json.loads(IN1713.read_text(encoding="utf-8"))
    p1662 = json.loads(IN1662.read_text(encoding="utf-8"))

    fullL = p1662.get("full_lagrangian_density_explicit", {})

    # First explicit candidate in frozen convention (template-level symbolic split).
    E_munu_candidate = (
        "E_{μν} := (M_Pl^2/2)G_{μν} + Λ g_{μν} + c1*H^{(R2)}_{μν} + c2*H^{(Ric2)}_{μν} + c3*H^{(Riem2)}_{μν} "
        "- T^{gauge}_{μν} - T^{higgs}_{μν} - T^{fermion}_{μν} - T^{scalar}_{μν} - T^{mix}_{μν} = 0"
    )

    T_split = {
        "T_gauge": "-(2/sqrt(-g)) * δL_gauge/δg^{μν}",
        "T_higgs": "-(2/sqrt(-g)) * δL_higgs/δg^{μν}",
        "T_fermion": "-(2/sqrt(-g)) * δL_fermion/δg^{μν}",
        "T_scalar": "-(2/sqrt(-g)) * δL_scalar/δg^{μν}",
        "T_mix": "-(2/sqrt(-g)) * δL_mix/δg^{μν}",
    }

    payload = {
        "checkpoint": "P1714_S664",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": "K_strict -> coefficients -> full explicit L_total -> metric attack plan -> first explicit E_munu candidate",
        "attack_plan_anchor": p1713.get("metric_residual_zero_attack_plan", {}),
        "full_lagrangian_explicit_anchor": fullL,
        "frozen_metric_convention": {
            "signature": "(-,+,+,+)",
            "einstein_tensor": "G_{μν}=R_{μν}-(1/2)g_{μν}R",
            "stress_tensor_definition": "T_{μν}=-(2/sqrt(-g))δS_matter/δg^{μν}",
        },
        "metric_candidate_export": {
            "E_munu_candidate": E_munu_candidate,
            "matter_split": T_split,
            "source_terms_from_full_lagrangian": {
                "L_gravity": fullL.get("L_gravity", "MISSING"),
                "L_gauge": fullL.get("L_gauge", "MISSING"),
                "L_higgs": fullL.get("L_higgs", "MISSING"),
                "L_fermion": fullL.get("L_fermion", "MISSING"),
                "L_scalar": fullL.get("L_scalar", "MISSING"),
                "L_mix": fullL.get("L_mix", "MISSING"),
            },
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": [
            "symbolic_expansion_of_H_R2_H_Ric2_H_Riem2",
            "explicit_metric_sector_EL_residual_zero_certificate",
            "bianchi_ward_cross_consistency_certificate",
            "global_helmholtz_integrability_nonproxy",
            "brst_nilpotency_and_cohomology_nonproxy",
            "cutkosky_unitarity_full_sector",
            "counterterm_flow_renormalization_closure",
            "background_independence_family_theorem",
            "qw2191_selector_source_or_nonclosure_theorem",
        ],
        "next_honest_step": "Rozwinąć E_munu_candidate do jawnej postaci składnikowej H^(R2), H^(Ric2), H^(Riem2) i uruchomić test residual EL_g - E_munu w tej samej konwencji.",
        "lay_summary": "Zaczęliśmy najtrudniejszą część: zapisaliśmy pierwszy jawny kandydat równania grawitacyjnego z rozbiciem wkładów materii. Teraz trzeba go rozpisać do pełnej postaci i sprawdzić testem zero-reszty.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
