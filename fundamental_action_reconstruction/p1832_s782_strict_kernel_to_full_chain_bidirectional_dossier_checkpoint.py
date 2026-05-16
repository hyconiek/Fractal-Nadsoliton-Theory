#!/usr/bin/env python3
"""P1832 S782 strict kernel->full chain bidirectional dossier checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def status_of(d: dict, fallback: str = "OPEN_UNKNOWN") -> str:
    return d.get("status", fallback)


def main() -> None:
    p1562 = load("p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_checkpoint.json")
    p1662 = load("p1662_s612_strict_full_lagrangian_explicit_density_checkpoint.json")
    p1563 = load("p1563_s513_strict_kernel_to_euler_lagrange_eom_export_checkpoint.json")
    p1706 = load("p1706_s656_strict_full_chain_kernel_to_full_lagrangian_to_eom_dossier_checkpoint.json")
    p1822 = load("p1822_s772_strict_s1_evidence_pack_readiness_checkpoint.json")
    p1831 = load("p1831_s781_strict_tg1_release_readiness_contract_checkpoint.json")

    forward = {
        "kernel_to_coefficients": status_of(p1562),
        "coefficients_to_full_lagrangian": status_of(p1662),
        "full_lagrangian_to_covariant_eom": status_of(p1563),
        "full_chain_dossier": status_of(p1706),
    }

    reverse = {
        "eom_to_variational_origin": "OPEN_OBSTRUCTION_WITH_TRACE",
        "helmholtz_integrability_global": "OPEN_OBSTRUCTION_WITH_TRACE",
        "kernel_identifiability_global": "OPEN_OBSTRUCTION_WITH_TRACE",
        "current_s1_readiness": status_of(p1822),
        "tg1_release_contract": status_of(p1831),
    }

    out = {
        "packet_id": "P1832",
        "stage_id": "S782",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "physics_chain": {
            "forward": forward,
            "reverse": reverse,
        },
        "full_lagrangian_scope_required": [
            "gravity_sector",
            "gauge_sector",
            "fermion_sector",
            "higgs_sector",
            "scalar_mix_sector",
            "curvature_corrections",
            "interaction_terms",
            "covariant_structures",
        ],
        "qg_closure_gates": {
            "renormalization_counterterm_closure": "OPEN_OBSTRUCTION_WITH_TRACE",
            "cutkosky_unitarity": "OPEN_OBSTRUCTION_WITH_TRACE",
            "background_independence": "OPEN_OBSTRUCTION_WITH_TRACE",
            "brst_nilpotency": "OPEN_OBSTRUCTION_WITH_TRACE",
        },
        "technical_progress": "Kernel->coefficients->full Lagrangian->EOM and reverse obligations are unified in one strict bidirectional dossier.",
        "proven": "Forward-chain artifacts exist as explicit checkpoints, but theorem-level reverse/QG closures remain open.",
        "open": "Global reverse variational closure and QG gates remain unresolved; TG1 release contract is still open.",
        "false_pass_risk": "Treating forward export presence as theorem-level QG closure would be a false promotion.",
        "next_honest_step": "Deliver missing nonproxy S1 artifacts and TG1 witness trace, then lift reverse chain (Helmholtz/BRST/Cutkosky/BI) with theorem objects.",
        "lay_explanation": "Mamy mapę całej drogi fizycznej w obie strony: przód jest opisany, ale końcowe testy kwantowej grawitacji nadal czekają na dowody.",
    }

    path = GEN / "p1832_s782_strict_kernel_to_full_chain_bidirectional_dossier_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
