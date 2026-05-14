#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1643 = GEN / "p1643_s593_strict_full_lagrangian_and_bidirectional_eom_obligation_summary.json"


def main() -> None:
    s43 = json.loads(IN1643.read_text(encoding="utf-8"))

    helmholtz_conditions = [
        {
            "id": "H1_self_adjoint_frechet",
            "statement": "Linearized Euler operator for (φ,H,ψ_f,A^A_μ,W^I_μ,B_μ,g_{μν}) is formally self-adjoint on admissible strict domain.",
            "status": "OPEN",
            "evidence_now": "local/sector checks only",
            "missing_export": "global proof object with boundary/regularity assumptions",
        },
        {
            "id": "H2_cross_variation_symmetry",
            "statement": "Mixed second variational derivatives commute: δE_i/δq_j = δE_j/δq_i for full strict operator basis.",
            "status": "OPEN",
            "evidence_now": "symbolic proxies for reduced subsystems",
            "missing_export": "machine-checkable mixed-block symmetry witness",
        },
        {
            "id": "H3_gauge_bianchi_compatibility",
            "statement": "Gauge/Bianchi identities are compatible with strict source terms and mixed couplings.",
            "status": "PARTIAL",
            "evidence_now": "sector-level identities recorded",
            "missing_export": "global compatibility theorem across overlap patches",
        },
        {
            "id": "H4_boundary_term_control",
            "statement": "All integrations by parts produce controlled boundary terms under strict admissibility class.",
            "status": "OPEN",
            "evidence_now": "implicit assumptions, not exported as theorem",
            "missing_export": "explicit admissibility theorem and decay/regularity witness",
        },
    ]

    reverse_reconstruction_map = {
        "target": "EOM -> L_total",
        "strategy": [
            "Prove H1..H4 on full strict operator domain",
            "Construct variational potential L_total modulo total divergence",
            "Fix divergence class using canonical normalization constraints",
            "Export uniqueness class and handoff to B5/B6",
        ],
        "non_legacy_bridge": True,
    }

    summary = {
        "checkpoint": "P1644_S594",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "KEEP_OPEN_P1644_HELMHOLTZ_WITNESS_SCAFFOLD_EXPORTED",
        "route_target": s43["route_target"],
        "strict_only": True,
        "legacy_bridge_used": False,
        "full_lagrangian_density_reference": s43["full_lagrangian_density"],
        "eom_system_reference": s43["eom_system"],
        "helmholtz_theorem_requirements": helmholtz_conditions,
        "reverse_reconstruction_map": reverse_reconstruction_map,
        "strict_core_closure": {
            "status": "OPEN",
            "reason": "B4 theorem-level witness remains incomplete until H1..H4 are fully discharged",
        },
        "next_honest_step": "Zbudować machine-checkable proof object dla H2 (cross-variation symmetry) na pełnej bazie operatorów strict; to odblokuje konstrukcję L_total z EOM.",
        "lay_summary": "Aby odzyskać lagranżian z równań ruchu, trzeba formalnie sprawdzić cztery warunki zgodności. Mamy plan i listę warunków, ale dowód całościowy jeszcze trwa.",
    }

    out = GEN / "p1644_s594_strict_eom_to_lagrangian_helmholtz_witness_scaffold_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
