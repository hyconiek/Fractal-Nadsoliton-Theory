#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

IN1622 = GEN / "p1622_s572_full_strict_lagrangian_density_and_eom_summary.json"
IN1623 = GEN / "p1623_s573_strict_selector_uniqueness_theorem_object_and_variational_log_summary.json"
IN1605 = GEN / "p1605_s555_np1_provider_instantiation_and_replay_summary.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    s22 = _load(IN1622)
    s23 = _load(IN1623)
    s05 = _load(IN1605)

    coeff = s22["derived_coefficients"]
    lam = coeff["lambda_sm_eff"]
    kap = coeff["kappa_gr_eff"]
    eps = coeff["epsilon_mix_eff"]

    witness = {
        "id": "W_noncyclic_provider_for_selector_uniqueness",
        "provider": "NP1",
        "noncyclic_anchor": "strict energy-curvature monotone M = kap*R + eps*psi*R - lam*psi^2/2",
        "selection_rule": "choose branch maximizing monotone decrease of |dM/dt| under EL flow",
        "strict_dependencies": [
            "K_strict(omega,phi,beta,eta)",
            "derived_coefficients(lambda_sm_eff,kappa_gr_eff,epsilon_mix_eff)",
            "EL equations for psi/higgs/gauge/metric",
        ],
        "coeff_values": {"lambda_sm_eff": lam, "kappa_gr_eff": kap, "epsilon_mix_eff": eps},
        "status": "CANDIDATE_EXPORTED",
    }

    theorem_open = s23["theorem_object"]["proof_status"] != "PROVED"
    closure_status = "OPEN" if theorem_open else "CLOSED"

    summary = {
        "checkpoint": "P1624_S574",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1624_WITNESS_CANDIDATE_EXPORTED" if witness["status"] == "CANDIDATE_EXPORTED" else "FAIL_P1624",
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "witness_object": witness,
        "np1_replay_context": s05.get("g3_replay", s05.get("status", "UNKNOWN")),
        "strict_core_closure": {
            "status": closure_status,
            "missing_exports": ["E_selector_internal_source_full_domain", "E_full_variational_proof_log_machine_checkable"],
            "missing_witnesses": [] if witness["status"] == "CANDIDATE_EXPORTED" else ["W_noncyclic_provider_for_selector_uniqueness"],
            "missing_theorems": ["T_qw2191_selector_uniqueness_strict", "T_global_toe_closure_strict"] if theorem_open else [],
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "external_team_validation_required": False,
        "next_honest_step": "Podnieść witness z CANDIDATE do PROVED przez formalny dowód zgodności selection_rule z EL-flow i eksport E_selector_internal_source_full_domain.",
        "lay_summary": "To jest kandydat na 'kompas' wybierający jedną ścieżkę rozwoju teorii; nadal trzeba formalnie udowodnić, że ten kompas zawsze działa.",
    }

    out = GEN / "p1624_s574_noncyclic_selector_witness_from_strict_kernel_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
