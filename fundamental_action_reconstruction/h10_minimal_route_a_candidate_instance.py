#!/usr/bin/env python3
import json
from pathlib import Path

candidate = {
    "id": "H10",
    "carrier": "A_1_cand",
    "status": "PASS_PARTIAL_MINIMAL_ROUTE_A_CANDIDATE_INSTANCE_CREATED",
    "as_of": "2026-03-06",
    "pair": {
        "label": "pair1",
        "basis": ["c1", "s1"],
        "plane": "V_1 = span{c1, s1}"
    },
    "matrix": {
        "form": [["a_1", "b_1"], ["b_1", "d_1"]],
        "symmetry": "real_symmetric",
        "entries": {
            "a_1": "UNRESOLVED_FROM_OPERATOR_EXPORT_CHAIN",
            "b_1": "UNRESOLVED_FROM_OPERATOR_EXPORT_CHAIN",
            "d_1": "UNRESOLVED_FROM_OPERATOR_EXPORT_CHAIN"
        }
    },
    "provenance": {
        "lane": "hypothesis_extension_only",
        "derived_from_operator_chain": False,
        "strict_core_component": False,
        "required_future_source": "explicit export from E_1, G_light, R_mat, O_obs or exported composite A_1"
    },
    "hard_limits": {
        "theorem_level_pass": False,
        "full_closure_pass": False,
        "coefficients_computed": False,
        "route_A_discharged": False,
        "q_w_2191_discharged": False
    }
}

summary = {
    "id": "H10",
    "title": "minimal route A candidate instance",
    "status": "PASS_PARTIAL_MINIMAL_ROUTE_A_CANDIDATE_INSTANCE_CREATED",
    "as_of": "2026-03-06",
    "artifact": "generated/route_a_minimal_candidate_instance.json",
    "frontier": {
        "H10_B1": "a minimal Route A candidate instance exists, but no provenance-valid exported A_1 derived from the operator chain exists yet",
        "H9_B1": "no actual Route A instance and no actual Route B instance exists for pair1 in the current repository exports is reduced to provenance-validity level for Route A",
        "T12_B1": "strict-core typing judgment with totality and uniqueness remains undischarged",
        "T2_B1": "bridge theorem still specified but not discharged",
        "C32_B2": "raw cross-pair overlap route remains degenerate"
    },
    "hard_limits": candidate["hard_limits"]
}

cand_out = Path("fundamental_action_reconstruction/generated/route_a_minimal_candidate_instance.json")
summary_out = Path("fundamental_action_reconstruction/generated/h10_minimal_route_a_candidate_instance_summary.json")
cand_out.write_text(json.dumps(candidate, indent=2) + "\n", encoding="utf-8")
summary_out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(cand_out)
print(summary_out)
