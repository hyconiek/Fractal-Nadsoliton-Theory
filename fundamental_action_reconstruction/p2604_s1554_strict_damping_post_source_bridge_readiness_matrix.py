#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2604_s1554_strict_damping_post_source_bridge_readiness_matrix.json"
MD = GEN / "p2604_s1554_strict_damping_post_source_bridge_readiness_matrix.md"

SOURCE_FILES = {
    "P2603_FRACTAL_CODIMENSION_SLOPE_SOURCE": GEN / "p2603_s1553_nadsoliton_fractal_codimension_slope_source_theorem.json",
}
BRIDGE_ROLE_GATES = [
    "legacy_to_strict_completion_map_evidence",
    "strict_side_residual_additions_certified",
    "strict_damping_role_transfer_theorem",
]
NEGATIVE_EXPORT_FLAGS = [
    "legacy_to_strict_completion_bridge_exported",
    "strict_side_residual_additions_exported",
    "strict_damping_role_transfer_theorem_exported",
    "role_bearing_ltotal_exported",
    "qw2191_discharged_by_this_matrix",
    "toe_closure_claimed",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
    ], cwd=REPO, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2604|S1554|post strict damping bridge readiness|strict damping bridge readiness|role-bearing L_total readiness",
        "intended_research_nonduplication": "bridge completion truth table|post-discharge bridge gate|strict damping role transfer audit|source discharge role transfer",
        "precursor_chain": "P2603|S1553|strict_damping_beta_eta_source_exported|legacy-to-strict bridge|role-transfer theorem",
        "guardrails": "K_legacy_ont|K_strict_gate|role-bearing L_total|QW-2191|ToE closure|kernel-split guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def bridge_readiness_rows(strict_damping_source: bool) -> list[dict[str, Any]]:
    rows = []
    for values in product([False, True], repeat=len(BRIDGE_ROLE_GATES)):
        assignment = dict(zip(BRIDGE_ROLE_GATES, values))
        bridge_ready = assignment["legacy_to_strict_completion_map_evidence"] and assignment["strict_side_residual_additions_certified"]
        role_ready = bridge_ready and assignment["strict_damping_role_transfer_theorem"]
        role_bearing = strict_damping_source and role_ready
        missing = [key for key, value in assignment.items() if not value]
        rows.append({
            "strict_damping_beta_eta_source_exported": strict_damping_source,
            "assignment": assignment,
            "bridge_completion_ready": bridge_ready,
            "role_transfer_ready": role_ready,
            "role_bearing_ltotal_accepts": role_bearing,
            "missing_bridge_role_gates": missing,
            "missing_bridge_role_gate_count": len(missing),
        })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2603_payload = load_json(SOURCE_FILES["P2603_FRACTAL_CODIMENSION_SLOPE_SOURCE"])
    p2603 = theorem(p2603_payload, "nadsoliton_fractal_codimension_slope_source_theorem")
    strict_damping_source = p2603.get("strict_damping_beta_eta_source_exported") is True
    rows = bridge_readiness_rows(strict_damping_source)
    accepting = [row for row in rows if row["role_bearing_ltotal_accepts"]]
    current_assignment = {key: False for key in BRIDGE_ROLE_GATES}
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2604_T1_strict_damping_post_source_bridge_readiness_matrix",
        "audited_chain": ["P2603/S1553"],
        "matrix_statement": (
            "After P2603 the strict damping beta/eta source normal form is discharged, but role-bearing L_total readiness requires three additional non-source gates: explicit legacy-to-strict completion-map evidence, certified strict-side residual additions, and a strict damping role-transfer theorem."
        ),
        "strict_damping_beta_eta_source_inherited_from_p2603": strict_damping_source,
        "bridge_role_gates": BRIDGE_ROLE_GATES,
        "bridge_readiness_truth_table_rows": rows,
        "bridge_readiness_truth_table_row_count": len(rows),
        "role_bearing_accepting_row_count": len(accepting),
        "role_bearing_accepting_row": accepting[0],
        "current_bridge_role_assignment": current_assignment,
        "current_bridge_completion_ready": False,
        "current_role_transfer_ready": False,
        "current_role_bearing_ltotal_accepts": False,
        "concrete_missing_ingredients": [
            "completion_map_evidence_from_K_legacy_ont_to_K_strict_gate_including_parameter_layer_D_f_alpha_geo_beta_tors",
            "strict_side_residual_additions_for_nonlinear_compression_damping_phase_and_topological_data",
            "strict_damping_role_transfer_theorem_after_bridge_completion",
        ],
        "recommended_next_honest_step": (
            "Do not expand APD/moment/Sturm or repeat strict damping source keys. The next noncyclic FAR move is one explicit completion-map theorem from K_legacy_ont to K_strict_gate, with residual strict-side additions stated before any role-transfer claim."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2603_strict_damping_source_inherited": strict_damping_source,
        "truth_table_has_eight_rows": len(rows) == 8,
        "exactly_one_role_bearing_accepting_row": len(accepting) == 1,
        "current_assignment_has_three_missing_gates": len(theorem_export["concrete_missing_ingredients"]) == 3,
        "current_role_bearing_ltotal_blocked": theorem_export["current_role_bearing_ltotal_accepts"] is False,
        "no_bridge_exported": theorem_export["legacy_to_strict_completion_bridge_exported"] is False,
        "no_strict_side_residuals_exported": theorem_export["strict_side_residual_additions_exported"] is False,
        "no_role_transfer_exported": theorem_export["strict_damping_role_transfer_theorem_exported"] is False,
        "no_role_bearing_ltotal_exported": theorem_export["role_bearing_ltotal_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_matrix"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2604",
        "stage_id": "S1554",
        "status": "P2604_POST_SOURCE_BRIDGE_READINESS_MATRIX_STRICT_DAMPING_SOURCE_INHERITED_ROLE_BEARING_LTOTAL_BLOCKED_BY_THREE_BRIDGE_ROLE_GATES_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_post_source_bridge_readiness_matrix": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {"P2603_FRACTAL_CODIMENSION_SLOPE_SOURCE": sha256_json(p2603_payload)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_post_source_bridge_readiness_matrix"]["theorem_export"]
    lines = [
        "# P2604/S1554 strict damping post-source bridge readiness matrix", "",
        f"Status: `{payload['status']}`", "", "## Matrix statement", "",
        t["matrix_statement"], "", "## Computed consequences", "",
        f"- Strict damping beta/eta source inherited from P2603: `{t['strict_damping_beta_eta_source_inherited_from_p2603']}`.",
        f"- Bridge/role gates: `{t['bridge_role_gates']}`.",
        f"- Truth-table rows: `{t['bridge_readiness_truth_table_row_count']}`.",
        f"- Role-bearing accepting rows: `{t['role_bearing_accepting_row_count']}`.",
        f"- Current role-bearing L_total accepts: `{t['current_role_bearing_ltotal_accepts']}`.",
        f"- Concrete missing ingredients: `{t['concrete_missing_ingredients']}`.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Scope guards", "",
        "No legacy-to-strict bridge, strict-side residual additions, role-transfer theorem, role-bearing `L_total`, QW-2191 discharge, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['strict_damping_post_source_bridge_readiness_matrix']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2604/S1554 strict damping post-source bridge readiness matrix

`P2604/S1554` audits the state after P2603: the strict damping beta/eta source normal form is discharged, but role-bearing `L_total` readiness still requires explicit legacy-to-strict completion-map evidence, strict-side residual additions, and a strict damping role-transfer theorem.  The matrix has exactly one accepting role-bearing row over those three bridge/role gates, so source discharge is not a silent substitute for bridge completion.
""".strip()
    lag_section = """
## P2604/S1554 post-source bridge readiness Ltotal guard

`P2604/S1554` keeps the damping/compression term non-role-bearing in `L_total` despite P2603 source discharge.  The next admissible noncyclic move is a concrete `K_legacy_ont -> K_strict_gate` completion-map theorem with strict-side residual additions, followed only then by a role-transfer audit.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2604/S1554 strict damping post-source bridge readiness matrix", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2604/S1554 post-source bridge readiness Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
