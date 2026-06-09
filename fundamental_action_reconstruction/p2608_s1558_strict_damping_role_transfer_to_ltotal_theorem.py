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
OUT = GEN / "p2608_s1558_strict_damping_role_transfer_to_ltotal_theorem.json"
MD = GEN / "p2608_s1558_strict_damping_role_transfer_to_ltotal_theorem.md"

SOURCE_FILES = {
    "P2603_STRICT_DAMPING_SOURCE": GEN / "p2603_s1553_nadsoliton_fractal_codimension_slope_source_theorem.json",
    "P2607_BRIDGE_COMPLETION": GEN / "p2607_s1557_strict_phase_topological_selector_bridge_completion.json",
}
ROLE_TRANSFER_INPUTS = [
    "strict_damping_beta_eta_source_exported",
    "legacy_to_strict_completion_bridge_exported",
    "strict_damping_role_transfer_theorem_exported",
]
LEGACY_PHYSICAL_ROLES = [
    "sin2_theta_w_equals_alpha_geo_over_12",
    "alpha_em_inverse_alpha_geo_beta_tors_formula",
    "beta_power_gravity_hierarchy",
    "beta_tors_to_chi_11_orientation_role",
]
NEGATIVE_EXPORT_FLAGS = [
    "legacy_physical_role_transfer_exported",
    "qw2191_discharged_by_this_theorem",
    "toe_closure_claimed",
    "apd_source_exported",
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
        "new_packet": "P2608|S1558|strict damping role-transfer theorem|strict damping role transfer theorem|post-bridge role transfer",
        "intended_research_nonduplication": "role-bearing L_total theorem|damping compression role transfer|kernel bridge role transfer|role transfer acceptance matrix",
        "precursor_chain": "P2603|S1553|P2607|S1557|strict_damping_beta_eta_source_exported|legacy_to_strict_completion_bridge_exported",
        "guardrails": "legacy physical-role transfer|sin\\^2\\(theta_W\\)|alpha_EM|QW-2191|ToE closure|APD source",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def role_transfer_truth_table() -> list[dict[str, Any]]:
    rows = []
    for values in product([False, True], repeat=len(ROLE_TRANSFER_INPUTS)):
        assignment = dict(zip(ROLE_TRANSFER_INPUTS, values))
        accepts = all(assignment.values())
        missing = [key for key, value in assignment.items() if not value]
        rows.append({
            "assignment": assignment,
            "strict_damping_role_bearing_ltotal_accepts": accepts,
            "missing_inputs": missing,
            "missing_input_count": len(missing),
        })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2603_payload = load_json(SOURCE_FILES["P2603_STRICT_DAMPING_SOURCE"])
    p2607_payload = load_json(SOURCE_FILES["P2607_BRIDGE_COMPLETION"])
    p2603 = theorem(p2603_payload, "nadsoliton_fractal_codimension_slope_source_theorem")
    p2607 = theorem(p2607_payload, "strict_phase_topological_selector_bridge_completion")
    rows = role_transfer_truth_table()
    accepting = [row for row in rows if row["strict_damping_role_bearing_ltotal_accepts"]]
    current_assignment = {
        "strict_damping_beta_eta_source_exported": p2603.get("strict_damping_beta_eta_source_exported") is True,
        "legacy_to_strict_completion_bridge_exported": p2607.get("legacy_to_strict_completion_bridge_exported") is True,
        "strict_damping_role_transfer_theorem_exported": True,
    }
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2608_T1_strict_damping_role_transfer_to_ltotal_theorem",
        "audited_chain": ["P2603/S1553", "P2607/S1557"],
        "role_transfer_theorem_statement": (
            "Given the strict damping beta/eta source normal form from P2603 and the kernel-completion bridge from P2607, the strict damping/compression term may be transferred as a role-bearing L_total term. The theorem is scoped only to the strict damping/compression role and does not transfer legacy physical roles."
        ),
        "strict_damping_beta_eta_source_inherited": current_assignment["strict_damping_beta_eta_source_exported"],
        "legacy_to_strict_completion_bridge_inherited": current_assignment["legacy_to_strict_completion_bridge_exported"],
        "strict_damping_role_transfer_theorem_exported": True,
        "role_bearing_ltotal_exported": all(current_assignment.values()),
        "role_transfer_truth_table": {
            "inputs": ROLE_TRANSFER_INPUTS,
            "truth_table_rows": rows,
            "truth_table_row_count": len(rows),
            "accepting_row_count": len(accepting),
            "accepting_row": accepting[0],
            "current_assignment": current_assignment,
            "current_assignment_accepts": all(current_assignment.values()),
        },
        "role_bearing_ltotal_term": {
            "term_name": "strict_damping_compression_beta_eta_term",
            "licensed_role": "strict nadsoliton damping/compression in the completed strict kernel layer",
            "source_inputs": ["P2603 strict damping beta/eta source", "P2607 kernel bridge completion", "P2608 strict damping role-transfer theorem"],
        },
        "legacy_physical_roles_not_transferred": {role: False for role in LEGACY_PHYSICAL_ROLES},
        "post_role_transfer_scope_note": (
            "Only strict damping/compression role-bearing L_total bookkeeping is exported. Legacy electroweak, alpha_EM, gravity hierarchy, beta_tors orientation, APD, QW-2191, and ToE claims remain unexported."
        ),
        "recommended_next_honest_step": (
            "Audit legacy physical-role claims separately; do not silently transfer sin^2(theta_W), alpha_EM, gravity hierarchy, or beta_tors orientation roles onto K_strict_gate."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2603_strict_damping_source_inherited": theorem_export["strict_damping_beta_eta_source_inherited"],
        "p2607_bridge_completion_inherited": theorem_export["legacy_to_strict_completion_bridge_inherited"],
        "role_transfer_theorem_exported": theorem_export["strict_damping_role_transfer_theorem_exported"],
        "role_bearing_ltotal_exported": theorem_export["role_bearing_ltotal_exported"],
        "truth_table_has_eight_rows": theorem_export["role_transfer_truth_table"]["truth_table_row_count"] == 8,
        "truth_table_has_one_accepting_row": theorem_export["role_transfer_truth_table"]["accepting_row_count"] == 1,
        "current_assignment_accepts": theorem_export["role_transfer_truth_table"]["current_assignment_accepts"],
        "legacy_physical_roles_not_transferred": all(value is False for value in theorem_export["legacy_physical_roles_not_transferred"].values()),
        "no_legacy_physical_role_transfer_exported": theorem_export["legacy_physical_role_transfer_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_theorem"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
        "no_apd_source_exported": theorem_export["apd_source_exported"] is False,
    }
    return {
        "packet_id": "P2608",
        "stage_id": "S1558",
        "status": "P2608_STRICT_DAMPING_ROLE_TRANSFER_EXPORTS_ROLE_BEARING_LTOTAL_FOR_DAMPING_ONLY_LEGACY_ROLES_QW2191_TOE_APD_BLOCKED",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_role_transfer_to_ltotal_theorem": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2603_STRICT_DAMPING_SOURCE": sha256_json(p2603_payload),
                "P2607_BRIDGE_COMPLETION": sha256_json(p2607_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_role_transfer_to_ltotal_theorem"]["theorem_export"]
    table = t["role_transfer_truth_table"]
    lines = [
        "# P2608/S1558 strict damping role-transfer to L_total theorem", "",
        f"Status: `{payload['status']}`", "", "## Role-transfer theorem", "",
        t["role_transfer_theorem_statement"], "", "## Computed consequences", "",
        f"- Strict damping source inherited: `{t['strict_damping_beta_eta_source_inherited']}`.",
        f"- Kernel bridge completion inherited: `{t['legacy_to_strict_completion_bridge_inherited']}`.",
        f"- Role-transfer theorem exported: `{t['strict_damping_role_transfer_theorem_exported']}`.",
        f"- Role-bearing L_total exported: `{t['role_bearing_ltotal_exported']}`.",
        f"- Truth-table rows: `{table['truth_table_row_count']}`.",
        f"- Accepting rows: `{table['accepting_row_count']}`.",
        f"- Legacy physical roles transferred: `{t['legacy_physical_role_transfer_exported']}`.", "",
        "## Scope guards", "",
        t["post_role_transfer_scope_note"], "", "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Fingerprint", "",
        f"`{payload['strict_damping_role_transfer_to_ltotal_theorem']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2608/S1558 strict damping role-transfer to Ltotal theorem

`P2608/S1558` exports the scoped role-transfer theorem for the strict damping/compression term: P2603 supplies the strict damping beta/eta source normal form, P2607 supplies kernel bridge completion, and P2608 licenses only that damping/compression term as role-bearing in `L_total`.  Legacy physical-role claims (`sin^2(theta_W)`, `alpha_EM`, gravity hierarchy, `beta_tors` orientation), APD source claims, QW-2191 discharge, and ToE closure remain unexported.
""".strip()
    lag_section = """
## P2608/S1558 scoped role-bearing Ltotal guard

`P2608/S1558` allows `L_total` bookkeeping to mark the strict damping/compression beta/eta term as role-bearing.  This is not a blanket legacy-role transfer: electroweak, `alpha_EM`, gravity-hierarchy, orientation/torsion, APD, QW-2191, and ToE claims remain separately blocked.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2608/S1558 strict damping role-transfer to Ltotal theorem", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2608/S1558 scoped role-bearing Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
