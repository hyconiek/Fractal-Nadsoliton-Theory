#!/usr/bin/env python3
"""Scratch probe: legacy physical-role transfer pre-audit certificate.

The bridge-completion reports now give finite and symbolic witnesses for the
kernel-comparison ansatz.  This probe performs the next mandatory guardrail
step without overclaiming it: a claim-by-claim pre-audit of legacy physical-role
formulas that would require a separate role-transfer theorem before they could
be attached to K_strict_gate.

It is computational in a narrow sense: each legacy role is encoded as a Boolean
dependency row over GF(2), the dependency matrix rank is computed, and every
claim is classified as blocked rather than transferred.  The point is to make
clear which missing theorem each legacy role would need after bridge completion.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_legacy_role_transfer_pre_audit_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_legacy_role_transfer_pre_audit_certificate_report.md"
S2_SPEC = ROOT / "fundamental_action_reconstruction" / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md"
T15_SPEC = ROOT / "fundamental_action_reconstruction" / "T15_LEGACY_TO_STRICT_KERNEL_BRIDGE_THEOREM_SPEC.md"

SOURCE_REPORTS = {
    "legacy_bridge_guardrail": HERE / "bridge_strict_completion_legacy_kernel_intermediate_bridge_guardrail_certificate_report.json",
    "component_gap_matrix": HERE / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_report.json",
    "symbolic_cancellation": HERE / "bridge_strict_completion_legacy_to_strict_symbolic_cancellation_certificate_report.json",
    "amplitude_scalar_normalization": HERE / "bridge_strict_completion_legacy_to_strict_amplitude_scalar_normalization_certificate_report.json",
    "damping_parameter_identifiability": HERE / "bridge_strict_completion_legacy_to_strict_damping_parameter_identifiability_certificate_report.json",
    "anchor_h1_classification": HERE / "bridge_strict_completion_anchor_h1_generator_classification_certificate_report.json",
}

MATRIX_COLUMNS = [
    "depends_alpha_geo",
    "depends_beta_tors",
    "depends_beta_power_hierarchy",
    "depends_chi11_or_orientation",
    "strict_successor_theorem_available_now",
    "unchanged_transfer_allowed_now",
    "modified_transfer_allowed_now",
    "rejected_by_current_bridge_data",
]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def gf2_rank(rows: list[list[int]]) -> int:
    work = [row[:] for row in rows if any(row)]
    if not work:
        return 0
    rank = 0
    col_count = len(work[0])
    for col in range(col_count):
        pivot = next((idx for idx in range(rank, len(work)) if work[idx][col] & 1), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for idx, row in enumerate(work):
            if idx != rank and row[col] & 1:
                work[idx] = [a ^ b for a, b in zip(row, work[rank])]
        rank += 1
        if rank == len(work):
            break
    return rank


def build_role_rows() -> list[dict[str, Any]]:
    rows = [
        {
            "role_id": "legacy_weak_mixing_angle",
            "legacy_formula": "sin^2(theta_W)=alpha_geo/12",
            "depends_alpha_geo": 1,
            "depends_beta_tors": 0,
            "depends_beta_power_hierarchy": 0,
            "depends_chi11_or_orientation": 0,
            "needed_transfer_theorem": "strict alpha_geo physical-role theorem after scalar normalization",
            "current_verdict": "blocked_not_transferred",
            "reason": "The bridge cancels alpha_geo as a kernel amplitude scalar; no theorem says the legacy electroweak role survives that normalization.",
        },
        {
            "role_id": "legacy_inverse_alpha_em",
            "legacy_formula": "alpha_EM^-1=alpha_geo/(2*beta_tors)*(1-beta_tors)",
            "depends_alpha_geo": 1,
            "depends_beta_tors": 1,
            "depends_beta_power_hierarchy": 0,
            "depends_chi11_or_orientation": 0,
            "needed_transfer_theorem": "joint strict alpha_geo and beta_tors role theorem",
            "current_verdict": "blocked_not_transferred",
            "reason": "The current bridge separates legacy linear beta_tors from strict nonlinear beta*d^eta and exports no beta_tors role theorem.",
        },
        {
            "role_id": "legacy_beta_power_gravity_hierarchy",
            "legacy_formula": "beta^N gravity hierarchy / beta_tors-power hierarchy",
            "depends_alpha_geo": 0,
            "depends_beta_tors": 1,
            "depends_beta_power_hierarchy": 1,
            "depends_chi11_or_orientation": 0,
            "needed_transfer_theorem": "strict successor hierarchy theorem across nonlinear damping/compression",
            "current_verdict": "blocked_not_transferred",
            "reason": "No report maps the legacy beta-power hierarchy through the strict beta=1, eta=9/5 nonlinear compression data.",
        },
        {
            "role_id": "legacy_torsion_to_chi11_orientation",
            "legacy_formula": "candidate beta_tors -> chi_11 orientation/torsion role",
            "depends_alpha_geo": 0,
            "depends_beta_tors": 1,
            "depends_beta_power_hierarchy": 0,
            "depends_chi11_or_orientation": 1,
            "needed_transfer_theorem": "selector/source theorem beta_tors -> chi_11 or replacement source",
            "current_verdict": "blocked_not_transferred",
            "reason": "GF(2)/H1 reports locate the bit but do not source it, and the anchor/H1 audit keeps selector source open.",
        },
    ]
    for row in rows:
        row["strict_successor_theorem_available_now"] = 0
        row["unchanged_transfer_allowed_now"] = 0
        row["modified_transfer_allowed_now"] = 0
        row["rejected_by_current_bridge_data"] = 0
        row["role_transfer_allowed_now"] = False
    return rows


def build_payload() -> dict[str, Any]:
    loaded = {name: load_json(path) for name, path in SOURCE_REPORTS.items()}
    guardrail = loaded["legacy_bridge_guardrail"]
    component_gap = loaded["component_gap_matrix"]
    symbolic = loaded["symbolic_cancellation"]
    amplitude = loaded["amplitude_scalar_normalization"]
    damping = loaded["damping_parameter_identifiability"]
    anchor = loaded["anchor_h1_classification"]
    s2_text = S2_SPEC.read_text(encoding="utf-8")
    t15_text = T15_SPEC.read_text(encoding="utf-8")

    role_rows = build_role_rows()
    matrix_rows = [[row[column] for column in MATRIX_COLUMNS] for row in role_rows]
    rank = gf2_rank(matrix_rows)
    blocked_count = sum(row["current_verdict"] == "blocked_not_transferred" for row in role_rows)

    summary = {
        "role_claim_count": len(role_rows),
        "dependency_matrix_columns": MATRIX_COLUMNS,
        "dependency_matrix_rank_gf2": rank,
        "all_roles_blocked_now": blocked_count == len(role_rows),
        "roles_transferred_now": sum(row["role_transfer_allowed_now"] for row in role_rows),
        "s2_lists_required_role_transfer_claims": all(term in s2_text for term in ["sin^2(theta_W)=alpha_geo/12", "alpha_EM^-1 = alpha_geo/(2*beta_tors)*(1-beta_tors)", "beta^N", "beta_tors -> chi_11"]),
        "t15_keeps_role_transfer_separate": "legacy physical-role transfer" in t15_text and "separate" in t15_text,
        "guardrail_requires_role_transfer_audit": guardrail["legacy_kernel_intermediate_bridge_summary"]["role_transfer_audit_required_after_full_bridge"],
        "symbolic_bridge_is_formula_only": not symbolic["symbolic_cancellation_summary"]["legacy_role_transfer_exported"],
        "component_gap_blocks_role_transfer": component_gap["completion_gap_summary"]["role_transfer_blocked_until_full_bridge"],
        "alpha_role_source_missing": not amplitude["amplitude_scalar_normalization_summary"]["legacy_role_transfer_allowed"],
        "beta_tors_role_source_missing": not damping["damping_parameter_identifiability_summary"]["legacy_beta_tors_to_beta_eta_theorem_exported"],
        "chi11_source_missing": anchor["classification_summary"]["selector_source_remains_open"],
        "strict_dynamic_source_exported": False,
        "role_transfer_theorem_exported": False,
        "q2191_discharged": False,
        "toe_closure_claimed": False,
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_LEGACY_ROLE_TRANSFER_PRE_AUDIT_CERTIFICATE__NO_TRANSFER_THEOREM",
        "status": "legacy-physical-role-claims-audited-all-blocked-until-separate-role-transfer-theorem",
        "source_reports": {name: rel(path) for name, path in SOURCE_REPORTS.items()},
        "source_specs": {
            "s2_priority_packet": rel(S2_SPEC),
            "t15_bridge_theorem_spec": rel(T15_SPEC),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "role-transfer",
                "legacy physical-role",
                "sin^2(theta_W)",
                "alpha_EM^-1",
                "beta^N gravity hierarchy",
                "beta_tors -> chi_11",
            ],
            "finding": "Existing guardrails require role-transfer auditing after bridge completion; this report performs a finite pre-audit and blocks every listed legacy role pending a separate theorem.",
        },
        "role_dependency_columns": MATRIX_COLUMNS,
        "role_transfer_rows": role_rows,
        "role_dependency_matrix_gf2": matrix_rows,
        "role_transfer_pre_audit_summary": summary,
        "cross_checks": {
            "source_reports_present": set(loaded) == set(SOURCE_REPORTS),
            "s2_and_t15_require_separate_role_transfer": summary["s2_lists_required_role_transfer_claims"] and summary["t15_keeps_role_transfer_separate"],
            "all_roles_blocked_with_nonzero_dependency_rank": summary["all_roles_blocked_now"] and summary["dependency_matrix_rank_gf2"] >= 3,
            "component_and_symbolic_reports_do_not_transfer_roles": summary["component_gap_blocks_role_transfer"] and summary["symbolic_bridge_is_formula_only"] and summary["guardrail_requires_role_transfer_audit"],
            "alpha_beta_chi11_sources_missing": summary["alpha_role_source_missing"] and summary["beta_tors_role_source_missing"] and summary["chi11_source_missing"],
            "closure_limits_preserved": not summary["strict_dynamic_source_exported"] and not summary["role_transfer_theorem_exported"] and not summary["q2191_discharged"] and not summary["toe_closure_claimed"],
        },
        "proof_certificate": {
            "grep_step": "rg was used to locate existing role-transfer guardrails and avoid treating the symbolic bridge as a physical-role theorem.",
            "matrix_step": f"The four audited legacy roles give a GF(2) dependency matrix of rank {rank} over columns {MATRIX_COLUMNS}.",
            "alpha_step": "sin^2(theta_W)=alpha_geo/12 is blocked because alpha_geo is only normalized as a kernel scalar, with no strict physical-role theorem.",
            "beta_step": "alpha_EM^-1 and beta^N hierarchy claims are blocked because beta_tors is not mapped to a strict physical role through beta=1, eta=9/5 compression.",
            "chi11_step": "beta_tors -> chi_11 is blocked because the anchor/H1 certificates locate the bit but keep its selector/source theorem open.",
            "scope_step": "This is a pre-audit: it classifies transfer obligations but exports zero role-transfer permissions and no QW-2191 or ToE closure.",
        },
        "hard_limits": [
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No alpha_geo electroweak role theorem is exported.",
            "No beta_tors electromagnetic/gravity hierarchy role theorem is exported.",
            "No beta_tors -> chi_11 theorem is exported.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Legacy physical-role transfer pre-audit certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "The listed legacy physical-role claims are dependency-audited and all remain",
        "blocked until a separate role-transfer theorem is supplied.  The symbolic",
        "bridge/cancellation certificates are kernel-comparison evidence only.",
        "",
        "## Summary",
        "",
    ]
    for key, value in payload["role_transfer_pre_audit_summary"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Role rows", ""])
    for row in payload["role_transfer_rows"]:
        lines.append(f"- `{row['role_id']}`: `{row['current_verdict']}` — {row['reason']}")
    lines.extend(["", "## Cross-checks", ""])
    for key, value in payload["cross_checks"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Proof certificate", ""])
    for key, value in payload["proof_certificate"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Hard limits", ""])
    for limit in payload["hard_limits"]:
        lines.append(f"- {limit}")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(payload)
    print(json.dumps(payload["role_transfer_pre_audit_summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
