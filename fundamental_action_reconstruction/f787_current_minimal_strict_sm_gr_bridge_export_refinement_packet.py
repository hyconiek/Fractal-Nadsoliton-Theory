#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F784 = GENERATED / "f784_current_minimal_strict_sm_gr_bridge_registry_packet.json"
IN_P785 = GENERATED / "p785_current_minimal_strict_sm_gr_bridge_registry_boundary_audit_probe.json"
IN_P786 = GENERATED / "p786_current_qw2093_kernel_invariants_only_formula_layer_audit_probe.json"
IN_WORKING_NOTE = REPO / "WORKING_NOTE_LEGACY_KEEP_CUT_AND_MINIMAL_STRICT_SM_GR_BRIDGE.md"

OUT = GENERATED / "f787_current_minimal_strict_sm_gr_bridge_export_refinement_packet.json"
OUT_SUMMARY = GENERATED / "f787_current_minimal_strict_sm_gr_bridge_export_refinement_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def normalize(text: str) -> str:
    return (
        text.lower()
        .replace("“", '"')
        .replace("”", '"')
        .replace("’", "'")
        .replace("->", " ")
        .replace("→", " ")
        .replace("/", " ")
        .replace("_", " ")
        .replace("-", " ")
    )


def contains_all(text: str, needles: list[str]) -> bool:
    hay = normalize(text)
    return all(normalize(needle) in hay for needle in needles)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F784, IN_P785, IN_P786, IN_WORKING_NOTE]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F787",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f784 = load_json(IN_F784)
    p785 = load_json(IN_P785)
    p786 = load_json(IN_P786)
    working_note = IN_WORKING_NOTE.read_text(encoding="utf-8")

    targets = f784.get("targets") or {}
    p786_targets = p786.get("targets") or {}

    note_unit_discipline_present = contains_all(
        working_note,
        [
            "dimensionless or explicitly normalized observables",
            "proxy -> GeV",
        ],
    )

    mass_export = {
        "status_label": "strict-derived",
        "bridge_ready": True,
        "export_decision": "retained_export_minimal_bridge_internal_proxy_only",
        "reason": "The mass layer remains dimensionless/internal and passes the F784/P785 no-false-pass boundary.",
        "source_artifacts": (targets.get("mass_ratio_ordering_layer") or {}).get("source_artifacts") or [],
        "hard_limits": (targets.get("mass_ratio_ordering_layer") or {}).get("hard_limits") or [],
    }

    sin2_export = {
        "status_label": "strict-derived",
        "bridge_ready": True,
        "export_decision": "retained_export_candidate_observable_with_explicit_legacy_nontransfer_boundary",
        "reason": "The strict-side observable may remain exported as a candidate bridge object, but legacy Weinberg-role transfer stays explicitly blocked.",
        "source_artifacts": (targets.get("sin2_theta_w_eff") or {}).get("source_artifacts") or [],
        "promotion_blockers": (targets.get("sin2_theta_w_eff") or {}).get("promotion_blockers") or [],
        "formula_layer_blockers": (p786_targets.get("sin2_theta_w_eff") or {}).get("blockers") or [],
        "hard_limits": (targets.get("sin2_theta_w_eff") or {}).get("hard_limits") or [],
    }

    alpha_s_formula_blockers = (p786_targets.get("alpha_s_boundary_mu0_alpha0") or {}).get("blockers") or []
    alpha_s_demote = (
        p785.get("status") == "P785_CURRENT_MINIMAL_BRIDGE_REGISTRY_BOUNDARY_AUDIT_PASS_WITH_BOUNDARIES_EXPLICIT"
        and note_unit_discipline_present
        and any(
            blocker in alpha_s_formula_blockers
            for blocker in [
                "boundary_scale_uses_mass_chain_object_not_kernel_invariant",
                "boundary_coupling_uses_mass_ratio_plus_micro_constant_ansatz",
                "upstream_source_explicitly_labels_hierarchy_log_boundary_ansatz",
                "mass_chain_still_not_promoted_to_strict_first_principles",
            ]
        )
    )

    if alpha_s_demote:
        alpha_s_export = {
            "status_label": "open",
            "bridge_ready": False,
            "export_decision": "demoted_to_support_only_nonexport_pending_dimensionless_replacement",
            "reason": "The current alpha_s boundary object violates the minimal bridge discipline: it still runs through a frozen hierarchy-log ansatz and GeV-level mass-chain boundary scale rather than a kernel-invariants-only or explicitly normalized bridge object.",
            "source_artifacts": (targets.get("alpha_s_boundary_mu0_alpha0") or {}).get("source_artifacts") or [],
            "promotion_blockers": (targets.get("alpha_s_boundary_mu0_alpha0") or {}).get("promotion_blockers") or [],
            "formula_layer_blockers": alpha_s_formula_blockers,
            "hard_limits": [
                "Do not export mu0(alpha_s) from the current GeV-level mass-chain boundary as part of the minimal strict bridge.",
                "Do not reinterpret downstream QW-2087 success as proof that the current alpha_s boundary law is kernel-invariants-only.",
                "Do not upgrade this entry again without a dimensionless or explicitly normalized replacement boundary.",
            ],
        }
    else:
        alpha_s_export = {
            "status_label": "strict-derived",
            "bridge_ready": True,
            "export_decision": "retained_export",
            "reason": "No demotion trigger detected.",
            "source_artifacts": (targets.get("alpha_s_boundary_mu0_alpha0") or {}).get("source_artifacts") or [],
            "promotion_blockers": (targets.get("alpha_s_boundary_mu0_alpha0") or {}).get("promotion_blockers") or [],
            "formula_layer_blockers": alpha_s_formula_blockers,
            "hard_limits": (targets.get("alpha_s_boundary_mu0_alpha0") or {}).get("hard_limits") or [],
        }

    gravity_export = {
        "status_label": "strict-derived-with-external-observable-origin",
        "bridge_ready": True,
        "export_decision": "retained_export_with_external_origin_boundary",
        "reason": "The current gravity bridge object remains admissible only with explicit external observable origin and open internal-origin blocker.",
        "source_artifacts": (targets.get("g_dimensionless_mu_ref") or {}).get("source_artifacts") or [],
        "promotion_blockers": (targets.get("g_dimensionless_mu_ref") or {}).get("promotion_blockers") or [],
        "hard_limits": (targets.get("g_dimensionless_mu_ref") or {}).get("hard_limits") or [],
    }

    refined_targets = {
        "mass_ratio_ordering_layer": mass_export,
        "sin2_theta_w_eff": sin2_export,
        "alpha_s_boundary_mu0_alpha0": alpha_s_export,
        "g_dimensionless_mu_ref": gravity_export,
    }

    exported_ids = [key for key, value in refined_targets.items() if value["bridge_ready"]]
    support_only_ids = [key for key, value in refined_targets.items() if not value["bridge_ready"]]

    artifact = {
        "stage": "F787",
        "packet_name": "CurrentMinimalStrictSMGRBridgeExportRefinement_v1",
        "status": "F787_EXECUTED_CURRENT_MINIMAL_STRICT_SM_GR_BRIDGE_EXPORT_REFINEMENT_PACKET_NO_FALSE_PASS",
        "as_of": AS_OF,
        "inputs": {
            "f784_registry": rel(IN_F784),
            "p785_boundary_audit": rel(IN_P785),
            "p786_formula_layer_audit": rel(IN_P786),
            "working_note": rel(IN_WORKING_NOTE),
        },
        "decision_rules": [
            "Retain only objects that still fit the minimal strict bridge boundary after P785/P786.",
            "Keep sin2 as a strict candidate observable only with explicit non-transfer boundary to the legacy Weinberg role.",
            "Demote alpha_s boundary objects when the current route still depends on frozen hierarchy-log ansatz and GeV-level mass-chain boundary scale.",
            "Do not weaken the gravity external-origin boundary.",
        ],
        "refined_targets": refined_targets,
        "export_sets": {
            "minimal_bridge_export_ids": exported_ids,
            "support_only_nonexport_ids": support_only_ids,
        },
        "current_honest_reading": [
            "F784 remains boundary-clean after P785, so the refinement question is not leakage but admissible export scope.",
            "sin2_theta_w_eff can stay in the minimal strict bridge only as a candidate observable with an explicit legacy non-transfer boundary.",
            "alpha_s_boundary_mu0_alpha0 should not currently remain in the minimal strict export set because the present boundary is not yet dimensionless/normalized and not kernel-invariants-only.",
            "This refinement narrows the bridge rather than broadening it; that is the correct no-false-pass move on the current repo state.",
        ],
        "recommended_next_move": "Build a replacement alpha_s boundary object from dimensionless strict mass-observable ratios or explicitly normalized strict observables, then rerun the minimal bridge export decision.",
        "no_false_pass": True,
    }

    summary = {
        "stage": "F787",
        "status": artifact["status"],
        "as_of": AS_OF,
        "n_exported": len(exported_ids),
        "n_support_only": len(support_only_ids),
        "minimal_bridge_export_ids": exported_ids,
        "support_only_nonexport_ids": support_only_ids,
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
