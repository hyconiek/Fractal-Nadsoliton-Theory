#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

P475_SUMMARY = (
    GENERATED / "p475_current_strict_projective_only_continuation_decision_packet_summary.json"
)
P630_SUMMARY = (
    GENERATED
    / "p630_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_r83_vacuum_eom_yukawa_elimination_packet_summary.json"
)
P434_SUMMARY = (
    GENERATED
    / "p434_current_strict_canonical_local_diagonal_mode2_defect_value_instantiation_evaluation_probe_summary.json"
)
P477_SUMMARY = (
    GENERATED
    / "p477_current_strict_r18_pair1_residual_zero_equations_value_instantiation_probe_summary.json"
)

OUT = (
    GENERATED
    / "p631_current_strict_direct_formal_c1s1_route_negative_freeze_decision_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p631_current_strict_direct_formal_c1s1_route_negative_freeze_decision_packet_summary.json"
)


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _p(path: Path) -> str:
    try:
        return str(path.relative_to(REPO))
    except ValueError:
        return str(path)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    required = {
        "P475_summary": P475_SUMMARY,
        "P630_summary": P630_SUMMARY,
        "P434_summary": P434_SUMMARY,
    }
    missing_required = [k for k, p in required.items() if not p.is_file()]
    if missing_required:
        payload = {
            "stage": "P631",
            "date": datetime.now(timezone.utc).date().isoformat(),
            "goal": "declare_professorial_negative_freeze_for_direct_formal_c1s1_residual_cancellation_route_under_t166_nonzero_decision",
            "status": "FAIL_MISSING_REQUIRED_INPUTS",
            "missing_required_inputs": missing_required,
            "required_paths": {k: _p(v) for k, v in required.items()},
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(
            json.dumps(
                {
                    "stage": payload["stage"],
                    "status": payload["status"],
                    "decision": None,
                    "recommended_next_strict_target": None,
                    "no_false_pass": True,
                },
                indent=2,
                ensure_ascii=True,
            )
            + "\n",
            encoding="ascii",
        )
        print(OUT_SUMMARY)
        return

    p475 = _read_json(P475_SUMMARY)
    projective_only = (p475.get("selected_continuation") == "projective_only") or (
        str(p475.get("decision") or "") == "PROJECTIVE_ONLY_CONTINUATION_SELECTED"
    )
    if not projective_only:
        payload = {
            "stage": "P631",
            "date": datetime.now(timezone.utc).date().isoformat(),
            "goal": "declare_professorial_negative_freeze_for_direct_formal_c1s1_residual_cancellation_route_under_t166_nonzero_decision",
            "status": "FAIL_PROJECTIVE_ONLY_CONTINUATION_NOT_SELECTED",
            "p475_summary": _p(P475_SUMMARY),
            "p475_selected_continuation": p475.get("selected_continuation"),
            "decision": None,
            "no_false_pass": True,
            "hard_limits": [
                "Decision applies only under projective-only continuation; it does not claim directed sign lift.",
            ],
        }
        OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(
            json.dumps(
                {
                    "stage": payload["stage"],
                    "status": payload["status"],
                    "decision": None,
                    "recommended_next_strict_target": "H37",
                    "no_false_pass": True,
                },
                indent=2,
                ensure_ascii=True,
            )
            + "\n",
            encoding="ascii",
        )
        print(OUT_SUMMARY)
        return

    p630 = _read_json(P630_SUMMARY)
    remaining = list(p630.get("remaining_missing_upstream_objects") or [])
    expected = {
        "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation",
        "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    }
    if not expected.issubset(set(remaining)):
        payload = {
            "stage": "P631",
            "date": datetime.now(timezone.utc).date().isoformat(),
            "goal": "declare_professorial_negative_freeze_for_direct_formal_c1s1_residual_cancellation_route_under_t166_nonzero_decision",
            "status": "FAIL_P630_REMAINING_MISSING_LIST_NOT_AT_EXPECTED_FRONTIER",
            "p630_summary": _p(P630_SUMMARY),
            "p630_status": p630.get("status"),
            "p630_remaining_missing_upstream_objects": remaining,
            "expected_subset": sorted(expected),
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(
            json.dumps(
                {
                    "stage": payload["stage"],
                    "status": payload["status"],
                    "decision": None,
                    "recommended_next_strict_target": "P630",
                    "no_false_pass": True,
                },
                indent=2,
                ensure_ascii=True,
            )
            + "\n",
            encoding="ascii",
        )
        print(OUT_SUMMARY)
        return

    p434 = _read_json(P434_SUMMARY)
    cuts_o2 = bool(p434.get("cuts_O2_on_pair1_by_N466"))
    f2_abs = p434.get("F2_abs")

    if not cuts_o2:
        payload = {
            "stage": "P631",
            "date": datetime.now(timezone.utc).date().isoformat(),
            "goal": "declare_professorial_negative_freeze_for_direct_formal_c1s1_residual_cancellation_route_under_t166_nonzero_decision",
            "status": "FAIL_T166_NONZERO_DECISION_NOT_PRESENT",
            "p434_summary": _p(P434_SUMMARY),
            "p434_status": p434.get("status"),
            "cuts_O2_on_pair1_by_N466": cuts_o2,
            "F2_abs": f2_abs,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(
            json.dumps(
                {
                    "stage": payload["stage"],
                    "status": payload["status"],
                    "decision": None,
                    "recommended_next_strict_target": "P630",
                    "no_false_pass": True,
                },
                indent=2,
                ensure_ascii=True,
            )
            + "\n",
            encoding="ascii",
        )
        print(OUT_SUMMARY)
        return

    evidence: dict[str, Any] = {
        "P630": {"summary_path": _p(P630_SUMMARY), "status": p630.get("status"), "remaining_missing": remaining},
        "P434": {
            "summary_path": _p(P434_SUMMARY),
            "status": p434.get("status"),
            "cuts_O2_on_pair1_by_N466": cuts_o2,
            "F2_abs": f2_abs,
            "theta_star_by_N468": p434.get("theta_star_by_N468"),
        },
        "theorem_links": [
            "N473 Corollary 1: (pair1 residual block vanishes) => (F2(d)=0) => (no diagonal/local O(2) cut on pair1).",
            "N482: packages the current value-instantiation decision |F2(d)| != 0 (T166) from P434.",
        ],
    }

    if P477_SUMMARY.exists():
        try:
            p477 = _read_json(P477_SUMMARY)
            evidence["P477"] = {
                "summary_path": _p(P477_SUMMARY),
                "status": p477.get("status"),
                "violated_equations": p477.get("violated_equations"),
                "all_zero_equations_satisfied": p477.get("all_zero_equations_satisfied"),
            }
        except Exception:
            evidence["P477"] = {"summary_path": _p(P477_SUMMARY), "parse_error": True}
    else:
        evidence["P477"] = {"summary_path": _p(P477_SUMMARY), "missing": True}

    status = "PASS_DECISION_DECLARED_DIRECT_FORMAL_C1S1_ROUTE_NEGATIVE_FREEZE_SELECTED"
    decision = "DIRECT_FORMAL_C1S1_ROUTE_NEGATIVE_FREEZE_SELECTED"

    payload = {
        "stage": "P631",
        "date": datetime.now(timezone.utc).date().isoformat(),
        "goal": "declare_professorial_negative_freeze_for_direct_formal_c1s1_residual_cancellation_route_under_t166_nonzero_decision",
        "status": status,
        "decision": decision,
        "decision_basis": {
            "continuation_premise": {
                "projective_only_selected": True,
                "p475_summary": _p(P475_SUMMARY),
            },
            "strict_diagonal_lane_fact": (
                "T166 is discharged by value-instantiation: P434 evaluates |F2(d)| != 0, hence the diagonal/local lane cuts O(2) on pair1 by N466 (packaged by N482)."
            ),
            "structural_incompatibility": (
                "By N473 (Corollary 1), if the declared pair1 residual block vanished (a1=b1=d1=0), then F2(d)=0 and therefore the diagonal/local O(2) cut would fail. "
                "So, on the strict branch that keeps the T166 nonzero decision active, the remaining direct-formal residual-cancellation targets (pair1 c1c1/s1s1 zero equations) cannot be treated as next strict bottlenecks without undoing that branch."
            ),
            "evidence_chain": evidence,
        },
        "continuation": {
            "selected": "post_projective_frontier_decision",
            "meaning": "Projective-only strict closure tasks on the direct-formal residual-cancellation lane are frozen negative on this branch; the next honest question is whether to attempt a directed/sign-sensitive selector state datum (H37/T171) or remain projective-only.",
        },
        "recommended_next_strict_target": {"target": "H37", "note": "post-projective directed selector state frontier (see T171)"},
        "hard_limits": [
            "no theorem-level pass",
            "no full-closure pass",
            "does not claim any strict zero witness for the declared pair1 residual c1c1/s1s1 equations",
            "does not claim any discharge of QW-2191",
            "does not claim H37/T171 is discharged",
            "does not claim ToE closure",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": payload["stage"],
        "status": payload["status"],
        "decision": payload["decision"],
        "recommended_next_strict_target": "H37",
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

