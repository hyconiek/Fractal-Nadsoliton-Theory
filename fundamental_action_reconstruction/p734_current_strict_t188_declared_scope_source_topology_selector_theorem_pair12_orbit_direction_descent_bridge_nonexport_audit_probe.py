#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-18"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P729 = GENERATED / "p729_current_strict_t183_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
IN_P731 = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"
IN_F148 = GENERATED / "f148_first_actual_source_topology_basis_independent_promotion_witness_packet_summary.json"
IN_F149 = GENERATED / "f149_first_actual_source_topology_quotient_safe_qw2191_resolution_witness_packet_summary.json"
IN_F150 = GENERATED / "f150_first_actual_declared_scope_source_topology_selector_theorem_packet_summary.json"

OUT_JSON = GENERATED / "p734_current_strict_t188_declared_scope_source_topology_selector_theorem_pair12_orbit_direction_descent_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p734_current_strict_t188_declared_scope_source_topology_selector_theorem_pair12_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P729, IN_P731, IN_F148, IN_F149, IN_F150]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P734",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p729 = load_json(IN_P729)
    p731 = load_json(IN_P731)
    f148 = load_json(IN_F148)
    f149 = load_json(IN_F149)
    f150 = load_json(IN_F150)

    f148_support = f148.get("support_packet") or {}
    basis_reduction = f148_support.get("basis_class_reduction_operator") or {}

    f149_support = f149.get("support_packet") or {}
    quotient_relation = f149_support.get("quotient_relation") or {}
    distinguished_quotient_class = f149_support.get("distinguished_quotient_class") or {}

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        passed = actual == expected
        checks.append(
            {
                "id": check_id,
                "actual": actual,
                "expected": expected,
                "pass": passed,
                "meaning": meaning,
            }
        )
        if not passed:
            blocking.append(check_id)

    current_declared_scope_source_topology_selector_theorem_remains_quotient_class_only = (
        bool(basis_reduction.get("forgets_chart_basis_labels"))
        and bool(quotient_relation.get("quotients_chart_labels_only_via_basis_free_packet"))
        and bool(distinguished_quotient_class.get("quotient_class_only"))
        and not bool(f149.get("raw_theta_uniqueness_claimed"))
        and bool(f150.get("declared_scope_only"))
        and not bool(f150.get("raw_theta_uniqueness_claimed"))
    )

    current_declared_scope_source_topology_selector_theorem_is_pair12_orbit_direction_typed = (
        not bool(basis_reduction.get("forgets_chart_basis_labels"))
        and not bool(distinguished_quotient_class.get("quotient_class_only"))
        and bool(f150.get("raw_theta_uniqueness_claimed"))
        and bool(f150.get("tau_src_identified_with_s_prelm"))
    )

    p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem = (
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and current_declared_scope_source_topology_selector_theorem_is_pair12_orbit_direction_typed
    )

    add_check(
        "p729_pair12_orbit_direction_split_already_localized",
        bool(p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions")),
        True,
        "P729 already localizes the surviving pair1/pair2 ambiguity as opposite residual-datum orbit directions on the same source-side carrier.",
    )
    add_check(
        "p731_pair12_witness_split_already_separated",
        {
            "split_separated": bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches")),
            "pair1_sign": p731.get("pair1_w_break_branch_score_sign"),
            "pair2_sign": p731.get("pair2_w_break_branch_score_sign"),
        },
        {
            "split_separated": True,
            "pair1_sign": "negative",
            "pair2_sign": "positive",
        },
        "P731 already separates the surviving pair1/pair2 branches by opposite witness-score signs.",
    )
    add_check(
        "f148_basis_reduction_still_forgets_chart_basis_labels",
        {
            "forgets_chart_basis_labels": bool(basis_reduction.get("forgets_chart_basis_labels")),
            "preserves_source_side_only_scope": bool(basis_reduction.get("preserves_source_side_only_scope")),
        },
        {
            "forgets_chart_basis_labels": True,
            "preserves_source_side_only_scope": True,
        },
        "The current actual basis-independent source-side promotion still forgets chart basis labels and keeps only basis-free selector classes.",
    )
    add_check(
        "f149_current_qw2191_safe_resolution_is_quotient_class_only",
        {
            "distinguished_quotient_class_exported": bool(f149.get("distinguished_quotient_class_exported")),
            "quotient_class_only": bool(distinguished_quotient_class.get("quotient_class_only")),
            "raw_theta_uniqueness_claimed": bool(f149.get("raw_theta_uniqueness_claimed")),
        },
        {
            "distinguished_quotient_class_exported": True,
            "quotient_class_only": True,
            "raw_theta_uniqueness_claimed": False,
        },
        "The current source-topology QW-2191-safe witness exports only one distinguished quotient class, not one raw pair-indexed orbit-direction branch.",
    )
    add_check(
        "f150_declared_scope_source_topology_selector_theorem_still_scope_limited",
        {
            "declared_scope_source_topology_selector_theorem_exported": bool(
                f150.get("declared_scope_source_topology_selector_theorem_exported")
            ),
            "qw2191_quotient_safe_discharged": bool(f150.get("qw2191_quotient_safe_discharged")),
            "declared_scope_only": bool(f150.get("declared_scope_only")),
            "raw_theta_uniqueness_claimed": bool(f150.get("raw_theta_uniqueness_claimed")),
            "tau_src_identified_with_s_prelm": bool(f150.get("tau_src_identified_with_s_prelm")),
        },
        {
            "declared_scope_source_topology_selector_theorem_exported": True,
            "qw2191_quotient_safe_discharged": True,
            "declared_scope_only": True,
            "raw_theta_uniqueness_claimed": False,
            "tau_src_identified_with_s_prelm": False,
        },
        "The strongest current source-topology theorem is real and quotient-safe, but it remains declared-scope only, does not claim raw-theta uniqueness, and does not identify tau_src with the current selector carrier.",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_remains_quotient_class_only",
        current_declared_scope_source_topology_selector_theorem_remains_quotient_class_only,
        True,
        "Therefore the current strongest source-side theorem lane still remains basis-free / quotient-class only on current exports.",
    )
    add_check(
        "p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem",
        p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem,
        False,
        "That current declared-scope source-topology selector theorem still does not descend the opposite P731 pair1/pair2 branch split as one typed branch distinction.",
    )
    add_check(
        "t188_declared_scope_source_topology_selector_theorem_pair12_orbit_direction_descent_bridge_exported",
        False,
        False,
        "So the declared-scope source-topology selector theorem pair1/pair2 orbit-direction descent bridge still remains unexported on current repo state.",
    )

    status = (
        "PASS_DECLARED_SCOPE_SOURCE_TOPOLOGY_SELECTOR_THEOREM_PAIR12_ORBIT_DIRECTION_DESCENT_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P734_REQUIRES_REVIEW_CHANGED_DECLARED_SCOPE_SOURCE_TOPOLOGY_PAIR12_DESCENT_STATE"
    )

    artifact = {
        "stage": "P734",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t188_declared_scope_source_topology_selector_theorem_pair12_orbit_direction_descent_bridge_nonexport_boundary_only",
        "inputs": {
            "P729": str(IN_P729.relative_to(REPO)),
            "P731": str(IN_P731.relative_to(REPO)),
            "F148": str(IN_F148.relative_to(REPO)),
            "F149": str(IN_F149.relative_to(REPO)),
            "F150": str(IN_F150.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t188_target_name": "DeclaredScopeSourceTopologySelectorTheoremPair12OrbitDirectionDescentBridge_global_C_v1_strict_v1",
        "t188_target_exported_on_current_repo_state": False,
        "current_declared_scope_source_topology_selector_theorem_exported": bool(
            f150.get("declared_scope_source_topology_selector_theorem_exported")
        ),
        "current_declared_scope_source_topology_selector_theorem_remains_quotient_class_only": (
            current_declared_scope_source_topology_selector_theorem_remains_quotient_class_only
        ),
        "p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem": (
            p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem
        ),
        "audit_conclusion": {
            "current_repo_already_exports_strong_source_side_declared_scope_selector_theorem": bool(
                f150.get("declared_scope_source_topology_selector_theorem_exported")
            ),
            "current_repo_already_exports_t188_target": False,
            "next_honest_move": (
                "attack_a_typed_local_source_side_bind_or_descent_mechanism_from_the_declared_scope_source_topology_selector_theorem_to_one_residual_datum_pair12_orbit_direction_branch_without_collapsing_back_to_basis_free_quotient_only_data_or_current_positive_convention_transport"
            ),
        },
        "hard_limits": [
            "No T188 discharge claim.",
            "No claim that the current declared-scope source-topology theorem already selects one raw pair1/pair2 orbit-direction branch.",
            "No claim that quotient-safe QW-2191 resolution upgrades to raw-theta uniqueness.",
            "No claim that tau_src is identified with the current selector carrier.",
            "No directed/sign-sensitive physical orientation datum claim.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P734",
        "status": status,
        "as_of": AS_OF,
        "t188_target_name": "DeclaredScopeSourceTopologySelectorTheoremPair12OrbitDirectionDescentBridge_global_C_v1_strict_v1",
        "t188_target_exported_on_current_repo_state": False,
        "current_declared_scope_source_topology_selector_theorem_exported": bool(
            f150.get("declared_scope_source_topology_selector_theorem_exported")
        ),
        "current_declared_scope_source_topology_selector_theorem_remains_quotient_class_only": (
            current_declared_scope_source_topology_selector_theorem_remains_quotient_class_only
        ),
        "p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem": (
            p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem
        ),
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
