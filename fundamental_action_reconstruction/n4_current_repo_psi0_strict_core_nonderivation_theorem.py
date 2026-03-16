#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent


def load_json(relative_path: str) -> Any:
    return json.loads((ROOT / relative_path).read_text(encoding="utf-8"))


def get_by_path(payload: Any, path: tuple[str, ...]) -> Any:
    current = payload
    for key in path:
        if not isinstance(current, dict) or key not in current:
            raise KeyError(".".join(path))
        current = current[key]
    return current


def main() -> None:
    sources = {
        "H30": load_json("generated/h30_kernel_invariant_psi0_anchor_audit.json"),
        "H31": load_json("generated/h31_psi0_to_pair1_reduction_audit.json"),
        "H34": load_json("generated/h34_basis_covariance_target_independence_audit_summary.json"),
        "H35": load_json("generated/h35_pair1_axis_selection_audit_summary.json"),
        "H36": load_json("generated/h36_directed_axis_orientation_audit.json"),
        "H37": load_json("generated/h37_sign_distinction_state_audit.json"),
        "H42": load_json("generated/h42_c_based_retardation_operator_on_pair1_audit.json"),
        "P1": load_json("generated/pair1_operator_probe_report.json"),
    }

    checks_spec = [
        {
            "id": "h30_deterministic_kernel_invariant_candidate",
            "source": "H30",
            "path": ("classification", "deterministic_from_kernel_invariants"),
            "expected": True,
            "meaning": "psi0 is deterministic from kernel invariants",
        },
        {
            "id": "h30_no_strict_core_selector_export",
            "source": "H30",
            "path": ("classification", "strict_core_selector_export"),
            "expected": False,
            "meaning": "psi0 is not exported as a strict-core selector datum",
        },
        {
            "id": "h31_coordinate_embedding_present",
            "source": "H31",
            "path": ("classification", "coordinate_level_embedding_present"),
            "expected": True,
            "meaning": "psi0 admits a coordinate embedding into pair1",
        },
        {
            "id": "h31_no_strict_core_selector_reduction",
            "source": "H31",
            "path": ("classification", "strict_core_selector_reduction_present"),
            "expected": False,
            "meaning": "the embedding is not a strict-core selector reduction",
        },
        {
            "id": "h31_pair1_not_proven_selector_target",
            "source": "H31",
            "path": ("classification", "pair1_as_selector_target_proven"),
            "expected": False,
            "meaning": "pair1 is not proven to be the selector-relevant target",
        },
        {
            "id": "h34_no_basis_covariance_target_independence",
            "source": "H34",
            "path": ("result",),
            "expected": "strict_core_supports_only_local_chart_embeddings_for_psi0_and_not_a_basis_covariant_or_target_independent_selector_reduction",
            "meaning": "no basis-covariant / target-independent selector reduction is exported",
        },
        {
            "id": "h35_no_strict_axis_selection",
            "source": "H35",
            "path": ("result",),
            "expected": "strict_core_supports_only_a_coordinate_level_direction_u_psi0_pair1_inside_pair1_and_not_a_strict_physical_axis_selection",
            "meaning": "no strict physical axis selection is exported on pair1",
        },
        {
            "id": "h42_bare_c_operator_is_trivial",
            "source": "H42",
            "path": ("case_without_psi0", "breaks_O2"),
            "expected": False,
            "meaning": "bare c-based retardation does not break O(2)",
        },
        {
            "id": "h42_nontrivial_split_requires_imported_anchor",
            "source": "H42",
            "path": ("case_with_psi0", "anchor_origin"),
            "expected": "imported_psi0",
            "meaning": "the nontrivial c-based split uses imported psi0",
        },
        {
            "id": "h42_psi0_plus_c_not_strict_core",
            "source": "H42",
            "path": ("hard_limits", "psi0_plus_c_is_strict_core"),
            "expected": False,
            "meaning": "psi0 plus c is not classified as strict core",
        },
        {
            "id": "p1_extension_lane_only",
            "source": "P1",
            "path": ("lane",),
            "expected": "hypothesis_extension_only",
            "meaning": "the computed pair1 split is extension-only",
        },
        {
            "id": "p1_anchor_imported_classifier",
            "source": "P1",
            "path": ("selected_operator", "classifier"),
            "expected": "ANCHOR_IMPORTED_SPLIT",
            "meaning": "the computed pair1 split is classified as anchor-imported",
        },
        {
            "id": "p1_no_strict_core_promotion",
            "source": "P1",
            "path": ("strict_core_promotion",),
            "expected": False,
            "meaning": "the computed pair1 split is not promoted to strict core",
        },
    ]

    checks: list[dict[str, Any]] = []
    mismatches: list[str] = []
    for item in checks_spec:
        try:
            actual = get_by_path(sources[item["source"]], item["path"])
        except KeyError as exc:
            actual = f"MISSING:{exc}"
        ok = actual == item["expected"]
        checks.append(
            {
                "id": item["id"],
                "source": item["source"],
                "path": ".".join(item["path"]),
                "actual": actual,
                "expected": item["expected"],
                "pass": ok,
                "meaning": item["meaning"],
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        summary = {
            "step": "N4",
            "status": "N4_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FRONTIER",
            "goal": "Check whether the current repo state already exports psi0 as a strict-core selector source without import.",
            "scope": "current_repo_state_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "current_repo_theorem_result": {
                "discharged": False,
                "reason": "the expected current-repo frontier has changed or a required evidence object is missing",
            },
            "required_next_step": "REVIEW_THE_CHANGED_FRONTIER_BEFORE_CLAIMING_CURRENT_REPO_PSI0_NONDERIVATION",
            "hard_limits": [
                "no theorem-level pass",
                "no full-closure pass",
                "no claim that a future strict-core source cannot be added",
            ],
        }
    else:
        summary = {
            "step": "N4",
            "status": "N4_DISCHARGED_CURRENT_REPO_PSI0_STRICT_CORE_NONDERIVATION_NO_FALSE_PASS",
            "goal": "Discharge a current-repo theorem: psi0 is not currently exported as a strict-core selector source, and every current computable selector split still requires extension.",
            "scope": "current_repo_state_only",
            "checks": checks,
            "current_repo_theorem_result": {
                "discharged": True,
                "current_repo_strict_core_psi0_selector_source_present": False,
                "current_repo_strict_core_map_to_nontrivial_A1_pair1_present": False,
                "current_repo_computable_selector_split_requires_extension": True,
            },
            "missing_strict_core_obligations": [
                "strict_core_selector_datum_export_for_psi0",
                "strict_core_selector_reduction_from_psi0_to_pair1_or_equivalent_target",
                "basis_covariance_or_target_independence_for_selector_reduction",
                "strict_physical_axis_selection_or_equivalent_orientation_object",
                "strict_core_operator_map_producing_nontrivial_A1_pair1_without_imported_anchor",
            ],
            "hard_limits": [
                "no theorem-level pass beyond current repo state",
                "no full-closure pass",
                "no claim that future strict-core work cannot add a new source",
                "no claim that QW-2191 is discharged",
                "no claim that ToE is closed",
            ],
            "required_next_step": "ADD_ONE_REAL_STRICT_CORE_ANCHOR_SOURCE_OBJECT_OR_PROVE_A_STRONGER_IMPOSSIBILITY_THEOREM_BEYOND_CURRENT_REPO_SCOPE",
        }

    out = ROOT / "generated" / "n4_current_repo_psi0_strict_core_nonderivation_theorem_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
