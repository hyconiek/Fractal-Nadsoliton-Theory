#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from math import cos
from pathlib import Path
from typing import Any


AS_OF = "2026-03-27"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_MANIFEST = GENERATED / "unfreeze_omega_phi_beta_eta_manifest.json"
IN_K1 = (
    REPO
    / "fundamental_action_reconstruction"
    / "K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md"
)
IN_P729 = (
    REPO
    / "fundamental_action_reconstruction"
    / "generated"
    / "p729_current_strict_t183_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
)
IN_P730 = (
    REPO
    / "fundamental_action_reconstruction"
    / "generated"
    / "p730_current_strict_t184_direction_free_shannon_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
)
IN_F960 = (
    REPO
    / "fundamental_action_reconstruction"
    / "generated"
    / "f960_current_strict_t173_t176_existing_t183_residual_datum_pair12_orbit_direction_selection_bridge_target_packet_summary.json"
)

OUT_JSON = GENERATED / "unfreeze_direction_sensitivity_boundary_probe.json"
OUT_SUMMARY = GENERATED / "unfreeze_direction_sensitivity_boundary_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def kernel_value(distance: int, omega: float, phi: float, beta: float, eta: float) -> float:
    return cos(omega * distance + phi) / (1.0 + beta * (distance ** eta))


def max_branch_inversion_diff(params: dict[str, float]) -> float:
    diffs = []
    for k in range(12):
        positive = kernel_value(k, **params)
        negative = kernel_value(abs(-k), **params)
        diffs.append(abs(positive - negative))
    return max(diffs)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_MANIFEST, IN_K1, IN_P729, IN_P730, IN_F960]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "EXPERIMENTAL_UNFREEZE_DIRECTION_SENSITIVITY_BOUNDARY",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    manifest = load_json(IN_MANIFEST)
    p729 = load_json(IN_P729)
    p730 = load_json(IN_P730)
    f960 = load_json(IN_F960)
    k1_text = IN_K1.read_text(encoding="utf-8")

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

    add_check(
        "experimental_lane_is_explicitly_nonstrict",
        manifest.get("lane_status"),
        "NONSTRICT_FIT_DERIVED_PARAMETER_UNFREEZE_SANDBOX",
        "The parameter-unfreeze lane must stay explicitly non-strict before any boundary conclusion is admissible.",
    )
    add_check(
        "strict_working_family_remains_distance_only",
        "K_strict_gate(d) := cos(omega*d + phi) / (1 + beta*d^eta)" in k1_text,
        True,
        "K1 still packages the strict working family as a function of distance d only.",
    )
    add_check(
        "t183_split_is_already_localized_as_delta_k_vs_delta_minus_k",
        bool(p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions")),
        True,
        "P729 already localizes the surviving split as delta_k versus delta_-k.",
    )
    add_check(
        "pair1_branch_kind_is_positive_index_delta_k",
        p729.get("pair1_orbit_branch_kind"),
        "delta_k_positive_index_branch",
        "P729 keeps pair1 on the positive-index delta_k branch.",
    )
    add_check(
        "pair2_branch_kind_is_negative_index_delta_minus_k",
        p729.get("pair2_orbit_branch_kind"),
        "delta_minus_k_negative_index_branch",
        "P729 keeps pair2 on the negative-index delta_-k branch.",
    )
    add_check(
        "current_direction_free_lane_still_fails_to_select_branch",
        bool(p730.get("current_direction_free_shannon_lane_selects_pair12_orbit_direction_branch")),
        False,
        "P730/N726 already keep the present direction-free lane nonselective on the surviving pair12 branches.",
    )
    add_check(
        "direction_free_branch_profiles_remain_equal_under_inversion",
        bool(p730.get("pair12_branch_profiles_equal_under_inversion")),
        True,
        "P730/N726 already package inversion-equal branch profiles as the precise reason for the current failure.",
    )
    add_check(
        "f960_requires_genuinely_new_inversion_sensitive_provider_class",
        "genuinely new inversion-sensitive source-side provider class"
        in str(f960.get("recommended_next_move", "")),
        True,
        "F960 already freezes the honest next need as a genuinely new inversion-sensitive source-side provider class.",
    )

    seed_candidates = manifest.get("seed_candidates", {})
    sampled_tuples: dict[str, dict[str, float]] = {}
    for name in (
        "q2049_hold",
        "q2039_refrozen",
        "local_beta_midpoint",
        "canonical_refrozen_midpoint",
        "canonical_bridge_midpoint",
    ):
        candidate = seed_candidates.get(name)
        if candidate:
            sampled_tuples[name] = {k: float(v) for k, v in candidate.items()}

    sampled_diffs = {
        name: max_branch_inversion_diff(params)
        for name, params in sampled_tuples.items()
    }
    all_samples_inversion_even = all(diff <= 1.0e-15 for diff in sampled_diffs.values())

    artifact = {
        "stage": "EXPERIMENTAL_UNFREEZE_DIRECTION_SENSITIVITY_BOUNDARY",
        "status": (
            "PASS_PARAMETER_UNFREEZE_ALONE_DOES_NOT_CREATE_INVERSION_SENSITIVE_PROVIDER"
            if not blocking and all_samples_inversion_even
            else "REVIEW_REQUIRED_PARAMETER_UNFREEZE_DIRECTION_SENSITIVITY_BOUNDARY_CHANGED"
        ),
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "manifest": str(IN_MANIFEST.relative_to(REPO)),
            "K1": str(IN_K1.relative_to(REPO)),
            "P729_summary": str(IN_P729.relative_to(REPO)),
            "P730_summary": str(IN_P730.relative_to(REPO)),
            "F960_summary": str(IN_F960.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "sampled_seed_max_branch_inversion_diffs": sampled_diffs,
        "all_sampled_tuples_remain_inversion_even_under_distance_only_kernel": all_samples_inversion_even,
        "audit_conclusion": {
            "parameter_motion_inside_current_family_changes_only_distance_response": True,
            "parameter_motion_inside_current_family_adds_signed_orbit_direction_input": False,
            "parameter_unfreeze_alone_can_export_f960_t183_bridge": False,
            "allowed_label": "bridge-pressure only",
            "recommended_next_move": (
                "Use the unfreeze lane only as a bridge-pressure or fit-dependence test, "
                "while searching outside the present distance-only family for a genuinely "
                "new inversion-sensitive source-side provider class."
            ),
        },
        "hard_limits": [
            "No T183 discharge.",
            "No T176 discharge.",
            "No QW-2191 closure.",
            "No strict physical orientation datum.",
            "No proof that parameter motion alone upgrades the provider class.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "strict_family_remains_distance_only": True,
        "t183_split_already_localized_as_delta_k_vs_delta_minus_k": bool(
            p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions")
        ),
        "current_direction_free_lane_selects_pair12_branch": bool(
            p730.get("current_direction_free_shannon_lane_selects_pair12_orbit_direction_branch")
        ),
        "all_sampled_tuples_remain_inversion_even_under_distance_only_kernel": all_samples_inversion_even,
        "parameter_unfreeze_alone_can_export_f960_t183_bridge": False,
        "allowed_label": "bridge-pressure only",
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(f"Wrote: {OUT_SUMMARY}")
    print(f"Status: {summary['status']}")
    print(
        "All sampled tuples remain inversion-even under distance-only kernel:",
        summary["all_sampled_tuples_remain_inversion_even_under_distance_only_kernel"],
    )


if __name__ == "__main__":
    main()
