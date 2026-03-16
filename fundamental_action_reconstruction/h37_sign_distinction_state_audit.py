#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    root = Path(__file__).resolve().parent
    repo_root = root.parent
    out = root / "generated"
    out.mkdir(parents=True, exist_ok=True)

    def _p(path: Path) -> str:
        try:
            return str(path.relative_to(repo_root))
        except ValueError:
            return str(path)

    required_paths = {
        "selector_state_global_projective": out
        / "selector_state_global_c_v1_projective_strict_v1.json",
        "p472_reflection_breaking_weight_scan_summary": out
        / "p472_current_strict_pair1_reflection_breaking_weight_repo_scan_probe_summary.json",
    }
    optional_paths = {
        "p473_extension_lane_projector_consistency_summary": out
        / "p473_current_strict_extension_lane_global_oriented_selector_state_projector_consistency_audit_probe_summary.json",
        "extension_lane_oriented_selector_state": out
        / "strict_extension_lane_selector_state_global_c_v1_oriented_vector_v1.json",
    }

    missing_required = [
        key for key, path in required_paths.items() if not path.is_file()
    ]

    checked_inputs: dict[str, dict] = {}
    evidence: dict[str, object] = {}
    status = "PASS_PARTIAL_SIGN_TRACKED_CONVENTION_LAYER_PRESENT_BUT_NO_STRICT_PHYSICAL_SIGN_DISTINCTION_OBSERVABLE"

    if missing_required:
        status = "FAIL_MISSING_REQUIRED_INPUTS"
    else:
        projective_state = _read_json(required_paths["selector_state_global_projective"])
        p472_summary = _read_json(
            required_paths["p472_reflection_breaking_weight_scan_summary"]
        )

        sign_gauge = (
            projective_state.get("state_type", {}).get("sign_gauge")
            if isinstance(projective_state.get("state_type"), dict)
            else None
        )

        outside_k_total = p472_summary.get(
            "candidates_weight_like_strictish_reflection_breaking_and_dot_nonzero_outside_k_total_rows"
        )

        checked_inputs["selector_state_global_projective"] = {
            "path": _p(required_paths["selector_state_global_projective"]),
            "no_false_pass": projective_state.get("no_false_pass"),
            "sign_gauge": sign_gauge,
        }
        checked_inputs["p472_reflection_breaking_weight_scan_summary"] = {
            "path": _p(required_paths["p472_reflection_breaking_weight_scan_summary"]),
            "no_false_pass": p472_summary.get("no_false_pass"),
            "outside_k_total_rows": outside_k_total,
        }

        evidence["strict_global_state_sign_gauge"] = sign_gauge
        evidence["p472_weight_like_reflection_breaking_candidates_outside_k_total_rows"] = (
            outside_k_total
        )

        if projective_state.get("no_false_pass") is not True:
            status = "FAIL_INPUT_CONTRADICTION_NO_FALSE_PASS_EXPECTED"
        if p472_summary.get("no_false_pass") is not True:
            status = "FAIL_INPUT_CONTRADICTION_NO_FALSE_PASS_EXPECTED"
        if not sign_gauge:
            status = "FAIL_MISSING_EXPECTED_SIGN_GAUGE_FIELD_IN_PROJECTIVE_STATE_OBJECT"
        if isinstance(outside_k_total, int) and outside_k_total > 0:
            status = "FAIL_P472_REPORTS_STRICTISH_WEIGHT_LIKE_REFLECTION_BREAKING_CANDIDATES_OUTSIDE_K_TOTAL_ROWS"

    for key, path in optional_paths.items():
        if not path.is_file():
            checked_inputs[key] = {"path": _p(path), "present": False}
            continue
        try:
            content = _read_json(path)
        except json.JSONDecodeError:
            checked_inputs[key] = {
                "path": _p(path),
                "present": True,
                "json_parse_error": True,
            }
            continue
        checked_inputs[key] = {
            "path": _p(path),
            "present": True,
            "no_false_pass": content.get("no_false_pass"),
            "status": content.get("status"),
            "overall_pass": content.get("overall_pass"),
        }
        if key == "p473_extension_lane_projector_consistency_summary":
            evidence["p473_extension_lane_projector_consistency_status"] = content.get(
                "status"
            )

    data = {
        "id": "H37",
        "date": "2026-03-16",
        "status": status,
        "checked_inputs": checked_inputs,
        "evidence": evidence,
        "missing_required_inputs": missing_required,
        "result": "strict_core_exports_sign_tracked_convention_layer_and_global_projective_selector_state_object_but_contains_no_strict_sign_sensitive_physical_observable_distinguishing_u_from_minus_u_on_pair1",
        "frontier": "H37_B2",
        "frontier_text": "strict core exports sign-tracked convention-layer oriented vectors and a global projective selector state object, but still contains no strict sign-sensitive physical state object or observable distinguishing u from -u on pair1",
        "hard_limits": [
            "no theorem-level pass",
            "no full-closure pass",
            "no claim that sign reversal inside pair1 is physically meaningful in strict core",
            "no claim that u and -u are physically distinct selector states",
            "no claim that QW-2191 is discharged",
        ],
    }

    summary = {
        "id": data["id"],
        "status": data["status"],
        "frontier": data["frontier"],
        "result": data["result"],
        "required_inputs_missing": missing_required,
    }

    (out / "h37_sign_distinction_state_audit.json").write_text(
        json.dumps(data, indent=2) + "\n", encoding="utf-8"
    )
    (out / "h37_sign_distinction_state_audit_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
