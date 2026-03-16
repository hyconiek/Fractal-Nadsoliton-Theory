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
        "t164_fixing_datum": out
        / "kappa_z12_generator_orientation_canonical_fixing_datum_strict_provenance_v1.json",
        "s_dir_pair1": out / "s_dir_pair1_strict_v1.json",
        "selector_state_global_directed": out
        / "selector_state_global_c_v1_directed_strict_v1.json",
        "f474_summary": out
        / "f474_current_strict_t171_global_directed_selector_state_datum_export_packet_summary.json",
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
    status = (
        "PASS_PARTIAL_SIGN_TRACKED_CONVENTION_LAYER_PRESENT_BUT_NO_STRICT_PHYSICAL_SIGN_DISTINCTION_OBSERVABLE"
    )
    result = "strict_core_contains_no_sign_sensitive_state_object_or_observable_on_pair1_and_therefore_does_not_distinguish_u_from_minus_u_as_physically_different_selector_states"
    frontier = "H37_B2"
    frontier_text = (
        "strict core exports sign-tracked convention-layer oriented vectors and a global projective selector state object, but still contains no strict sign-sensitive physical state object or observable distinguishing u from -u on pair1"
    )

    if missing_required:
        status = "FAIL_MISSING_REQUIRED_INPUTS"
    else:
        projective_state = _read_json(
            required_paths["selector_state_global_projective"]
        )
        t164 = _read_json(required_paths["t164_fixing_datum"])
        s_dir = _read_json(required_paths["s_dir_pair1"])
        directed_state = _read_json(required_paths["selector_state_global_directed"])
        f474 = _read_json(required_paths["f474_summary"])
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
        checked_inputs["t164_fixing_datum"] = {
            "path": _p(required_paths["t164_fixing_datum"]),
            "object": t164.get("object"),
            "no_false_pass": t164.get("no_false_pass"),
        }
        checked_inputs["s_dir_pair1"] = {
            "path": _p(required_paths["s_dir_pair1"]),
            "object": s_dir.get("object"),
            "no_false_pass": s_dir.get("no_false_pass"),
            "odd_under_sign": (s_dir.get("definition") or {}).get("odd_under_sign")
            if isinstance(s_dir.get("definition"), dict)
            else None,
        }
        checked_inputs["selector_state_global_directed"] = {
            "path": _p(required_paths["selector_state_global_directed"]),
            "object": directed_state.get("object"),
            "no_false_pass": directed_state.get("no_false_pass"),
            "state_level": (directed_state.get("state_type") or {}).get("level")
            if isinstance(directed_state.get("state_type"), dict)
            else None,
        }
        checked_inputs["f474_summary"] = {
            "path": _p(required_paths["f474_summary"]),
            "no_false_pass": f474.get("no_false_pass"),
            "t171_discharged": f474.get("t171_discharged"),
            "h37_discharged": f474.get("h37_discharged"),
        }
        checked_inputs["p472_reflection_breaking_weight_scan_summary"] = {
            "path": _p(required_paths["p472_reflection_breaking_weight_scan_summary"]),
            "no_false_pass": p472_summary.get("no_false_pass"),
            "outside_k_total_rows": outside_k_total,
        }

        evidence["strict_global_state_sign_gauge"] = sign_gauge
        evidence["t164_object"] = t164.get("object")
        evidence["s_dir_object"] = s_dir.get("object")
        evidence["directed_state_object"] = directed_state.get("object")
        evidence["t171_discharged"] = f474.get("t171_discharged")
        evidence["h37_discharged"] = f474.get("h37_discharged")
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

        # Directed sign distinction is present if T164 + S_dir + directed global state are exported,
        # and F474 explicitly marks the discharge.
        directed_ok = True
        directed_ok = directed_ok and (projective_state.get("no_false_pass") is True)
        directed_ok = directed_ok and bool(sign_gauge)
        directed_ok = directed_ok and (
            t164.get("object")
            == "Kappa_Z12_generator_orientation_canonical_fixing_datum_strict_provenance_v1"
        )
        directed_ok = directed_ok and (t164.get("no_false_pass") is True)
        directed_ok = directed_ok and (s_dir.get("object") == "S_dir_pair1_strict_v1")
        directed_ok = directed_ok and (s_dir.get("no_false_pass") is True)
        directed_ok = directed_ok and (
            isinstance(s_dir.get("definition"), dict)
            and bool(s_dir["definition"].get("odd_under_sign"))
        )
        directed_ok = directed_ok and (
            directed_state.get("object")
            == "SelectorState_global_C_v1_directed_strict_v1"
        )
        directed_ok = directed_ok and (directed_state.get("no_false_pass") is True)
        directed_ok = directed_ok and (
            isinstance(directed_state.get("state_type"), dict)
            and directed_state["state_type"].get("level") == "directed_vector_state"
        )
        directed_ok = directed_ok and (f474.get("no_false_pass") is True)
        directed_ok = directed_ok and (f474.get("t171_discharged") is True)
        directed_ok = directed_ok and (f474.get("h37_discharged") is True)
        directed_ok = directed_ok and (p472_summary.get("no_false_pass") is True)
        directed_ok = directed_ok and (
            not (isinstance(outside_k_total, int) and outside_k_total > 0)
        )

        if directed_ok:
            status = "PASS_STRICT_DIRECTED_SIGN_SENSITIVE_DISTINCTION_OBSERVABLE_EXPORTED_AND_GLOBAL_DIRECTED_SELECTOR_STATE_OBJECT_EXPORTED__PREMISE_BASED_T164"
            result = "strict_core_exports_an_explicit_sign_sensitive_directed_observable_on_pair1_and_a_global_directed_selector_state_object_on_C_v1__premise_based_via_T164__u_and_minus_u_are_distinguishable_in_declared_scope"
            frontier = "H37_B3"
            frontier_text = (
                "strict core exports an explicit sign-sensitive directed observable on pair1 and a global directed selector state object on C_v1 (premise-based via T164), descending to the projective state; Aut(Z_12)-invariant canonicity remains undischarged and QW-2191 remains open"
            )

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
        "result": result,
        "frontier": frontier,
        "frontier_text": frontier_text,
        "hard_limits": [
            "no theorem-level pass",
            "no full-closure pass",
            "no claim of Aut(Z_12)-invariant canonicity (premise-based via T164)",
            "no claim that QW-2191 is discharged",
        ],
        "no_false_pass": True,
    }

    summary = {
        "id": data["id"],
        "status": data["status"],
        "frontier": data["frontier"],
        "result": data["result"],
        "required_inputs_missing": missing_required,
        "no_false_pass": True,
    }

    (out / "h37_sign_distinction_state_audit.json").write_text(
        json.dumps(data, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    (out / "h37_sign_distinction_state_audit_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )


if __name__ == "__main__":
    main()
