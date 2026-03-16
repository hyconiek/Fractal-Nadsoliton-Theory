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
        "projective_state": out / "selector_state_global_c_v1_projective_strict_v1.json",
        "directed_state": out / "selector_state_global_c_v1_directed_strict_v1.json",
        "p632_directed_decision_summary": out
        / "p632_current_strict_directed_continuation_decision_packet_summary.json",
    }

    missing_required = [
        key for key, path in required_paths.items() if not path.is_file()
    ]

    checked_inputs: dict[str, dict] = {}
    evidence: dict[str, object] = {}

    status = (
        "PASS_PARTIAL_PROJECTIVE_STATE_ONLY__GLOBAL_PROJECTIVE_STATE_OBJECT_EXPORTED_BUT_NO_SIGN_SENSITIVE_ORIENTATION_DATUM"
    )
    frontier = "H38_B1"
    result = "strict_core_supports_at_most_a_local_projective_or_ray_level_selector_representative_on_pair1_and_does_not_furnish_a_physically_individuated_directed_selector_state"
    frontier_text = "strict core supports at most a local projective or ray-level selector representative on pair1 and does not furnish a physically individuated directed selector state"

    if missing_required:
        status = "FAIL_MISSING_REQUIRED_INPUTS"
    else:
        projective_state = _read_json(required_paths["projective_state"])
        directed_state = _read_json(required_paths["directed_state"])
        p632 = _read_json(required_paths["p632_directed_decision_summary"])

        sign_gauge_projective = (
            projective_state.get("state_type", {}).get("sign_gauge")
            if isinstance(projective_state.get("state_type"), dict)
            else None
        )
        sign_gauge_directed = (
            directed_state.get("state_type", {}).get("sign_gauge")
            if isinstance(directed_state.get("state_type"), dict)
            else None
        )

        checked_inputs["projective_state"] = {
            "path": _p(required_paths["projective_state"]),
            "object": projective_state.get("object"),
            "no_false_pass": projective_state.get("no_false_pass"),
            "sign_gauge": sign_gauge_projective,
        }
        checked_inputs["directed_state"] = {
            "path": _p(required_paths["directed_state"]),
            "object": directed_state.get("object"),
            "no_false_pass": directed_state.get("no_false_pass"),
            "sign_gauge": sign_gauge_directed,
        }
        checked_inputs["p632_directed_decision_summary"] = {
            "path": _p(required_paths["p632_directed_decision_summary"]),
            "no_false_pass": p632.get("no_false_pass"),
            "selected_continuation": p632.get("selected_continuation"),
            "decision": p632.get("decision"),
        }

        evidence["projective_sign_gauge"] = sign_gauge_projective
        evidence["directed_sign_gauge"] = sign_gauge_directed
        evidence["directed_continuation_selected"] = (
            p632.get("selected_continuation") == "directed"
        )

        ok = True
        ok = ok and (projective_state.get("no_false_pass") is True)
        ok = ok and bool(sign_gauge_projective)
        ok = ok and (directed_state.get("no_false_pass") is True)
        ok = ok and (
            directed_state.get("object") == "SelectorState_global_C_v1_directed_strict_v1"
        )
        ok = ok and bool(sign_gauge_directed)
        ok = ok and (p632.get("no_false_pass") is True)
        ok = ok and (p632.get("selected_continuation") == "directed")

        if ok:
            status = "PASS_DIRECTED_CONTINUATION_SELECTED__GLOBAL_DIRECTED_STATE_OBJECT_EXPORTED__PROJECTIVE_STATE_RETAINED_AS_QUOTIENT"
            frontier = "H38_B2"
            result = "strict_core_exports_both_a_global_projective_selector_state_object_and_a_global_directed_selector_state_object__directed_state_descends_to_projective_state__premise_based_via_T164"
            frontier_text = "strict core exports both a global projective selector state object (quotient/ray layer) and a global directed selector state object (vector layer, premise-based via T164), with the directed layer descending to the projective one"
        else:
            status = "FAIL_INPUT_CONTRADICTION_OR_UNEXPECTED_OBJECTS"

    data = {
        "id": "H38",
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

    (out / "h38_projective_selector_state_audit.json").write_text(
        json.dumps(data, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    (out / "h38_projective_selector_state_audit_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )


if __name__ == "__main__":
    main()
