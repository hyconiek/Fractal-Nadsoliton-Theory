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
        "t164_fixing_datum": out
        / "kappa_z12_generator_orientation_canonical_fixing_datum_strict_provenance_v1.json",
        "s_dir_pair1": out / "s_dir_pair1_strict_v1.json",
        "global_directed_state": out / "selector_state_global_c_v1_directed_strict_v1.json",
        "f474_summary": out
        / "f474_current_strict_t171_global_directed_selector_state_datum_export_packet_summary.json",
    }

    missing_required = [
        key for key, path in required_paths.items() if not path.is_file()
    ]

    checked_inputs: dict[str, dict] = {}
    evidence: dict[str, object] = {}

    status = "PASS_PARTIAL_NO_STRICT_DIRECTED_AXIS_SELECTION"
    frontier = "H36_B1"
    result = "strict_core_supports_only_a_coordinate_level_undirected_axis_representative_u_psi0_pair1_inside_pair1_and_not_a_strict_directed_orientation_selection"
    frontier_text = (
        "strict core supports only a coordinate-level undirected axis representative u_psi0_pair1 inside pair1 and contains no strict argument selecting a directed orientation on that axis"
    )

    if missing_required:
        status = "FAIL_MISSING_REQUIRED_INPUTS"
    else:
        t164 = _read_json(required_paths["t164_fixing_datum"])
        s_dir = _read_json(required_paths["s_dir_pair1"])
        directed_state = _read_json(required_paths["global_directed_state"])
        f474 = _read_json(required_paths["f474_summary"])

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
        checked_inputs["global_directed_state"] = {
            "path": _p(required_paths["global_directed_state"]),
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

        evidence["t164_object"] = t164.get("object")
        evidence["s_dir_object"] = s_dir.get("object")
        evidence["directed_state_object"] = directed_state.get("object")
        evidence["t171_discharged"] = f474.get("t171_discharged")

        ok = True
        ok = ok and (t164.get("object") == "Kappa_Z12_generator_orientation_canonical_fixing_datum_strict_provenance_v1")
        ok = ok and (t164.get("no_false_pass") is True)
        ok = ok and (s_dir.get("object") == "S_dir_pair1_strict_v1")
        ok = ok and (s_dir.get("no_false_pass") is True)
        ok = ok and (
            isinstance(s_dir.get("definition"), dict)
            and bool(s_dir["definition"].get("odd_under_sign"))
        )
        ok = ok and (
            directed_state.get("object")
            == "SelectorState_global_C_v1_directed_strict_v1"
        )
        ok = ok and (directed_state.get("no_false_pass") is True)
        ok = ok and (
            isinstance(directed_state.get("state_type"), dict)
            and directed_state["state_type"].get("level") == "directed_vector_state"
        )
        ok = ok and (f474.get("no_false_pass") is True)
        ok = ok and (f474.get("t171_discharged") is True)

        if ok:
            status = "PASS_DIRECTED_AXIS_ORIENTATION_PRESENT__PREMISE_BASED_T164"
            frontier = "H36_B2"
            result = "strict_core_exports_a_directed_sign_sensitive_orientation_selection_mechanism_on_pair1_in_the_declared_premise_based_scope_T164_T171__psi0_lane_remains_axis_only"
            frontier_text = (
                "strict core exports a directed/sign-sensitive orientation selection mechanism on pair1 in the declared premise-based scope (T164+T171), while psi0 alone still supplies only an undirected axis representative"
            )
        else:
            status = "FAIL_INPUT_CONTRADICTION_OR_UNEXPECTED_OBJECTS"

    data = {
        "id": "H36",
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
            "no claim that psi0 selects a physically directed axis",
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

    (out / "h36_directed_axis_orientation_audit.json").write_text(
        json.dumps(data, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    (out / "h36_directed_axis_orientation_audit_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )


if __name__ == "__main__":
    main()
