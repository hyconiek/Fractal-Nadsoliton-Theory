#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

AS_OF = "2026-03-16"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_A1 = GENERATED / "a_1_pair1_orientation_projector_operator_strict_core_v1.json"
IN_A2 = GENERATED / "a_2_pair2_orientation_projector_operator_strict_core_v1.json"
IN_A3 = GENERATED / "a_3_pair3_orientation_projector_operator_strict_core_v1.json"
IN_A4 = GENERATED / "a_4_pair4_orientation_projector_operator_strict_core_v1.json"
IN_A5 = GENERATED / "a_5_pair5_orientation_projector_operator_strict_core_v1.json"

IN_O12 = GENERATED / "o12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1.json"

IN_O13 = GENERATED / "o13_pair1_pair3_selector_chart_transport_operator_oriented_alpha13_mod_2pi_strict_convention_v1.json"
IN_O14 = GENERATED / "o14_pair1_pair4_selector_chart_transport_operator_oriented_alpha14_mod_2pi_strict_convention_v1.json"
IN_O15 = GENERATED / "o15_pair1_pair5_selector_chart_transport_operator_oriented_alpha15_mod_2pi_strict_convention_v1.json"
IN_O23 = GENERATED / "o23_pair2_pair3_selector_chart_transport_operator_oriented_alpha23_mod_2pi_strict_convention_v1.json"
IN_O24 = GENERATED / "o24_pair2_pair4_selector_chart_transport_operator_oriented_alpha24_mod_2pi_strict_convention_v1.json"
IN_O25 = GENERATED / "o25_pair2_pair5_selector_chart_transport_operator_oriented_alpha25_mod_2pi_strict_convention_v1.json"
IN_O34 = GENERATED / "o34_pair3_pair4_selector_chart_transport_operator_oriented_alpha34_mod_2pi_strict_convention_v1.json"
IN_O35 = GENERATED / "o35_pair3_pair5_selector_chart_transport_operator_oriented_alpha35_mod_2pi_strict_convention_v1.json"
IN_O45 = GENERATED / "o45_pair4_pair5_selector_chart_transport_operator_oriented_alpha45_mod_2pi_strict_convention_v1.json"

OUT_JSON = GENERATED / "p471_current_strict_pair12345_oriented_transport_operator_level_cocycle_audit_probe.json"
OUT_SUMMARY = GENERATED / "p471_current_strict_pair12345_oriented_transport_operator_level_cocycle_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def max_abs(a: np.ndarray) -> float:
    return float(np.max(np.abs(a)))


def as_vector(x: Any, n: int, label: str) -> np.ndarray:
    if not (isinstance(x, list) and len(x) == n and all(isinstance(v, (int, float)) for v in x)):
        raise ValueError(f"{label} must be a length-{n} numeric list")
    v = np.array([float(v) for v in x], dtype=float)
    if not np.isfinite(v).all():
        raise ValueError(f"{label} must contain only finite numbers")
    return v


def as_square_matrix(x: Any, n: int, label: str) -> np.ndarray:
    if not (isinstance(x, list) and len(x) == n and all(isinstance(row, list) and len(row) == n for row in x)):
        raise ValueError(f"{label} must be a {n}x{n} nested list")
    a = np.array(x, dtype=float)
    if not np.isfinite(a).all():
        raise ValueError(f"{label} must contain only finite numbers")
    return a


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    required = (
        IN_A1,
        IN_A2,
        IN_A3,
        IN_A4,
        IN_A5,
        IN_O12,
        IN_O13,
        IN_O14,
        IN_O15,
        IN_O23,
        IN_O24,
        IN_O25,
        IN_O34,
        IN_O35,
        IN_O45,
    )
    missing = [str(p.relative_to(REPO)) for p in required if not p.exists()]
    if missing:
        summary = {
            "stage": "P471",
            "status": "NOT_COMPUTABLE_MISSING_INPUTS",
            "missing": missing,
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    a1_obj = load_json(IN_A1)
    a2_obj = load_json(IN_A2)
    a3_obj = load_json(IN_A3)
    a4_obj = load_json(IN_A4)
    a5_obj = load_json(IN_A5)

    o12_obj = load_json(IN_O12)
    o13_obj = load_json(IN_O13)
    o14_obj = load_json(IN_O14)
    o15_obj = load_json(IN_O15)
    o23_obj = load_json(IN_O23)
    o24_obj = load_json(IN_O24)
    o25_obj = load_json(IN_O25)
    o34_obj = load_json(IN_O34)
    o35_obj = load_json(IN_O35)
    o45_obj = load_json(IN_O45)

    try:
        u1 = as_vector((a1_obj.get("data") or {}).get("u_1"), n=12, label="A1.u_1")
        u2 = as_vector((a2_obj.get("data") or {}).get("u_2"), n=12, label="A2.u_2")
        u3 = as_vector((a3_obj.get("data") or {}).get("u_3"), n=12, label="A3.u_3")
        u4 = as_vector((a4_obj.get("data") or {}).get("u_4"), n=12, label="A4.u_4")
        u5 = as_vector((a5_obj.get("data") or {}).get("u_5"), n=12, label="A5.u_5")
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_INPUT_VECTORS", "error": str(exc)}, ensure_ascii=True))

    try:
        O12 = as_square_matrix((o12_obj.get("outputs") or {}).get("O12"), n=12, label="O12.outputs.O12")
        O13 = as_square_matrix((o13_obj.get("outputs") or {}).get("O"), n=12, label="O13.outputs.O")
        O14 = as_square_matrix((o14_obj.get("outputs") or {}).get("O"), n=12, label="O14.outputs.O")
        O15 = as_square_matrix((o15_obj.get("outputs") or {}).get("O"), n=12, label="O15.outputs.O")
        O23 = as_square_matrix((o23_obj.get("outputs") or {}).get("O"), n=12, label="O23.outputs.O")
        O24 = as_square_matrix((o24_obj.get("outputs") or {}).get("O"), n=12, label="O24.outputs.O")
        O25 = as_square_matrix((o25_obj.get("outputs") or {}).get("O"), n=12, label="O25.outputs.O")
        O34 = as_square_matrix((o34_obj.get("outputs") or {}).get("O"), n=12, label="O34.outputs.O")
        O35 = as_square_matrix((o35_obj.get("outputs") or {}).get("O"), n=12, label="O35.outputs.O")
        O45 = as_square_matrix((o45_obj.get("outputs") or {}).get("O"), n=12, label="O45.outputs.O")
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_TRANSPORT_OPERATOR_SHAPE", "error": str(exc)}, ensure_ascii=True))

    def op_res(path: np.ndarray, direct: np.ndarray) -> float:
        return max_abs(path - direct)

    def v_res(path: np.ndarray, direct: np.ndarray, u: np.ndarray) -> float:
        return float(np.linalg.norm((path @ u) - (direct @ u)))

    # Triple list matches the previously exported section-level cocycle audits (P470/F467).
    triples = [
        ("pair1-pair2-pair3", O23 @ O12, O13, u1, "O_23 O_12", "O_13", "u_1"),
        ("pair1-pair2-pair4", O24 @ O12, O14, u1, "O_24 O_12", "O_14", "u_1"),
        ("pair1-pair2-pair5", O25 @ O12, O15, u1, "O_25 O_12", "O_15", "u_1"),
        ("pair1-pair3-pair4", O34 @ O13, O14, u1, "O_34 O_13", "O_14", "u_1"),
        ("pair1-pair3-pair5", O35 @ O13, O15, u1, "O_35 O_13", "O_15", "u_1"),
        ("pair1-pair4-pair5", O45 @ O14, O15, u1, "O_45 O_14", "O_15", "u_1"),
        ("pair2-pair3-pair4", O34 @ O23, O24, u2, "O_34 O_23", "O_24", "u_2"),
        ("pair2-pair3-pair5", O35 @ O23, O25, u2, "O_35 O_23", "O_25", "u_2"),
        ("pair2-pair4-pair5", O45 @ O24, O25, u2, "O_45 O_24", "O_25", "u_2"),
        ("pair3-pair4-pair5", O45 @ O34, O35, u3, "O_45 O_34", "O_35", "u_3"),
    ]

    audits = []
    for triple, path, direct, u_base, path_label, direct_label, u_label in triples:
        audits.append(
            {
                "triple": triple,
                "path": path_label,
                "direct": direct_label,
                "u_base": u_label,
                "operator_level_max_abs_residual": op_res(path, direct),
                "vector_section_level_l2_residual_on_u_base": v_res(path, direct, u_base),
            }
        )

    op_max = max(float(item["operator_level_max_abs_residual"]) for item in audits)
    v_max = max(float(item["vector_section_level_l2_residual_on_u_base"]) for item in audits)

    vector_section_cocycle_tolerance = 1e-12
    operator_level_equality_tolerance = 1e-12

    vector_section_cocycle_holds = bool(v_max <= vector_section_cocycle_tolerance)
    operator_level_cocycle_holds = bool(op_max <= operator_level_equality_tolerance)

    artifact = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P471",
        "goal": "audit_operator_level_vs_vector_section_level_cocycle_for_pair12345_oriented_transport",
        "as_of": AS_OF,
        "inputs": {
            "o12": str(IN_O12.relative_to(REPO)),
            "o13": str(IN_O13.relative_to(REPO)),
            "o14": str(IN_O14.relative_to(REPO)),
            "o15": str(IN_O15.relative_to(REPO)),
            "o23": str(IN_O23.relative_to(REPO)),
            "o24": str(IN_O24.relative_to(REPO)),
            "o25": str(IN_O25.relative_to(REPO)),
            "o34": str(IN_O34.relative_to(REPO)),
            "o35": str(IN_O35.relative_to(REPO)),
            "o45": str(IN_O45.relative_to(REPO)),
            "a1_pair1": str(IN_A1.relative_to(REPO)),
            "a2_pair2": str(IN_A2.relative_to(REPO)),
            "a3_pair3": str(IN_A3.relative_to(REPO)),
            "a4_pair4": str(IN_A4.relative_to(REPO)),
            "a5_pair5": str(IN_A5.relative_to(REPO)),
        },
        "audits": {
            "triple_cocycle_operator_vs_section_level": audits,
            "tolerances": {
                "vector_section_level_l2": vector_section_cocycle_tolerance,
                "operator_level_max_abs": operator_level_equality_tolerance,
            },
        },
        "result": {
            "vector_section_level_cocycle_holds_on_exported_u_section": vector_section_cocycle_holds,
            "operator_level_cocycle_holds_as_matrix_equality": operator_level_cocycle_holds,
            "operator_level_cocycle_expected": False,
            "meaning": (
                "The exported oriented transport lift is a chart-gluing ingredient for the exported ray/vector section; "
                "it is not a global transition groupoid on the full carrier. Operator-level cocycle equalities O_jk O_ij = O_ik "
                "must not be assumed for arbitrary vectors."
            ),
        },
        "hard_limits": [
            "Probe-level numeric audit only; does not export a global selector atlas nor global transition object on C_v1.",
            "Does not discharge QW-2191 nor imply strict-core selector closure.",
            "Does not promote the oriented lift to a sign-sensitive physical orientation datum.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P471",
        "status": "PASS_PROBE_READY",
        "generated_utc": artifact["generated_utc"],
        "audits": {
            "vector_section_level_cocycle_max_l2_residual": v_max,
            "operator_level_cocycle_max_abs_residual": op_max,
        },
        "result": artifact["result"],
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

