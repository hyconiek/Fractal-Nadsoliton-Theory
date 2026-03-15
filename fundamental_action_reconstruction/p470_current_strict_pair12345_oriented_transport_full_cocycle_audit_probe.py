#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

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

OUT_JSON = GENERATED / "p470_current_strict_pair12345_oriented_transport_full_cocycle_audit_probe.json"
OUT_SUMMARY = GENERATED / "p470_current_strict_pair12345_oriented_transport_full_cocycle_audit_probe_summary.json"


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
            "stage": "P470",
            "status": "NOT_COMPUTABLE_MISSING_INPUTS",
            "missing": missing,
            "theorem_level_pass": False,
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

    n = 12
    I = np.eye(n, dtype=float)

    # Orthogonality / involution checks.
    orth = {k: max_abs(O.T @ O - I) for k, O in {"o12": O12, "o13": O13, "o14": O14, "o15": O15, "o23": O23, "o24": O24, "o25": O25, "o34": O34, "o35": O35, "o45": O45}.items()}
    invol = {k: max_abs(O @ O - I) for k, O in {"o12": O12, "o13": O13, "o14": O14, "o15": O15, "o23": O23, "o24": O24, "o25": O25, "o34": O34, "o35": O35, "o45": O45}.items()}

    # Vector-level transport checks.
    def v_res(O: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
        return float(np.linalg.norm((O @ a) - b))

    transport = {
        "12": {"u2_minus_O12_u1_l2": v_res(O12, u1, u2), "u1_minus_O12_u2_l2": v_res(O12, u2, u1)},
        "13": {"u3_minus_O13_u1_l2": v_res(O13, u1, u3), "u1_minus_O13_u3_l2": v_res(O13, u3, u1)},
        "14": {"u4_minus_O14_u1_l2": v_res(O14, u1, u4), "u1_minus_O14_u4_l2": v_res(O14, u4, u1)},
        "15": {"u5_minus_O15_u1_l2": v_res(O15, u1, u5), "u1_minus_O15_u5_l2": v_res(O15, u5, u1)},
        "23": {"u3_minus_O23_u2_l2": v_res(O23, u2, u3), "u2_minus_O23_u3_l2": v_res(O23, u3, u2)},
        "24": {"u4_minus_O24_u2_l2": v_res(O24, u2, u4), "u2_minus_O24_u4_l2": v_res(O24, u4, u2)},
        "25": {"u5_minus_O25_u2_l2": v_res(O25, u2, u5), "u2_minus_O25_u5_l2": v_res(O25, u5, u2)},
        "34": {"u4_minus_O34_u3_l2": v_res(O34, u3, u4), "u3_minus_O34_u4_l2": v_res(O34, u4, u3)},
        "35": {"u5_minus_O35_u3_l2": v_res(O35, u3, u5), "u3_minus_O35_u5_l2": v_res(O35, u5, u3)},
        "45": {"u5_minus_O45_u4_l2": v_res(O45, u4, u5), "u4_minus_O45_u5_l2": v_res(O45, u5, u4)},
    }

    # Full triple cocycle/path-independence checks on oriented vectors.
    def cocycle_residual(*, path: np.ndarray, direct: np.ndarray, u: np.ndarray) -> float:
        return float(np.linalg.norm((path @ u) - (direct @ u)))

    cocycle = {
        "pair1_pair2_pair3": cocycle_residual(path=(O23 @ O12), direct=O13, u=u1),
        "pair1_pair2_pair4": cocycle_residual(path=(O24 @ O12), direct=O14, u=u1),
        "pair1_pair2_pair5": cocycle_residual(path=(O25 @ O12), direct=O15, u=u1),
        "pair1_pair3_pair4": cocycle_residual(path=(O34 @ O13), direct=O14, u=u1),
        "pair1_pair3_pair5": cocycle_residual(path=(O35 @ O13), direct=O15, u=u1),
        "pair1_pair4_pair5": cocycle_residual(path=(O45 @ O14), direct=O15, u=u1),
        "pair2_pair3_pair4": cocycle_residual(path=(O34 @ O23), direct=O24, u=u2),
        "pair2_pair3_pair5": cocycle_residual(path=(O35 @ O23), direct=O25, u=u2),
        "pair2_pair4_pair5": cocycle_residual(path=(O45 @ O24), direct=O25, u=u2),
        "pair3_pair4_pair5": cocycle_residual(path=(O45 @ O34), direct=O35, u=u3),
        "meaning": "full triple-overlap path-independence on the exported oriented vector section (tracked sign convention)",
    }

    artifact = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P470",
        "goal": "audit_pair12345_oriented_transport_mod_2pi_vector_transport_and_full_triple_cocycle_path_independence",
        "inputs": {k: str(p.relative_to(REPO)) for k, p in {
            "a1_pair1": IN_A1,
            "a2_pair2": IN_A2,
            "a3_pair3": IN_A3,
            "a4_pair4": IN_A4,
            "a5_pair5": IN_A5,
            "o12": IN_O12,
            "o13": IN_O13,
            "o14": IN_O14,
            "o15": IN_O15,
            "o23": IN_O23,
            "o24": IN_O24,
            "o25": IN_O25,
            "o34": IN_O34,
            "o35": IN_O35,
            "o45": IN_O45,
        }.items()},
        "audits": {
            "orthogonality": orth,
            "involution": invol,
            "transport_vector_level": transport,
            "cocycle_vector_level": cocycle,
        },
        "hard_limits": [
            "Probe-level numeric audit only; does not export a global selector atlas nor global cocycle data on C_v1.",
            "Oriented transport is a tracked gauge/convention lift, not a sign-sensitive physical orientation datum.",
            "Does not discharge QW-2191 nor imply strict-core selector closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    cocycle_max = max(float(v) for k, v in cocycle.items() if k.startswith("pair"))
    transport_max = max(float(v) for edge in transport.values() for v in edge.values())

    summary = {
        "stage": "P470",
        "status": "PASS_PROBE_READY",
        "generated_utc": artifact["generated_utc"],
        "audits": {
            "orthogonality_max_abs_residual": max(float(v) for v in orth.values()),
            "involution_max_abs_residual": max(float(v) for v in invol.values()),
            "transport_vector_level_max_l2_residual": transport_max,
            "cocycle_vector_level_max_l2_residual": cocycle_max,
        },
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

