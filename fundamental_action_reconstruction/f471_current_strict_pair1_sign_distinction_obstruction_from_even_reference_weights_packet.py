#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-16"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_A1 = GENERATED / "a_1_pair1_orientation_projector_operator_strict_core_v1.json"
IN_R_ORD = GENERATED / "r_ord_z12_v1_reference_distribution.json"

OUT_OBJ = GENERATED / "sign_distinction_obstruction_pair1_even_reference_weights_strict_v1.json"
OUT_SUMMARY = (
    GENERATED
    / "f471_current_strict_pair1_sign_distinction_obstruction_from_even_reference_weights_packet_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing = [str(p.relative_to(REPO)) for p in (IN_A1, IN_R_ORD) if not p.exists()]
    if missing:
        raise SystemExit(
            json.dumps(
                {
                    "stage": "F471",
                    "status": "NOT_COMPUTABLE_MISSING_INPUTS",
                    "missing": missing,
                    "expected": [
                        "F456 exports A_1(pair1) with u_1",
                        "F446 exports r_ord_z12 reference distribution with ord_z12_by_x",
                    ],
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    a1 = load_json(IN_A1)
    r_ord = load_json(IN_R_ORD)

    try:
        u1 = (a1.get("data") or {}).get("u_1")
        if not (isinstance(u1, list) and len(u1) == 12 and all(isinstance(v, (int, float)) for v in u1)):
            raise ValueError("A1.data.u_1 must be a length-12 numeric list")
        u1 = [float(v) for v in u1]
        coords = (a1.get("data") or {}).get("u_1_coords_in_c1_s1")
        coords_ok = isinstance(coords, list) and len(coords) == 2 and all(isinstance(v, (int, float)) for v in coords)
        u1_coords = [float(v) for v in coords] if coords_ok else None
    except Exception as exc:
        raise SystemExit(json.dumps({"stage": "F471", "status": "INVALID_A1_SHAPE", "error": str(exc)}, ensure_ascii=True))

    try:
        ord_by_x = r_ord.get("ord_z12_by_x")
        if not (isinstance(ord_by_x, list) and len(ord_by_x) == 12 and all(isinstance(v, int) for v in ord_by_x)):
            raise ValueError("r_ord.ord_z12_by_x must be a length-12 int list")
        ord_by_x = [int(v) for v in ord_by_x]
    except Exception as exc:
        raise SystemExit(json.dumps({"stage": "F471", "status": "INVALID_R_ORD_SHAPE", "error": str(exc)}, ensure_ascii=True))

    def reflect(x: int) -> int:
        return (-x) % 12

    even_ord = all(ord_by_x[x] == ord_by_x[reflect(x)] for x in range(12))

    alpha_geo = 4.0 * math.log(2.0)
    r_weights = [math.exp(-alpha_geo * float(o)) for o in ord_by_x]
    even_r = all(abs(r_weights[x] - r_weights[reflect(x)]) <= 0.0 for x in range(12))

    odd_u1 = all(abs(u1[x] + u1[reflect(x)]) <= 1e-15 for x in range(12))

    dot_ord = sum(float(o) * v for o, v in zip(ord_by_x, u1))
    dot_r = sum(w * v for w, v in zip(r_weights, u1))

    generated_utc = datetime.now(timezone.utc).isoformat()

    obj = {
        "object": "SignDistinctionObstruction_pair1_even_reference_weights_strict_v1",
        "stage": "F471",
        "status": "actual_exported_pair1_sign_distinction_obstruction_for_even_reference_weights_on_current_exported_instance__no_false_pass",
        "as_of": AS_OF,
        "generated_utc": generated_utc,
        "intent": (
            "Record the narrow obstruction that even reference weights (ord_Z12 and r_ord derived from it) cannot distinguish "
            "u_1 from -u_1 on the current exported pair1 axis because u_1 is odd under reflection while the weights are even; "
            "therefore the weighted sum Σ_x w(x) u_1(x) cancels exactly on this instance."
        ),
        "inputs": {
            "A1_pair1_ref": str(IN_A1.relative_to(REPO)),
            "r_ord_z12_ref": str(IN_R_ORD.relative_to(REPO)),
        },
        "checks": {
            "ord_even_under_reflection": bool(even_ord),
            "r_ord_weights_even_under_reflection": bool(even_r),
            "u1_odd_under_reflection": bool(odd_u1),
            "u1_coords_in_c1_s1": u1_coords,
        },
        "computed": {
            "dot_u1_with_ord_z12_by_x": float(dot_ord),
            "dot_u1_with_r_ord_weights": float(dot_r),
        },
        "meaning": [
            "On this exported instance, any sign-sensitive scalar of the form Σ_x w(x) u_1(x) vanishes for the even weights w=ord_Z12 and w=r_ord.",
            "Therefore the ord-reference infrastructure (F446/N479/N480) does not by itself supply a sign-distinction observable on pair1 for this instance.",
        ],
        "hard_limits": [
            "Obstruction record only; does not discharge H37 and does not export a sign-sensitive physical orientation datum.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim global discharge of QW-2191.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    OUT_OBJ.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")

    summary = {
        "stage": "F471",
        "status": "PASS_EXPORTED_PAIR1_SIGN_DISTINCTION_OBSTRUCTION_OBJECT",
        "exported": {"obstruction_object": str(OUT_OBJ.relative_to(REPO))},
        "dot_u1_with_ord_z12_by_x": float(dot_ord),
        "dot_u1_with_r_ord_weights": float(dot_r),
        "no_false_pass": True,
    }
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")

    print(OUT_OBJ)


if __name__ == "__main__":
    main()

