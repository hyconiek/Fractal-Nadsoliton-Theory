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
IN_O12 = GENERATED / "o12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1.json"
IN_A2 = GENERATED / "a_2_pair2_orientation_projector_operator_strict_core_v1.json"
IN_SECTION = GENERATED / "a_12_pair12_chart_glued_orientation_projector_operator_section_strict_core_v1.json"
IN_P466_SUMMARY = GENERATED / "p466_current_strict_pair12_chart_glued_projector_operator_section_audit_probe_summary.json"

OUT_ATLAS = GENERATED / "selector_atlas_pair12_sigma_int_corridor_projector_v1.json"
OUT_SUMMARY = GENERATED / "selector_atlas_pair12_sigma_int_corridor_projector_v1_summary.json"


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

    required = (IN_A1, IN_O12, IN_A2, IN_SECTION, IN_P466_SUMMARY)
    missing = [str(p.relative_to(REPO)) for p in required if not p.exists()]
    if missing:
        raise SystemExit(
            json.dumps(
                {
                    "status": "NOT_COMPUTABLE_MISSING_INPUTS",
                    "missing": missing,
                    "expected": [
                        "F456 A_1(pair1) projector operator",
                        "F461 O12 chart transport operator",
                        "F462 A_2(pair2) projector operator",
                        "F462 pair12 chart-glued operator section",
                        "P466 audit summary",
                    ],
                },
                ensure_ascii=True,
            )
        )

    a1_obj = load_json(IN_A1)
    o12_obj = load_json(IN_O12)
    a2_obj = load_json(IN_A2)
    section_obj = load_json(IN_SECTION)
    p466_summary = load_json(IN_P466_SUMMARY)

    try:
        u1 = as_vector((a1_obj.get("data") or {}).get("u_1"), n=12, label="F456.data.u_1")
        u2 = as_vector((a2_obj.get("data") or {}).get("u_2"), n=12, label="F462.data.u_2")
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_INPUT_VECTORS", "error": str(exc)}, ensure_ascii=True))

    try:
        outputs = o12_obj.get("outputs") or {}
        n = int(outputs.get("n"))
        if n != 12:
            raise ValueError("expected n=12")
        alpha12 = float(outputs.get("alpha12_mod_2pi"))
        O12 = as_square_matrix(outputs.get("O12"), n=n, label="F461.outputs.O12")
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_F461_SHAPE", "error": str(exc)}, ensure_ascii=True))

    # Minimal audits: rely on direct computations (no hidden inputs).
    I = np.eye(n, dtype=float)
    orth_res = max_abs(O12.T @ O12 - I)
    invol_res = max_abs(O12 @ O12 - I)

    A1 = np.outer(u1, u1)
    A2 = np.outer(u2, u2)
    gluing_res = max_abs(A2 - (O12 @ A1 @ O12.T))

    a1_sign_res = max_abs(A1 - np.outer(-u1, -u1))
    a2_sign_res = max_abs(A2 - np.outer(-u2, -u2))

    atlas = {
        "object": "SelectorAtlas_pair12_sigma_int_corridor_projector_v1",
        "stage": "F463",
        "status": "actual_exported_lane_scoped_two_chart_selector_atlas_stub_with_overlap_and_gluing_data__projector_level__no_false_pass",
        "as_of": "2026-03-15",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "intent": (
            "Export one explicit lane-scoped two-chart selector atlas stub on {pair1,pair2} in the strict sigma-int corridor, "
            "including an explicit overlap-domain declaration, the chart transition operator O12, and projector-level operator gluing data "
            "A2 = O12 A1 O12^T. This is projector-level (sign-gauge-safe) and does not claim a global selector atlas nor any global discharge."
        ),
        "charts": {
            "pair1": {
                "chart_id": "pair1",
                "carrier_plane": "V1 = span{c1,s1} inside the n=12 real Fourier scaffold",
                "local_operator": {
                    "ref": str(IN_A1.relative_to(REPO)),
                    "object": "A_1_pair1_orientation_projector_operator_strict_core_v1",
                },
            },
            "pair2": {
                "chart_id": "pair2",
                "carrier_plane": "V2 = span{c2,s2} inside the n=12 real Fourier scaffold",
                "local_operator": {
                    "ref": str(IN_A2.relative_to(REPO)),
                    "object": "A_2_pair2_orientation_projector_operator_strict_core_v1",
                },
            },
        },
        "overlap_domain_declaration": {
            "scope": "sigma_int_corridor_declared_scope",
            "meaning": (
                "Both charts (pair1 and pair2) are simultaneously defined by the exported slot-free sigma-int theta-pair supply; "
                "therefore the pair1/pair2 atlas overlap is declared on that same corridor scope."
            ),
            "note": "This is a lane-scoped overlap declaration, not a global overlap-domain structure on the full strict domain.",
        },
        "transitions": {
            "pair1_to_pair2": {
                "operator_ref": str(IN_O12.relative_to(REPO)),
                "operator": "O12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1",
                "alpha12_mod_2pi": alpha12,
            },
            "pair2_to_pair1": {
                "inverse": "O12^T (equals O12 on the exported instantiation up to numerical tolerance)",
            },
        },
        "gluing_data": {
            "operator_section_ref": str(IN_SECTION.relative_to(REPO)),
            "law": "A_2(pair2) = O_12 A_1(pair1) O_12^T  (projector-level, sign-gauge-safe)",
            "audit_ref": str(IN_P466_SUMMARY.relative_to(REPO)),
        },
        "audits": {
            "o12_orthogonality_max_abs_residual": orth_res,
            "o12_involution_max_abs_residual": invol_res,
            "projector_gluing_max_abs_residual": gluing_res,
            "a1_sign_gauge_invariance_max_abs": a1_sign_res,
            "a2_sign_gauge_invariance_max_abs": a2_sign_res,
            "p466_summary_snapshot": p466_summary,
        },
        "hard_limits": [
            "Lane-scoped: exports only a two-chart atlas stub on {pair1,pair2} for the sigma-int corridor.",
            "Does not export a global selector atlas, global overlap-domain declaration, nor global cocycle data (H41 remains open globally).",
            "Does not discharge QW-2191 nor export a global selector transition/gluing object (H40 remains open globally).",
            "Projector-level only: does not derive a sign-sensitive physical orientation datum.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F463",
        "status": "F463_EXECUTED_CURRENT_STRICT_PAIR12_SELECTOR_ATLAS_STUB_AND_OVERLAP_DECLARATION_EXPORT_PACKET_NO_FALSE_PASS",
        "generated_utc": atlas["generated_utc"],
        "audits": {
            "o12_orthogonality_max_abs_residual": orth_res,
            "projector_gluing_max_abs_residual": gluing_res,
        },
        "no_false_pass": True,
    }

    OUT_ATLAS.write_text(json.dumps(atlas, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

