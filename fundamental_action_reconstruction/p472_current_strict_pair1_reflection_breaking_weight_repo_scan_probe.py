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

IN_U1 = GENERATED / "a_1_pair1_orientation_projector_operator_strict_core_v1.json"

OUT_JSON = GENERATED / "p472_current_strict_pair1_reflection_breaking_weight_repo_scan_probe.json"
OUT_SUMMARY = GENERATED / "p472_current_strict_pair1_reflection_breaking_weight_repo_scan_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_numeric_list_len(obj: Any, n: int) -> bool:
    return (
        isinstance(obj, list)
        and len(obj) == n
        and all(isinstance(v, (int, float)) and math.isfinite(float(v)) for v in obj)
    )


def walk_collect_len12_arrays(obj: Any, path: list[str], out: list[dict[str, Any]]) -> None:
    if is_numeric_list_len(obj, 12):
        out.append({"path": ".".join(path) if path else "<root>", "values": [float(v) for v in obj]})
        return
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(k, str):
                walk_collect_len12_arrays(v, path + [k], out)
            else:
                walk_collect_len12_arrays(v, path + [str(k)], out)
        return
    if isinstance(obj, list):
        for i, v in enumerate(obj):
            walk_collect_len12_arrays(v, path + [f"[{i}]"], out)


def is_weight_like_path(path: str) -> bool:
    p = path.lower()
    bad_tokens = [
        "eigen",
        "eigval",
        "eigvec",
        "eigensystem",
        "spectrum",
        "eigenvalues",
        "eigenvectors",
        "columns",
        "matrix",
        "projector",
        "operator",
        "transition",
        "atlas",
        "glued",
        "basis",
        "mode_index_assignment",
        "coords",
        "theta",
        "u_",
        "u1",
        "uplus",
        "uminus",
        "pair",
    ]
    if any(tok in p for tok in bad_tokens):
        return False
    good_tokens = [
        "vpsi",
        "g4",
        "g6",
        "d_local",
        "residual_profile",
        "ord",
        "r_ord",
        "weight",
        "weights",
        "reference",
        "distribution",
        "ktotal",
        "k_total",
        "distance",
    ]
    return any(tok in p for tok in good_tokens)


def classify_scope(meta: dict[str, Any], filename: str) -> str:
    classification = str(meta.get("classification") or "")
    status = str(meta.get("status") or "")
    object_id = str(meta.get("object") or meta.get("object_id") or "")
    s = " ".join([classification, status, object_id, filename]).lower()
    if "strict_extension" in s or filename.lower().startswith("ax") or "axiom_lane" in s:
        return "strict_extension_only_or_axiom_lane"
    if "strict_derived" in s or "strict-derived" in s or "strict derived" in s or "strict_derived" in classification.lower():
        return "strict_derived"
    if "strict_core" in s or "strict core" in s:
        return "strict_core"
    return "unknown_or_mixed"


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_U1.exists():
        raise SystemExit(
            json.dumps(
                {
                    "stage": "P472",
                    "status": "NOT_COMPUTABLE_MISSING_INPUTS",
                    "missing": [str(IN_U1.relative_to(REPO))],
                    "expected": ["F456 exports A_1(pair1) with u_1 (strict core)"],
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    a1 = load_json(IN_U1)
    try:
        u1_raw = (a1.get("data") or {}).get("u_1")
        if not is_numeric_list_len(u1_raw, 12):
            raise ValueError("A1.data.u_1 must be a length-12 numeric list")
        u1 = [float(v) for v in u1_raw]
    except Exception as exc:
        raise SystemExit(
            json.dumps(
                {"stage": "P472", "status": "INVALID_U1_INPUT_SHAPE", "error": str(exc), "no_false_pass": True},
                ensure_ascii=True,
            )
        )

    n = 12

    def reflect(x: int) -> int:
        return (-x) % n

    tol_even = 1e-12
    tol_dot = 1e-9

    json_files = sorted(p for p in GENERATED.glob("*.json") if p.is_file())

    candidates: list[dict[str, Any]] = []
    candidates_weight_like: list[dict[str, Any]] = []
    for path in json_files:
        try:
            obj = load_json(path)
        except Exception:
            continue

        arrays: list[dict[str, Any]] = []
        walk_collect_len12_arrays(obj, [], arrays)
        if not arrays:
            continue

        scope = classify_scope(obj if isinstance(obj, dict) else {}, filename=path.name)
        meta = {
            "file": str(path.relative_to(REPO)),
            "scope": scope,
            "object": (str(obj.get("object")) if isinstance(obj, dict) and obj.get("object") else None),
            "object_id": (str(obj.get("object_id")) if isinstance(obj, dict) and obj.get("object_id") else None),
            "status": (str(obj.get("status")) if isinstance(obj, dict) and obj.get("status") else None),
            "classification": (
                str(obj.get("classification")) if isinstance(obj, dict) and obj.get("classification") else None
            ),
        }

        for arr in arrays:
            w = arr["values"]
            max_abs_even = max(abs(w[x] - w[reflect(x)]) for x in range(n))
            dot = sum(w[x] * u1[x] for x in range(n))
            record = {
                **meta,
                "array_path": arr["path"],
                "reflection_even_max_abs": float(max_abs_even),
                "reflection_even": bool(max_abs_even <= tol_even),
                "dot_w_u1": float(dot),
                "dot_nonzero": bool(abs(float(dot)) > tol_dot),
            }
            candidates.append(record)
            if is_weight_like_path(arr["path"]):
                candidates_weight_like.append(record)

    def is_strictish(c: dict[str, Any]) -> bool:
        return c["scope"] in ("strict_core", "strict_derived")

    strictish = [c for c in candidates if is_strictish(c)]
    strictish_breaking = [c for c in strictish if not bool(c["reflection_even"])]
    strictish_breaking_nonzero = [c for c in strictish_breaking if bool(c["dot_nonzero"])]

    extension = [c for c in candidates if c["scope"] == "strict_extension_only_or_axiom_lane"]
    extension_breaking_nonzero = [c for c in extension if (not bool(c["reflection_even"])) and bool(c["dot_nonzero"])]

    strictish_w = [c for c in candidates_weight_like if is_strictish(c)]
    strictish_w_breaking = [c for c in strictish_w if not bool(c["reflection_even"])]
    strictish_w_breaking_nonzero = [c for c in strictish_w_breaking if bool(c["dot_nonzero"])]

    extension_w = [c for c in candidates_weight_like if c["scope"] == "strict_extension_only_or_axiom_lane"]
    extension_w_breaking_nonzero = [
        c for c in extension_w if (not bool(c["reflection_even"])) and bool(c["dot_nonzero"])
    ]

    strictish_w_breaking_nonzero_outside_k_total_rows = [
        c for c in strictish_w_breaking_nonzero if "outputs.k_total" not in str(c.get("array_path") or "")
    ]

    # Sort by absolute dot magnitude, descending (most promising sign-sensitive scalar candidates).
    strictish_breaking_nonzero_sorted = sorted(strictish_breaking_nonzero, key=lambda c: abs(float(c["dot_w_u1"])), reverse=True)
    extension_breaking_nonzero_sorted = sorted(extension_breaking_nonzero, key=lambda c: abs(float(c["dot_w_u1"])), reverse=True)
    strictish_w_breaking_nonzero_sorted = sorted(
        strictish_w_breaking_nonzero, key=lambda c: abs(float(c["dot_w_u1"])), reverse=True
    )
    extension_w_breaking_nonzero_sorted = sorted(
        extension_w_breaking_nonzero, key=lambda c: abs(float(c["dot_w_u1"])), reverse=True
    )

    artifact = {
        "stage": "P472",
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": (
            "scan_generated_repo_artifacts_for_site_indexed_length_12_numeric_arrays_that_break_reflection_symmetry "
            "and therefore could, in principle, distinguish u_1 from -u_1 on pair1 via a scalar Σ_x w(x)u_1(x), "
            "without smuggling a marked direction"
        ),
        "inputs": {
            "u1_ref": str(IN_U1.relative_to(REPO)),
            "generated_dir": str(GENERATED.relative_to(REPO)),
            "tolerances": {"reflection_even_max_abs_tol": tol_even, "dot_nonzero_abs_tol": tol_dot},
        },
        "scan": {
            "files_scanned": len(json_files),
            "candidates_total": len(candidates),
            "candidates_strictish_total": len(strictish),
            "candidates_strictish_reflection_breaking": len(strictish_breaking),
            "candidates_strictish_reflection_breaking_and_dot_nonzero": len(strictish_breaking_nonzero),
            "candidates_extension_reflection_breaking_and_dot_nonzero": len(extension_breaking_nonzero),
            "candidates_weight_like_total": len(candidates_weight_like),
            "candidates_weight_like_strictish_total": len(strictish_w),
            "candidates_weight_like_strictish_reflection_breaking": len(strictish_w_breaking),
            "candidates_weight_like_strictish_reflection_breaking_and_dot_nonzero": len(strictish_w_breaking_nonzero),
            "candidates_weight_like_strictish_reflection_breaking_and_dot_nonzero_outside_k_total_rows": len(
                strictish_w_breaking_nonzero_outside_k_total_rows
            ),
            "candidates_weight_like_extension_reflection_breaking_and_dot_nonzero": len(extension_w_breaking_nonzero),
        },
        "top_candidates": {
            "strictish_reflection_breaking_and_dot_nonzero": strictish_breaking_nonzero_sorted[:25],
            "extension_reflection_breaking_and_dot_nonzero": extension_breaking_nonzero_sorted[:25],
            "weight_like_strictish_reflection_breaking_and_dot_nonzero": strictish_w_breaking_nonzero_sorted[:25],
            "weight_like_extension_reflection_breaking_and_dot_nonzero": extension_w_breaking_nonzero_sorted[:25],
            "weight_like_strictish_reflection_breaking_and_dot_nonzero_outside_k_total_rows": sorted(
                strictish_w_breaking_nonzero_outside_k_total_rows, key=lambda c: abs(float(c["dot_w_u1"])), reverse=True
            )[:25],
        },
        "no_false_pass": True,
        "hard_limits": [
            "Probe-only: does not discharge H37 and does not export a sign-sensitive physical orientation datum.",
            "Heuristic scope classifier only; strict promotion requires checking theorem/packet provenance.",
            "A nonzero dot product here is only a necessary indicator for a sign-sensitive scalar; it is not a physical interpretation claim.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim global discharge of QW-2191.",
            "Does not claim ToE closure.",
        ],
    }

    summary = {
        "stage": "P472",
        "status": "PASS_SCAN_COMPLETE",
        "candidates_strictish_reflection_breaking_and_dot_nonzero": len(strictish_breaking_nonzero),
        "candidates_extension_reflection_breaking_and_dot_nonzero": len(extension_breaking_nonzero),
        "candidates_weight_like_strictish_reflection_breaking_and_dot_nonzero": len(strictish_w_breaking_nonzero),
        "candidates_weight_like_strictish_reflection_breaking_and_dot_nonzero_outside_k_total_rows": len(
            strictish_w_breaking_nonzero_outside_k_total_rows
        ),
        "candidates_weight_like_extension_reflection_breaking_and_dot_nonzero": len(extension_w_breaking_nonzero),
        "recommended_next_strict_target": "H37",
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
