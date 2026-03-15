#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F460 = GENERATED / "psi_hessian_diagonal_local_strict_derived_eigenprojectors_v1.json"
IN_F453 = GENERATED / "mode_index_assignment_canonical_local_diagonal_strict_derived_v1.json"

OUT_JSON = (
    GENERATED
    / "p464_current_strict_diagonal_local_psi_hessian_eigenprojectors_to_pair_plane_weight_profile_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p464_current_strict_diagonal_local_psi_hessian_eigenprojectors_to_pair_plane_weight_profile_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def max_abs(a: np.ndarray) -> float:
    return float(np.max(np.abs(a)))


def as_vector(x: Any, n: int, label: str) -> np.ndarray:
    if not (isinstance(x, list) and len(x) == n and all(isinstance(v, (int, float)) for v in x)):
        raise ValueError(f"{label} must be a length-{n} numeric list")
    return np.array([float(v) for v in x], dtype=float)


def as_square_matrix(x: Any, n: int, label: str) -> np.ndarray:
    if not (isinstance(x, list) and len(x) == n and all(isinstance(row, list) and len(row) == n for row in x)):
        raise ValueError(f"{label} must be a {n}x{n} nested list")
    return np.array(x, dtype=float)


def shannon_entropy(weights: list[float], eps: float = 0.0) -> float:
    h = 0.0
    for w in weights:
        ww = float(w)
        if ww <= eps:
            continue
        h -= ww * math.log(ww)
    return float(h)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing: list[str] = []
    if not IN_F460.exists():
        missing.append(str(IN_F460.relative_to(REPO)))
    if not IN_F453.exists():
        missing.append(str(IN_F453.relative_to(REPO)))
    if missing:
        raise SystemExit(
            json.dumps(
                {
                    "status": "NOT_COMPUTABLE_MISSING_INPUTS",
                    "missing": missing,
                    "expected": ["F460 eigenprojectors export", "F453 diagonal/local mode-index assignment export"],
                },
                ensure_ascii=True,
            )
        )

    f460 = load_json(IN_F460)
    f453 = load_json(IN_F453)

    try:
        projectors_raw = ((f460.get("outputs") or {}).get("eigenprojectors_rank1")) or []
        evals_raw = ((f460.get("outputs") or {}).get("eigenvalues_ascending")) or []
        if not (isinstance(projectors_raw, list) and len(projectors_raw) == 12):
            raise ValueError("F460.outputs.eigenprojectors_rank1 must be length-12 list")
        if not (isinstance(evals_raw, list) and len(evals_raw) == 12):
            raise ValueError("F460.outputs.eigenvalues_ascending must be length-12 list")
        p_list = [as_square_matrix(pj, n=12, label=f"F460.outputs.eigenprojectors_rank1[{j}]") for j, pj in enumerate(projectors_raw)]
        evals = [float(v) for v in evals_raw]
    except Exception as exc:
        raise SystemExit(
            json.dumps(
                {
                    "status": "INVALID_F460_SHAPE",
                    "expected": "F460.outputs.{eigenvalues_ascending,eigenprojectors_rank1}",
                    "error": str(exc),
                },
                ensure_ascii=True,
            )
        )

    try:
        outputs = f453.get("outputs") or {}
        e0 = as_vector(outputs.get("e0"), n=12, label="F453.outputs.e0")
        e6 = as_vector(outputs.get("e6"), n=12, label="F453.outputs.e6")
        pairs = (outputs.get("pairs") or {})
        if not isinstance(pairs, dict):
            raise ValueError("F453.outputs.pairs must be an object")
        pair_labels = [f"pair{m}" for m in range(1, 6)]
        pair_basis: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for label in pair_labels:
            entry = pairs.get(label) or {}
            u_plus = as_vector(entry.get("u_plus"), n=12, label=f"F453.outputs.pairs.{label}.u_plus")
            u_minus = as_vector(entry.get("u_minus"), n=12, label=f"F453.outputs.pairs.{label}.u_minus")
            pair_basis[label] = (u_plus, u_minus)
    except Exception as exc:
        raise SystemExit(
            json.dumps(
                {
                    "status": "INVALID_F453_SHAPE",
                    "expected": "F453.outputs.{e0,e6,pairs.pairm.u_plus,u_minus}",
                    "error": str(exc),
                },
                ensure_ascii=True,
            )
        )

    # Build plane projectors Π_label using orthonormal basis vectors in each plane.
    labels = ["e0", "pair1", "pair2", "pair3", "pair4", "pair5", "e6"]
    plane_projectors: dict[str, np.ndarray] = {"e0": np.outer(e0, e0), "e6": np.outer(e6, e6)}
    plane_orthonormal_residuals: dict[str, float] = {
        "e0_norm_residual": abs(float(e0 @ e0) - 1.0),
        "e6_norm_residual": abs(float(e6 @ e6) - 1.0),
        "e0_dot_e6_abs": abs(float(e0 @ e6)),
    }

    for label in [f"pair{m}" for m in range(1, 6)]:
        u_plus, u_minus = pair_basis[label]
        B = np.column_stack([u_plus, u_minus])
        plane_projectors[label] = B @ B.T
        gram_res = B.T @ B - np.eye(2)
        plane_orthonormal_residuals[f"{label}_orthonormal_max_abs_residual"] = max_abs(gram_res)
        plane_orthonormal_residuals[f"{label}_uplus_norm_residual"] = abs(float(u_plus @ u_plus) - 1.0)
        plane_orthonormal_residuals[f"{label}_uminus_norm_residual"] = abs(float(u_minus @ u_minus) - 1.0)
        plane_orthonormal_residuals[f"{label}_uplus_dot_uminus_abs"] = abs(float(u_plus @ u_minus))

    projector_sum = np.zeros((12, 12), dtype=float)
    for label in labels:
        projector_sum += plane_projectors[label]
    projector_sum_residual = max_abs(projector_sum - np.eye(12))

    # Weights: w_{j,label} = tr(Π_label P_j).
    weights_by_j: list[dict[str, float]] = []
    per_eigenmode_compact: list[dict[str, Any]] = []
    weight_sum_residuals: list[float] = []

    for j in range(12):
        pj = p_list[j]
        row: dict[str, float] = {}
        for label in labels:
            w = float(np.trace(plane_projectors[label] @ pj))
            # Numeric hygiene: clip tiny negatives.
            if w < 0.0 and w > -1e-12:
                w = 0.0
            row[label] = w
        weights_by_j.append(row)
        s = float(sum(row[label] for label in labels))
        weight_sum_residuals.append(abs(s - 1.0))

        w_list = [row[label] for label in labels]
        top_label = max(labels, key=lambda lab: row[lab])
        top_w = float(row[top_label])
        per_eigenmode_compact.append(
            {
                "j": j,
                "lambda": float(evals[j]),
                "top_plane": {"label": top_label, "mass": top_w},
                "max_plane_weight": top_w,
                "entropy_nats": shannon_entropy(w_list, eps=0.0),
            }
        )

    artifact = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P464",
        "goal": "export_pair_plane_weight_profile_for_f460_psi_hessian_eigenprojectors_using_f453_pair_plane_projectors",
        "inputs": {
            "f460_eigenprojectors": str(IN_F460.relative_to(REPO)),
            "f453_mode_index_assignment": str(IN_F453.relative_to(REPO)),
        },
        "labels_order": labels,
        "audits": {
            "plane_projector_sum_max_abs_residual_vs_I12": projector_sum_residual,
            "per_eigenprojector_weight_sum_abs_residual_max": float(max(weight_sum_residuals)),
            "per_eigenprojector_weight_sum_abs_residual_mean": float(sum(weight_sum_residuals) / len(weight_sum_residuals)),
            "plane_basis_orthonormality": plane_orthonormal_residuals,
        },
        "weights_by_eigenmode": weights_by_j,
        "per_eigenmode_compact": per_eigenmode_compact,
        "hard_limits": [
            "Does not claim F453 diagonalizes the full H_psi; this exports only plane weights of the F459/F460 eigensystem.",
            "Does not export any global selector transition/gluing object.",
            "Does not claim global discharge of QW-2191 beyond the declared lane-scoped instantiations.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    max_max_weight = max(item["max_plane_weight"] for item in per_eigenmode_compact)
    min_max_weight = min(item["max_plane_weight"] for item in per_eigenmode_compact)
    mean_entropy = float(sum(item["entropy_nats"] for item in per_eigenmode_compact) / len(per_eigenmode_compact))

    summary = {
        "stage": "P464",
        "status": "PASS_PROBE_READY",
        "generated_utc": artifact["generated_utc"],
        "labels_order": labels,
        "audits": {
            "plane_projector_sum_max_abs_residual_vs_I12": projector_sum_residual,
            "per_eigenprojector_weight_sum_abs_residual_max": artifact["audits"]["per_eigenprojector_weight_sum_abs_residual_max"],
        },
        "mixing_metrics": {
            "max_over_eigenmodes_of_max_plane_weight": float(max_max_weight),
            "min_over_eigenmodes_of_max_plane_weight": float(min_max_weight),
            "mean_entropy_nats_over_eigenmodes": mean_entropy,
        },
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

