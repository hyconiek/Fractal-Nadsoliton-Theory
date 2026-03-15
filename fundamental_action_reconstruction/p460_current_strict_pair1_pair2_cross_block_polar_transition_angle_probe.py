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

IN_M_CONTROL = GENERATED / "m_control_residual_diag_declared_value_instantiated_v1.json"
IN_ALPHA12 = (
    GENERATED
    / "alpha12_pair1_pair2_transition_angle_strict_derived_from_sigma_int_slot_free_theta_pair_v1.json"
)

OUT_JSON = GENERATED / "p460_current_strict_pair1_pair2_cross_block_polar_transition_angle_probe.json"
OUT_SUMMARY = GENERATED / "p460_current_strict_pair1_pair2_cross_block_polar_transition_angle_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_finite_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def require_2x2(name: str, value: Any, missing: list[str]) -> np.ndarray | None:
    if not (
        isinstance(value, list)
        and len(value) == 2
        and all(isinstance(row, list) and len(row) == 2 and all(is_finite_number(x) for x in row) for row in value)
    ):
        missing.append(f"{name} (2x2 finite numeric matrix)")
        return None
    return np.array(value, dtype=float)


def wrap_angle_mod_2pi(theta: float) -> float:
    t = float(theta) % (2.0 * math.pi)
    if t < 0.0:
        t += 2.0 * math.pi
    return t


def min_angle_diff_mod(modulus: float, a: float, b: float) -> float:
    da = (float(a) % modulus) - (float(b) % modulus)
    da = (da + 0.5 * modulus) % modulus - 0.5 * modulus
    return float(abs(da))


def polar_orthogonal_factor(b: np.ndarray, *, tol: float) -> tuple[np.ndarray, dict[str, Any]]:
    # Compute Q in the polar decomposition B = Q S, via Q = B (B^T B)^(-1/2).
    m = b.T @ b
    evals, evecs = np.linalg.eigh(m)
    evals = np.array(evals, dtype=float)
    evecs = np.array(evecs, dtype=float)

    if evals.min() <= tol:
        raise ValueError("cross-block is singular/ill-conditioned; cannot compute (B^T B)^(-1/2) safely")

    inv_sqrt = np.diag(1.0 / np.sqrt(evals))
    m_inv_sqrt = evecs @ inv_sqrt @ evecs.T
    q = b @ m_inv_sqrt

    residual = float(np.linalg.norm(q.T @ q - np.eye(2)))
    det_q = float(np.linalg.det(q))
    meta = {
        "btb_eigvals": [float(x) for x in evals.tolist()],
        "q_orthonormal_residual": residual,
        "det_q": det_q,
    }
    return q, meta


def rotation_angle_if_proper(q: np.ndarray) -> float | None:
    det_q = float(np.linalg.det(q))
    if det_q < 0.0:
        return None
    return float(math.atan2(float(q[1, 0]), float(q[0, 0])))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing_files: list[str] = []
    for p in (IN_M_CONTROL, IN_ALPHA12):
        if not p.exists():
            missing_files.append(str(p.relative_to(REPO)))

    if missing_files:
        summary = {
            "stage": "P460",
            "status": "NOT_COMPUTABLE_MISSING_REQUIRED_FILES",
            "missing_files": missing_files,
            "theorem_level_pass": False,
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    mobj = load_json(IN_M_CONTROL)
    a12 = load_json(IN_ALPHA12)

    missing: list[str] = []

    cross_12 = require_2x2("M_control.cross_pair1_to_pair2_block", (mobj.get("computed") or {}).get("cross_pair1_to_pair2_block"), missing)
    alpha12_mod_2pi = (a12.get("outputs") or {}).get("alpha_12_mod_2pi")
    alpha12_mod_pi = (a12.get("outputs") or {}).get("alpha_12_mod_pi")
    if not is_finite_number(alpha12_mod_2pi):
        missing.append("F457.outputs.alpha_12_mod_2pi (finite number)")
    if not is_finite_number(alpha12_mod_pi):
        missing.append("F457.outputs.alpha_12_mod_pi (finite number)")

    if missing:
        summary = {
            "stage": "P460",
            "status": "NOT_COMPUTABLE_MISSING_OR_INVALID_INPUTS",
            "missing_or_invalid": sorted(set(missing)),
            "theorem_level_pass": False,
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    assert cross_12 is not None
    alpha12_mod_2pi_f = float(alpha12_mod_2pi)
    alpha12_mod_pi_f = float(alpha12_mod_pi)

    u, svals, vt = np.linalg.svd(cross_12)
    svals = np.array(svals, dtype=float)
    sigma_1 = float(max(svals))
    sigma_2 = float(min(svals))
    anisotropy_ratio = (sigma_1 / sigma_2) if sigma_2 != 0.0 else float("inf")

    tol = 1e-12
    q, qmeta = polar_orthogonal_factor(cross_12, tol=tol)
    alpha_cross = rotation_angle_if_proper(q)
    alpha_cross_mod_2pi = wrap_angle_mod_2pi(alpha_cross) if alpha_cross is not None else None
    alpha_cross_mod_pi = (float(alpha_cross_mod_2pi) % math.pi) if alpha_cross_mod_2pi is not None else None

    comparison = {
        "alpha12_mod_2pi_from_theta_supply": alpha12_mod_2pi_f,
        "alpha12_mod_pi_axis_only": alpha12_mod_pi_f,
        "alpha_cross_mod_2pi_from_polar_factor": alpha_cross_mod_2pi,
        "alpha_cross_mod_pi_axis_only": alpha_cross_mod_pi,
        "delta_mod_2pi_if_both": (
            min_angle_diff_mod(2.0 * math.pi, alpha12_mod_2pi_f, float(alpha_cross_mod_2pi))
            if alpha_cross_mod_2pi is not None
            else None
        ),
        "delta_mod_pi_if_both": (
            min_angle_diff_mod(math.pi, alpha12_mod_pi_f, float(alpha_cross_mod_pi))
            if alpha_cross_mod_pi is not None
            else None
        ),
    }

    artifact = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P460",
        "goal": "extract_pair1_pair2_cross_block_from_declared_control_pullback_and_compute_polar_orthogonal_factor_transition_angle",
        "inputs": {
            "m_control_value_instantiated": str(IN_M_CONTROL.relative_to(REPO)),
            "alpha12_from_theta_supply": str(IN_ALPHA12.relative_to(REPO)),
        },
        "extracted": {
            "cross_pair1_to_pair2_block": [[float(x) for x in row] for row in cross_12.tolist()],
            "svd": {
                "singular_values": [float(x) for x in svals.tolist()],
                "anisotropy_ratio": anisotropy_ratio,
            },
        },
        "polar": {
            "Q": [[float(x) for x in row] for row in q.tolist()],
            "btb_eigvals": qmeta["btb_eigvals"],
            "q_orthonormal_residual": qmeta["q_orthonormal_residual"],
            "det_q": qmeta["det_q"],
            "alpha_cross_mod_2pi_if_det_positive": alpha_cross_mod_2pi,
            "alpha_cross_mod_pi_axis_only_if_det_positive": alpha_cross_mod_pi,
            "tolerance": tol,
        },
        "comparison_to_theta_supply": comparison,
        "hard_limits": [
            "Consumes P459 which uses the conditional N477 Yukawa-free diagonal residual rewrite; no stationarity witness is exported.",
            "Does not claim host matching/cancellation.",
            "Does not claim the extracted polar factor defines a global selector transition/gluing object.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim global discharge of QW-2191.",
            "Does not assign a theorem-level physical interpretation to alpha_cross.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P460",
        "status": "PASS_PROBE_READY_CROSS_BLOCK_POLAR_FACTOR_COMPUTED" if alpha_cross_mod_2pi is not None else "PASS_PROBE_READY_DET_NEGATIVE_REFLECTION_CASE",
        "alpha_cross_mod_2pi": alpha_cross_mod_2pi,
        "alpha_cross_mod_pi": alpha_cross_mod_pi,
        "alpha12_mod_2pi": alpha12_mod_2pi_f,
        "alpha12_mod_pi": alpha12_mod_pi_f,
        "cross_block_singular_values": [sigma_1, sigma_2],
        "q_orthonormal_residual": qmeta["q_orthonormal_residual"],
        "det_q": qmeta["det_q"],
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

