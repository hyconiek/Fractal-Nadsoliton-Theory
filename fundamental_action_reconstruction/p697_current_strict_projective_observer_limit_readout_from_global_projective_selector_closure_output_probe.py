#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

AS_OF = "2026-03-17"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_PROJECTIVE_CLOSURE_F672 = GENERATED / "selector_closure_global_c_v1_projective_strict_v1.json"

OUT_JSON = (
    GENERATED
    / "p697_current_strict_projective_observer_limit_readout_from_global_projective_selector_closure_output_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p697_current_strict_projective_observer_limit_readout_from_global_projective_selector_closure_output_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_finite_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def is_numeric_matrix(obj: Any, n: int) -> bool:
    return (
        isinstance(obj, list)
        and len(obj) == n
        and all(isinstance(row, list) and len(row) == n and all(is_finite_number(v) for v in row) for row in obj)
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_PROJECTIVE_CLOSURE_F672]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P697",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    closure = load_json(IN_PROJECTIVE_CLOSURE_F672)
    out_obs = closure.get("output_observable")
    if not isinstance(out_obs, dict):
        artifact = {
            "stage": "P697",
            "status": "INVALID_F672_OUTPUT_OBSERVABLE_SHAPE",
            "as_of": AS_OF,
            "error": "F672 must export output_observable as an object",
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    basis = out_obs.get("basis")
    if basis != ["o_+", "o_-"]:
        artifact = {
            "stage": "P697",
            "status": "NOT_COMPUTABLE_UNEXPECTED_OUTPUT_BASIS",
            "as_of": AS_OF,
            "expected_basis": ["o_+", "o_-"],
            "actual_basis": basis,
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    proj = out_obs.get("output_projector_matrix_in_o_plus_o_minus")
    if not is_numeric_matrix(proj, 2):
        artifact = {
            "stage": "P697",
            "status": "INVALID_F672_OUTPUT_PROJECTOR_SHAPE",
            "as_of": AS_OF,
            "error": "F672 output_observable.output_projector_matrix_in_o_plus_o_minus must be a 2x2 numeric matrix",
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    P = np.array(proj, dtype=float)
    trace = float(np.trace(P))
    sym_residual_inf = float(np.max(np.abs(P - P.T)))
    eig = np.linalg.eigvalsh(0.5 * (P + P.T))
    eig_sorted = [float(eig[0]), float(eig[1])]

    p_plus = float(P[0, 0])
    p_minus = float(P[1, 1])
    coherence_offdiag_max_abs = float(np.max(np.abs(P - np.diag(np.diag(P)))))

    # Minimal observer-limit readout invariants from the projective output observable.
    y_bias = 0.5 * (p_plus - p_minus)
    y_total = 0.5 * (p_plus + p_minus)
    z_commit = p_plus
    z_residual = p_minus

    tol = 1e-12
    ok_projector_like = abs(trace - 1.0) <= 1e-9 and sym_residual_inf <= 1e-9 and eig_sorted[0] >= -1e-10
    ok_readout = (z_commit > 0.0) and (z_residual >= -tol)

    status = "PASS_PROJECTIVE_OBSERVER_LIMIT_READOUT_COMPUTABLE_FROM_GLOBAL_PROJECTIVE_SELECTOR_CLOSURE_OUTPUT_PROJECTOR"
    if not (ok_projector_like and ok_readout):
        status = "FAIL_PROJECTIVE_OUTPUT_OBSERVABLE_NOT_PROJECTOR_LIKE_OR_READOUT_NOT_WELL_FORMED"

    artifact = {
        "stage": "P697",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": "compute_minimal_observer_limit_readout_invariants_from_global_projective_selector_closure_output_projector__treat_sign_as_gauge__no_false_pass",
        "inputs": {
            "projective_selector_closure_ref": str(IN_PROJECTIVE_CLOSURE_F672.relative_to(REPO)),
        },
        "audits": {
            "output_basis": list(basis),
            "output_projector_trace": trace,
            "output_projector_symmetry_residual_inf_norm": sym_residual_inf,
            "output_projector_eigenvalues_symmetrized": eig_sorted,
            "output_projector_offdiag_max_abs": coherence_offdiag_max_abs,
        },
        "results": {
            "p_plus": p_plus,
            "p_minus": p_minus,
            "y_bias": y_bias,
            "y_total": y_total,
            "z_commit": z_commit,
            "z_residual": z_residual,
        },
        "hard_limits": [
            "Projective/gauge-safe downstream readout only; does not claim any directed/sign-sensitive physical orientation datum in strict core.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim strict-core ToE closure.",
            "Does not claim actual emergent observer construction/closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P697",
        "status": status,
        "p_plus": p_plus,
        "p_minus": p_minus,
        "y_bias": y_bias,
        "y_total": y_total,
        "z_commit": z_commit,
        "z_residual": z_residual,
        "output_projector_trace": trace,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

