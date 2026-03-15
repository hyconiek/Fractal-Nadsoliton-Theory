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

IN_F459 = GENERATED / "psi_hessian_diagonal_local_strict_derived_value_instantiated_v1.json"

OUT_PROJECTORS = GENERATED / "psi_hessian_diagonal_local_strict_derived_eigenprojectors_v1.json"
OUT_SUMMARY = GENERATED / "psi_hessian_diagonal_local_strict_derived_eigenprojectors_v1_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def max_abs(a: np.ndarray) -> float:
    return float(np.max(np.abs(a)))


def as_square_matrix(x: Any, n: int, label: str) -> np.ndarray:
    if not (isinstance(x, list) and len(x) == n and all(isinstance(row, list) and len(row) == n for row in x)):
        raise ValueError(f"{label} must be a {n}x{n} nested list")
    return np.array(x, dtype=float)


def as_vector_list(x: Any, n: int, label: str) -> list[np.ndarray]:
    if not (isinstance(x, list) and len(x) == n and all(isinstance(v, list) and len(v) == n for v in x)):
        raise ValueError(f"{label} must be a length-{n} list of length-{n} lists")
    return [np.array(v, dtype=float) for v in x]


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_F459.exists():
        raise SystemExit(
            json.dumps(
                {
                    "status": "NOT_COMPUTABLE_MISSING_INPUT",
                    "missing": str(IN_F459.relative_to(REPO)),
                    "expected": "F459 Psi Hessian eigensystem value instantiation export",
                },
                ensure_ascii=True,
            )
        )

    f459 = load_json(IN_F459)
    outputs = f459.get("outputs") or {}

    try:
        n = int(outputs.get("n"))
        if n != 12:
            raise ValueError("expected n=12")
        h_psi = as_square_matrix(outputs.get("h_psi"), n=n, label="F459.outputs.h_psi")
        eigenvalues = outputs.get("eigenvalues_ascending")
        if not (isinstance(eigenvalues, list) and len(eigenvalues) == n):
            raise ValueError("F459.outputs.eigenvalues_ascending must be length-12 list")
        evals = np.array([float(v) for v in eigenvalues], dtype=float)
        evec_columns = as_vector_list(outputs.get("eigenvectors_columns"), n=n, label="F459.outputs.eigenvectors_columns")
    except Exception as exc:
        raise SystemExit(
            json.dumps(
                {
                    "status": "INVALID_INPUT_SHAPE",
                    "expected": "F459.outputs.{n,h_psi,eigenvalues_ascending,eigenvectors_columns}",
                    "error": str(exc),
                },
                ensure_ascii=True,
            )
        )

    V = np.column_stack(evec_columns)  # columns are eigenvectors
    overlap = V.T @ V
    overlap_residual = overlap - np.eye(n)

    projectors: list[np.ndarray] = []
    idempotent_residuals: list[float] = []
    symmetry_residuals: list[float] = []
    trace_residuals: list[float] = []

    for j in range(n):
        vj = V[:, j]
        pj = np.outer(vj, vj)
        projectors.append(pj)
        idempotent_residuals.append(max_abs(pj @ pj - pj))
        symmetry_residuals.append(max_abs(pj - pj.T))
        trace_residuals.append(abs(float(np.trace(pj)) - 1.0))

    p_sum = np.zeros((n, n), dtype=float)
    for pj in projectors:
        p_sum += pj

    h_rec = np.zeros((n, n), dtype=float)
    for j in range(n):
        h_rec += float(evals[j]) * projectors[j]

    max_pair_projector_overlap = 0.0
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            max_pair_projector_overlap = max(max_pair_projector_overlap, max_abs(projectors[i] @ projectors[j]))

    artifact = {
        "object": "PsiHessian_eigenprojectors_diagonal_local_strict_derived_v1",
        "stage": "F460",
        "status": "actual_exported_strict_derived_eigenprojectors__psi_hessian_spectral_packaging__no_false_pass",
        "as_of": "2026-03-15",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": "export_sign_gauge_invariant_rank_one_spectral_projectors_for_the_f459_diagonal_local_psi_hessian",
        "inputs": {
            "f459_eigensystem_ref": str(IN_F459.relative_to(REPO)),
            "eigenvector_sign_note": "F459 eigenvectors are representatives with an explicit sign convention; this packet exports sign-invariant projectors P_j := |v_j><v_j|.",
        },
        "construction": {
            "projectors": "P_j := v_j v_j^T  (rank-one orthogonal projectors for orthonormal eigenvectors)",
            "spectral_reconstruction": "H_psi_reconstructed := Σ_j λ_j P_j",
        },
        "outputs": {
            "n": n,
            "eigenvalues_ascending": [float(v) for v in evals.tolist()],
            "eigenprojectors_rank1": [[[float(x) for x in row] for row in pj.tolist()] for pj in projectors],
            "checks": {
                "evecs_orthonormal_max_abs_residual": max_abs(overlap_residual),
                "projector_sum_max_abs_residual_vs_I": max_abs(p_sum - np.eye(n)),
                "spectral_reconstruction_max_abs_residual": max_abs(h_psi - h_rec),
                "projector_idempotent_max_abs_residual_max": float(max(idempotent_residuals)),
                "projector_symmetry_max_abs_residual_max": float(max(symmetry_residuals)),
                "projector_trace_residual_max": float(max(trace_residuals)),
                "projector_pairwise_overlap_max_abs": float(max_pair_projector_overlap),
            },
        },
        "hard_limits": [
            "This is downstream spectral packaging only: it does not claim host matching/cancellation.",
            "Does not claim any global discharge of QW-2191.",
            "Does not derive a sign-sensitive physical orientation convention; projectors are sign-gauge-invariant.",
            "Does not claim admissible S_sel_int nor strict-core selector closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F460",
        "status": "F460_EXECUTED_CURRENT_STRICT_DIAGONAL_LOCAL_PSI_HESSIAN_EIGENPROJECTORS_EXPORT_PACKET_NO_FALSE_PASS",
        "generated_utc": artifact["generated_utc"],
        "inputs": {"f459_ref": str(IN_F459.relative_to(REPO))},
        "outputs": {
            "n": n,
            "min_eigenvalue": float(np.min(evals)),
            "max_eigenvalue": float(np.max(evals)),
            "spectral_reconstruction_max_abs_residual": artifact["outputs"]["checks"]["spectral_reconstruction_max_abs_residual"],
            "projector_sum_max_abs_residual_vs_I": artifact["outputs"]["checks"]["projector_sum_max_abs_residual_vs_I"],
        },
        "no_false_pass": True,
    }

    OUT_PROJECTORS.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

