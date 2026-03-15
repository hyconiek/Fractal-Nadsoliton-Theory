#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

R14_K_TOTAL = GENERATED / "r14_explicit_frozen_kernel_specialization_packet_for_host_matching_route.json"
R15_M0 = GENERATED / "r15_explicit_host_scalar_floor_embedding_packet_for_host_matching_route.json"

P437_SCRIPT = ROOT / "p437_current_strict_r14_r15_n477_canonical_local_diagonal_sigma_opposite_pair_sums_value_evaluation_harness_probe.py"
P437_OUT = (
    GENERATED
    / "p437_current_strict_r14_r15_n477_canonical_local_diagonal_sigma_opposite_pair_sums_value_evaluation_harness_probe.json"
)
P437_SUMMARY = (
    GENERATED
    / "p437_current_strict_r14_r15_n477_canonical_local_diagonal_sigma_opposite_pair_sums_value_evaluation_harness_probe_summary.json"
)

OUT_JSON = GENERATED / "psi_hessian_diagonal_local_strict_derived_value_instantiated_v1.json"
OUT_SUMMARY = GENERATED / "psi_hessian_diagonal_local_strict_derived_value_instantiated_v1_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def run_script(path: Path) -> dict[str, Any]:
    proc = subprocess.run(
        ["python3", str(path)],
        cwd=str(REPO),
        capture_output=True,
        text=True,
        check=False,
    )
    return {
        "script": str(path.relative_to(REPO)),
        "returncode": proc.returncode,
        "stdout_tail": "\n".join(proc.stdout.splitlines()[-10:]),
        "stderr_tail": "\n".join(proc.stderr.splitlines()[-10:]),
    }


def extract_k_total_matrix(r14: dict[str, Any]) -> np.ndarray:
    rows = r14.get("specialization_rows")
    if not (isinstance(rows, list) and len(rows) == 12 and all(isinstance(r, list) and len(r) == 12 for r in rows)):
        raise ValueError("R14.specialization_rows must be a 12x12 list")
    m = np.zeros((12, 12), dtype=float)
    for i in range(12):
        for j in range(12):
            cell = rows[i][j]
            if not isinstance(cell, dict):
                raise ValueError("R14.specialization_rows entries must be dicts")
            v = cell.get("host_kernel_value", cell.get("specialized_value"))
            if not isinstance(v, (int, float)) or not math.isfinite(float(v)):
                raise ValueError("R14 matrix entry not finite numeric")
            m[i, j] = float(v)
    return m


def canonicalize_eigenvector_sign(v: np.ndarray, tol: float = 0.0) -> np.ndarray:
    if v.ndim != 1:
        raise ValueError("expected 1D eigenvector")
    idx = int(np.argmax(np.abs(v)))
    if float(abs(v[idx])) <= tol:
        return v
    if float(v[idx]) < 0.0:
        return -v
    return v


def group_degeneracies(evals: list[float], tol: float) -> list[dict[str, Any]]:
    groups: list[list[float]] = []
    for lam in evals:
        if not groups:
            groups.append([lam])
            continue
        if abs(lam - groups[-1][-1]) <= tol:
            groups[-1].append(lam)
        else:
            groups.append([lam])
    return [
        {"multiplicity": len(g), "lambda_min": float(min(g)), "lambda_max": float(max(g)), "span": float(max(g) - min(g))}
        for g in groups
    ]


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing: list[str] = []
    for p in (R14_K_TOTAL, R15_M0):
        if not p.exists():
            missing.append(str(p.relative_to(REPO)))
    if missing:
        raise SystemExit(json.dumps({"status": "MISSING_REQUIRED_INPUTS", "missing": missing}, ensure_ascii=True))

    # Ensure P437 is executed on the current strict-derived inputs (idempotent if unchanged).
    run = run_script(P437_SCRIPT)
    if int(run["returncode"]) != 0:
        raise SystemExit(json.dumps({"status": "DEPENDENCY_SCRIPT_FAILED", "dependency_run": run}, ensure_ascii=True))
    if not (P437_OUT.exists() and P437_SUMMARY.exists()):
        raise SystemExit(
            json.dumps(
                {
                    "status": "DEPENDENCY_OUTPUT_MISSING",
                    "missing": [
                        (str(P437_OUT.relative_to(REPO)) if not P437_OUT.exists() else None),
                        (str(P437_SUMMARY.relative_to(REPO)) if not P437_SUMMARY.exists() else None),
                    ],
                },
                ensure_ascii=True,
            )
        )

    r14 = load_json(R14_K_TOTAL)
    r15 = load_json(R15_M0)
    p437 = load_json(P437_OUT)
    p437_summary = load_json(P437_SUMMARY)

    k_total = extract_k_total_matrix(r14)

    m0_sq = r15.get("diagonal_decomposition", {}).get("host_diagonal_floor_value")
    if not isinstance(m0_sq, (int, float)) or not math.isfinite(float(m0_sq)):
        raise SystemExit(json.dumps({"status": "BAD_R15_M0_SQUARED", "path": str(R15_M0.relative_to(REPO))}, ensure_ascii=True))
    m0_sq_f = float(m0_sq)

    d = (p437.get("computed") or {}).get("d_local_residual_profile")
    if not (isinstance(d, list) and len(d) == 12 and all(isinstance(x, (int, float)) and math.isfinite(float(x)) for x in d)):
        raise SystemExit(
            json.dumps(
                {"status": "BAD_P437_DIAGONAL_PROFILE", "expected": "length-12 finite list", "path": str(P437_OUT.relative_to(REPO))},
                ensure_ascii=True,
            )
        )
    d_residual = np.array([float(x) for x in d], dtype=float)

    diag_full = m0_sq_f + d_residual
    h_psi = k_total.copy()
    h_psi[np.diag_indices_from(h_psi)] += diag_full

    # Symmetry audit
    sym_res = float(np.linalg.norm(h_psi - h_psi.T))

    evals, evecs = np.linalg.eigh(h_psi)
    evals_l = [float(x) for x in evals.tolist()]
    evecs_canon = np.zeros_like(evecs)
    for j in range(evecs.shape[1]):
        evecs_canon[:, j] = canonicalize_eigenvector_sign(evecs[:, j])

    orth_res = float(np.linalg.norm(evecs_canon.T @ evecs_canon - np.eye(12)))

    deg_tol = 1e-9
    degeneracy_groups = group_degeneracies(evals_l, tol=deg_tol)

    strict_inputs_marked = bool(p437_summary.get("input_marked_strict_derived") is True) and bool(
        p437_summary.get("provider_marked_strict_derived") is True
    )
    strict_provenance_chain = bool(p437_summary.get("provider_theorem_refs_present") is True) and bool(
        p437_summary.get("theorem_level_pass") is True
    )

    artifact = {
        "object": "PsiHessian_diagonal_local_strict_derived_value_instantiated_v1",
        "stage": "F459",
        "status": "actual_exported_strict_derived_value_instantiation__psi_hessian_eigensystem__no_false_pass",
        "as_of": "2026-03-15",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": "export_full_psi_sector_hessian_matrix_and_eigensystem_from_diagonal_local_strict_derived_lane",
        "inputs": {
            "r14_k_total": str(R14_K_TOTAL.relative_to(REPO)),
            "r15_m0_squared": str(R15_M0.relative_to(REPO)),
            "p437_diagonal_local_residual_profile": str(P437_OUT.relative_to(REPO)),
            "p437_summary": str(P437_SUMMARY.relative_to(REPO)),
        },
        "provenance_checks": {
            "p437_inputs_marked_strict_derived": strict_inputs_marked,
            "p437_theorem_chain_present": strict_provenance_chain,
            "note": "This packet is a strict-derived value instantiation on the diagonal/local lane; it does not claim host matching/cancellation.",
        },
        "construction": {
            "matrix_definition": "H_psi := K_total + diag(m0^2 + d_local_residual_profile)",
            "basis": "psi0..psi11 site basis (canonical 12-slot carrier)",
            "eigensystem": "numpy.linalg.eigh on the symmetric 12x12 matrix",
            "eigenvector_sign_convention": "flip v -> -v so that the component with maximal |v_i| is positive",
        },
        "outputs": {
            "n": 12,
            "m0_squared": m0_sq_f,
            "d_local_residual_profile": [float(x) for x in d_residual.tolist()],
            "diag_full_m0_plus_residual": [float(x) for x in diag_full.tolist()],
            "k_total": [[float(x) for x in row.tolist()] for row in k_total],
            "h_psi": [[float(x) for x in row.tolist()] for row in h_psi],
            "eigenvalues_ascending": evals_l,
            "eigenvectors_columns": [[float(x) for x in evecs_canon[:, j].tolist()] for j in range(12)],
        },
        "audits": {
            "symmetry_residual_l2": sym_res,
            "evecs_orthonormal_residual_l2": orth_res,
            "degeneracy_tolerance": deg_tol,
            "degeneracy_groups": degeneracy_groups,
        },
        "hard_limits": [
            "Does not claim host matching/cancellation for the diagonal/local residual sector.",
            "Does not claim global discharge of QW-2191 (kernel-alone obstruction remains; this is lane-scoped).",
            "Does not claim strict-core selector closure nor admissible S_sel_int.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F459",
        "status": "F459_EXECUTED_CURRENT_STRICT_DIAGONAL_LOCAL_PSI_HESSIAN_EIGENSYSTEM_VALUE_INSTANTIATION_PACKET_NO_FALSE_PASS",
        "generated_utc": artifact["generated_utc"],
        "outputs": {
            "n": 12,
            "symmetry_residual_l2": sym_res,
            "evecs_orthonormal_residual_l2": orth_res,
            "min_eigenvalue": float(min(evals_l)),
            "max_eigenvalue": float(max(evals_l)),
            "degeneracy_group_count": int(len(degeneracy_groups)),
            "max_multiplicity": int(max(g["multiplicity"] for g in degeneracy_groups)),
        },
        "inputs": artifact["inputs"],
        "dependency_run": run,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

