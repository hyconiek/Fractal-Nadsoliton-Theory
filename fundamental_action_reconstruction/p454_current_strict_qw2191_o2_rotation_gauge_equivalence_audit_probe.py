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

IN_ASSIGNMENT = GENERATED / "mode_index_assignment_canonical_local_diagonal_strict_derived_v1.json"
IN_QW2190 = REPO / "report_qw2190_kernel_mode_representation_emergence_gate.json"

OUT_JSON = GENERATED / "p454_current_strict_qw2191_o2_rotation_gauge_equivalence_audit_probe.json"
OUT_SUMMARY = GENERATED / "p454_current_strict_qw2191_o2_rotation_gauge_equivalence_audit_probe_summary.json"


def load_json_path(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def kernel_value(d: float, omega: float, phi: float, beta: float, eta: float) -> float:
    return float(np.cos(omega * d + phi) / (1.0 + beta * (d**eta)))


def build_ktotal_ring(n: int, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    m = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d = min(abs(i - j), n - abs(i - j))
            m[i, j] = kernel_value(float(d), omega, phi, beta, eta)
    return 0.5 * (m + m.T)


def pauli_generators() -> list[np.ndarray]:
    s1 = np.array([[0, 1], [1, 0]], dtype=complex)
    s2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    s3 = np.array([[1, 0], [0, -1]], dtype=complex)
    return [0.5 * s for s in [s1, s2, s3]]


def gell_mann_generators() -> list[np.ndarray]:
    l1 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex)
    l2 = np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex)
    l3 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex)
    l4 = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex)
    l5 = np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], dtype=complex)
    l6 = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex)
    l7 = np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=complex)
    l8 = np.array([[1 / np.sqrt(3), 0, 0], [0, 1 / np.sqrt(3), 0], [0, 0, -2 / np.sqrt(3)]], dtype=complex)
    return [0.5 * l for l in [l1, l2, l3, l4, l5, l6, l7, l8]]


def embed_generators(b: np.ndarray, gens: list[np.ndarray]) -> list[np.ndarray]:
    bc = b.astype(complex)
    return [bc @ g @ bc.T.conj() for g in gens]


def su2_closure_residual(ts: list[np.ndarray]) -> float:
    eps = np.zeros((3, 3, 3), dtype=float)
    eps[0, 1, 2] = eps[1, 2, 0] = eps[2, 0, 1] = 1.0
    eps[1, 0, 2] = eps[2, 1, 0] = eps[0, 2, 1] = -1.0
    r = 0.0
    for i in range(3):
        for j in range(3):
            lhs = ts[i] @ ts[j] - ts[j] @ ts[i]
            rhs = 1j * sum(eps[i, j, k] * ts[k] for k in range(3))
            r = max(r, float(np.linalg.norm(lhs - rhs)))
    return r


def su3_closure_residual(ts: list[np.ndarray]) -> float:
    r = 0.0
    for a in range(8):
        for b in range(8):
            comm = ts[a] @ ts[b] - ts[b] @ ts[a]
            coeff = [(-2j * np.trace(comm @ ts[c])) for c in range(8)]
            recon = 1j * sum(coeff[c] * ts[c] for c in range(8))
            r = max(r, float(np.linalg.norm(comm - recon)))
    return r


def inv_residual(k: np.ndarray, b: np.ndarray) -> float:
    p = b @ b.T
    n = k.shape[0]
    return float(np.linalg.norm((np.eye(n) - p) @ k @ p))


def orthonormal_residual(b: np.ndarray) -> float:
    gram = b.T @ b
    return float(np.linalg.norm(gram - np.eye(gram.shape[0])))


def lift_subspace_conjugator(b: np.ndarray, r: np.ndarray) -> np.ndarray:
    # For orthonormal columns b: P = b b^T is the orthogonal projector.
    # For orthogonal r: U = b r b^T + (I - P) acts as r on span(b) and identity on the complement.
    p = b @ b.T
    n = int(b.shape[0])
    return b @ r @ b.T + (np.eye(n) - p)


def max_conjugation_residual(u: np.ndarray, base: list[np.ndarray], alt: list[np.ndarray]) -> float:
    r = 0.0
    for g0, g1 in zip(base, alt):
        r = max(r, float(np.linalg.norm(g1 - (u @ g0 @ u.T.conj()))))
    return r


def rot2(theta: float) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([[c, -s], [s, c]], dtype=float)


def rot3_pair1(theta: float) -> np.ndarray:
    r = np.eye(3, dtype=float)
    r[1:, 1:] = rot2(theta)
    return r


def evaluate_case(
    case_id: str,
    ktotal: np.ndarray,
    b3: np.ndarray,
    b2: np.ndarray,
    r3: np.ndarray,
    r2: np.ndarray,
    su3_base: list[np.ndarray],
    su2_base: list[np.ndarray],
) -> dict[str, Any]:
    b3_alt = b3 @ r3
    b2_alt = b2 @ r2

    su3_alt = embed_generators(b3_alt, gell_mann_generators())
    su2_alt = embed_generators(b2_alt, pauli_generators())

    u3 = lift_subspace_conjugator(b3, r3)
    u2 = lift_subspace_conjugator(b2, r2)

    return {
        "case_id": case_id,
        "angles": {
            "theta_pair1": float(math.atan2(r3[2, 1], r3[1, 1])),
            "theta_pair2": float(math.atan2(r2[1, 0], r2[0, 0])),
        },
        "audits": {
            "b3_orthonormal_residual": orthonormal_residual(b3_alt),
            "b2_orthonormal_residual": orthonormal_residual(b2_alt),
            "inv_res_su3": inv_residual(ktotal, b3_alt),
            "inv_res_su2": inv_residual(ktotal, b2_alt),
            "su3_closure_res": su3_closure_residual(su3_alt),
            "su2_closure_res": su2_closure_residual(su2_alt),
            "u3_orthonormal_residual": float(np.linalg.norm(u3 @ u3.T - np.eye(u3.shape[0]))),
            "u2_orthonormal_residual": float(np.linalg.norm(u2 @ u2.T - np.eye(u2.shape[0]))),
            "su3_conjugation_residual": max_conjugation_residual(u3, su3_base, su3_alt),
            "su2_conjugation_residual": max_conjugation_residual(u2, su2_base, su2_alt),
        },
    }


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_ASSIGNMENT.exists():
        raise SystemExit(json.dumps({"status": "MISSING_F453_ASSIGNMENT_OBJECT", "missing": str(IN_ASSIGNMENT)}, ensure_ascii=True))

    assignment = load_json_path(IN_ASSIGNMENT)
    q2190 = load_json_path(IN_QW2190)

    n = int(q2190["mode_mapping"]["n_octaves"])
    kpars = {k: float(v) for k, v in q2190["kernel"].items()}
    ktotal = build_ktotal_ring(n, kpars["omega"], kpars["phi"], kpars["beta"], kpars["eta"])

    e0 = np.array([float(x) for x in assignment["outputs"]["e0"]], dtype=float)
    pair1 = assignment["outputs"]["pairs"]["pair1"]
    pair2 = assignment["outputs"]["pairs"]["pair2"]
    u1p = np.array([float(x) for x in pair1["u_plus"]], dtype=float)
    u1m = np.array([float(x) for x in pair1["u_minus"]], dtype=float)
    u2p = np.array([float(x) for x in pair2["u_plus"]], dtype=float)
    u2m = np.array([float(x) for x in pair2["u_minus"]], dtype=float)

    b3 = np.column_stack([e0, u1p, u1m])
    b2 = np.column_stack([u2p, u2m])

    base = {
        "b3_orthonormal_residual": orthonormal_residual(b3),
        "b2_orthonormal_residual": orthonormal_residual(b2),
        "inv_res_su3": inv_residual(ktotal, b3),
        "inv_res_su2": inv_residual(ktotal, b2),
    }

    su3_base = embed_generators(b3, gell_mann_generators())
    su2_base = embed_generators(b2, pauli_generators())

    base["su3_closure_res"] = su3_closure_residual(su3_base)
    base["su2_closure_res"] = su2_closure_residual(su2_base)

    angles = [
        0.0,
        math.pi / 12.0,
        math.pi / 8.0,
        math.pi / 7.0,
        math.pi / 6.0,
        math.pi / 5.0,
        math.pi / 4.0,
        math.pi / 3.0,
        math.pi / 2.0,
    ]

    cases: list[dict[str, Any]] = []
    for t1 in angles:
        for t2 in angles:
            r3 = rot3_pair1(t1)
            r2 = rot2(t2)
            case_id = f"pair1_theta={t1:.12g}__pair2_theta={t2:.12g}"
            cases.append(evaluate_case(case_id, ktotal, b3, b2, r3, r2, su3_base, su2_base))

    max_inv_su3_dev = float(max(abs(c["audits"]["inv_res_su3"] - base["inv_res_su3"]) for c in cases))
    max_inv_su2_dev = float(max(abs(c["audits"]["inv_res_su2"] - base["inv_res_su2"]) for c in cases))
    max_su3_cl_dev = float(max(abs(c["audits"]["su3_closure_res"] - base["su3_closure_res"]) for c in cases))
    max_su2_cl_dev = float(max(abs(c["audits"]["su2_closure_res"] - base["su2_closure_res"]) for c in cases))
    max_su3_conj = float(max(c["audits"]["su3_conjugation_residual"] for c in cases))
    max_su2_conj = float(max(c["audits"]["su2_conjugation_residual"] for c in cases))
    max_u3_orth = float(max(c["audits"]["u3_orthonormal_residual"] for c in cases))
    max_u2_orth = float(max(c["audits"]["u2_orthonormal_residual"] for c in cases))

    tol = 1e-12

    checks = [
        {
            "id": "baseline_bases_orthonormal",
            "pass": bool(base["b3_orthonormal_residual"] <= tol and base["b2_orthonormal_residual"] <= tol),
            "expected": f"<= {tol:g}",
            "actual": {
                "b3_orthonormal_residual": base["b3_orthonormal_residual"],
                "b2_orthonormal_residual": base["b2_orthonormal_residual"],
            },
            "meaning": "baseline embedded bases are numerically orthonormal (precondition for clean projector/lift arguments)",
        },
        {
            "id": "o2_rotation_family_conjugation_equivalence",
            "pass": bool(max_su3_conj <= tol and max_su2_conj <= tol),
            "expected": f"<= {tol:g}",
            "actual": {"max_su3_conjugation_residual": max_su3_conj, "max_su2_conjugation_residual": max_su2_conj},
            "meaning": "every tested O(2) rotation embedding is conjugate to the baseline embedding by a lifted subspace conjugator",
        },
        {
            "id": "o2_rotation_family_invariance_residual_unchanged",
            "pass": bool(max_inv_su3_dev <= tol and max_inv_su2_dev <= tol),
            "expected": f"<= {tol:g}",
            "actual": {"max_inv_res_su3_dev": max_inv_su3_dev, "max_inv_res_su2_dev": max_inv_su2_dev},
            "meaning": "kernel-subspace invariance residual is unchanged across tested O(2) rotations (as expected from projector invariance)",
        },
        {
            "id": "o2_rotation_family_lie_closure_residual_unchanged",
            "pass": bool(max_su3_cl_dev <= tol and max_su2_cl_dev <= tol),
            "expected": f"<= {tol:g}",
            "actual": {"max_su3_closure_res_dev": max_su3_cl_dev, "max_su2_closure_res_dev": max_su2_cl_dev},
            "meaning": "Lie-closure residual is unchanged across tested O(2) rotations (as expected from conjugation equivalence)",
        },
        {
            "id": "lifted_conjugators_orthonormal",
            "pass": bool(max_u3_orth <= tol and max_u2_orth <= tol),
            "expected": f"<= {tol:g}",
            "actual": {"max_u3_orthonormal_residual": max_u3_orth, "max_u2_orthonormal_residual": max_u2_orth},
            "meaning": "lifted conjugators are numerically orthonormal across tested O(2) rotations",
        },
    ]

    passed = bool(all(bool(c["pass"]) for c in checks))

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P454",
        "goal": (
            "audit_that_continuous_O2_basis_rotations_inside_pair1_pair2_planes_act_only_as_conjugations_for_QW2190_SU3_SU2_embeddings "
            "and_do_not_change_invariance_or_lie_closure_audits"
        ),
        "inputs": {
            "mode_index_assignment_object": str(IN_ASSIGNMENT.relative_to(REPO)),
            "qw2190_mode_scaffold": str(IN_QW2190.relative_to(REPO)),
            "angles_tested": angles,
        },
        "baseline": base,
        "family_max_deviations": {
            "max_inv_res_su3_dev": max_inv_su3_dev,
            "max_inv_res_su2_dev": max_inv_su2_dev,
            "max_su3_closure_res_dev": max_su3_cl_dev,
            "max_su2_closure_res_dev": max_su2_cl_dev,
            "max_su3_conjugation_residual": max_su3_conj,
            "max_su2_conjugation_residual": max_su2_conj,
            "max_u3_orthonormal_residual": max_u3_orth,
            "max_u2_orthonormal_residual": max_u2_orth,
            "tolerance": tol,
        },
        "cases": cases,
        "checks": checks,
        "verdict": {
            "audits_pass": passed,
            "statement": (
                "Across a grid of continuous O(2) rotations applied to the pair1/pair2 eigenbasis vectors exported by F453, "
                "the induced SU(3)/SU(2) embedded generators are conjugate to the baseline embedding by lifted subspace conjugators, "
                "and invariance/Lie-closure residuals remain unchanged at numerical tolerance. "
                "This supports treating the full O(2) basis-rotation freedom as a gauge/basis convention for the QW-2190 embedding audits "
                "(no strict promotion to global physical uniqueness)."
            ),
        },
        "hard_limits": [
            "Does not claim axiom-free global physical uniqueness beyond the declared QW-2190 embedding audits and the tested pair1/pair2 scope.",
            "Does not claim global discharge of QW-2191 (kernel-alone obstruction remains true as a canonical-representative obstruction).",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": report["stage"],
        "status": "PASS_O2_ROTATION_GAUGE_EQUIVALENCE_AUDIT" if passed else "FAIL_O2_ROTATION_GAUGE_EQUIVALENCE_AUDIT",
        "audits_pass": passed,
        "angles_tested_count": len(angles),
        "case_count": len(cases),
        "max_su3_conjugation_residual": max_su3_conj,
        "max_su2_conjugation_residual": max_su2_conj,
        "max_inv_res_su3_dev": max_inv_su3_dev,
        "max_inv_res_su2_dev": max_inv_su2_dev,
        "max_su3_closure_res_dev": max_su3_cl_dev,
        "max_su2_closure_res_dev": max_su2_cl_dev,
        "max_u3_orthonormal_residual": max_u3_orth,
        "max_u2_orthonormal_residual": max_u2_orth,
        "tolerance": tol,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()

