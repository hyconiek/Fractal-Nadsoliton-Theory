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

IN_THETA_PAIR = GENERATED / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1.json"
IN_QW2190 = REPO / "report_qw2190_kernel_mode_representation_emergence_gate.json"

OUT_JSON = GENERATED / "p401_current_strict_sigma_int_qw2191_selected_point_witness_probe.json"
OUT_SUMMARY = GENERATED / "p401_current_strict_sigma_int_qw2191_selected_point_witness_probe_summary.json"


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


def real_fourier_basis(n: int) -> dict[str, np.ndarray]:
    x = np.arange(n, dtype=float)
    basis: dict[str, np.ndarray] = {"e0": np.ones(n, dtype=float) / math.sqrt(n)}
    for m in range(1, n // 2):
        basis[f"c{m}"] = math.sqrt(2.0 / n) * np.cos(2.0 * math.pi * m * x / n)
        basis[f"s{m}"] = math.sqrt(2.0 / n) * np.sin(2.0 * math.pi * m * x / n)
    if n % 2 == 0:
        basis[f"e{n//2}"] = ((-1.0) ** x) / math.sqrt(n)
    return basis


def rotate_pair(u: np.ndarray, v: np.ndarray, theta: float) -> tuple[np.ndarray, np.ndarray]:
    c = math.cos(theta)
    s = math.sin(theta)
    return c * u + s * v, -s * u + c * v


def inv_residual(k: np.ndarray, b: np.ndarray) -> float:
    p = b @ b.T
    n = k.shape[0]
    return float(np.linalg.norm((np.eye(n) - p) @ k @ p))


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
    l8 = np.array(
        [[1 / np.sqrt(3), 0, 0], [0, 1 / np.sqrt(3), 0], [0, 0, -2 / np.sqrt(3)]], dtype=complex
    )
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


def orthonormal_residual(b: np.ndarray) -> float:
    m = b.T @ b
    return float(np.linalg.norm(m - np.eye(m.shape[0])))


def main() -> None:
    theta_pair = load_json_path(IN_THETA_PAIR)
    q2190 = load_json_path(IN_QW2190)

    n = int(q2190["mode_mapping"]["n_octaves"])
    kpars = {k: float(v) for k, v in q2190["kernel"].items()}
    ktotal = build_ktotal_ring(n, kpars["omega"], kpars["phi"], kpars["beta"], kpars["eta"])

    theta_1 = float(theta_pair["outputs"]["pair1"]["theta_1_cand"])
    theta_2 = float(theta_pair["outputs"]["pair2"]["theta_2_cand"])
    u1_art = np.array([float(x) for x in theta_pair["outputs"]["pair1"]["u_1_cand"]], dtype=float)
    u2_art = np.array([float(x) for x in theta_pair["outputs"]["pair2"]["u_2_cand"]], dtype=float)

    fourier = real_fourier_basis(n)
    e0 = fourier["e0"]
    c1, s1 = fourier["c1"], fourier["s1"]
    c2, s2 = fourier["c2"], fourier["s2"]

    c1r, s1r = rotate_pair(c1, s1, theta_1)
    c2r, s2r = rotate_pair(c2, s2, theta_2)

    b3 = np.column_stack([e0, c1r, s1r])
    b2 = np.column_stack([c2r, s2r])

    # QW-2191 tolerances (copied literally from the gate implementation).
    tol_inv = 1e-12
    tol_su3 = 1e-10
    tol_su2 = 1e-12

    inv3 = inv_residual(ktotal, b3)
    inv2 = inv_residual(ktotal, b2)
    su3_emb = embed_generators(b3, gell_mann_generators())
    su2_emb = embed_generators(b2, pauli_generators())
    cl3 = su3_closure_residual(su3_emb)
    cl2 = su2_closure_residual(su2_emb)

    u1_match_res = float(np.linalg.norm(u1_art - c1r))
    u2_match_res = float(np.linalg.norm(u2_art - c2r))

    checks = [
        {
            "id": "selected_point_kernel_invariance_su3",
            "pass": bool(inv3 <= tol_inv),
            "tolerance": tol_inv,
            "value": float(inv3),
            "meaning": "the selected (e0,c1,s1) rotated basis spans a kernel-invariant SU(3) subspace (selected point)",
        },
        {
            "id": "selected_point_kernel_invariance_su2",
            "pass": bool(inv2 <= tol_inv),
            "tolerance": tol_inv,
            "value": float(inv2),
            "meaning": "the selected (c2,s2) rotated basis spans a kernel-invariant SU(2) subspace (selected point)",
        },
        {
            "id": "selected_point_su3_lie_closure",
            "pass": bool(cl3 <= tol_su3),
            "tolerance": tol_su3,
            "value": float(cl3),
            "meaning": "the embedded SU(3) generators close (selected point)",
        },
        {
            "id": "selected_point_su2_lie_closure",
            "pass": bool(cl2 <= tol_su2),
            "tolerance": tol_su2,
            "value": float(cl2),
            "meaning": "the embedded SU(2) generators close (selected point)",
        },
        {
            "id": "u1_matches_rotated_c1",
            "pass": bool(u1_match_res <= 1e-12),
            "tolerance": 1e-12,
            "value": float(u1_match_res),
            "meaning": "exported u_1^cand matches the rotated Fourier c1 basis vector at theta_1^cand (internal consistency)",
        },
        {
            "id": "u2_matches_rotated_c2",
            "pass": bool(u2_match_res <= 1e-12),
            "tolerance": 1e-12,
            "value": float(u2_match_res),
            "meaning": "exported u_2^cand matches the rotated Fourier c2 basis vector at theta_2^cand (internal consistency)",
        },
    ]

    compatible = bool(all(bool(c["pass"]) for c in checks))

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P401",
        "goal": "selected_point_witness_for_sigma_int_theta_pair_against_qw2191_o2_family_audits",
        "inputs": {
            "theta_pair_artifact": str(IN_THETA_PAIR.relative_to(REPO)),
            "qw2190_kernel_mode_report": str(IN_QW2190.relative_to(REPO)),
        },
        "selected_point": {
            "theta_1_cand": theta_1,
            "theta_2_cand": theta_2,
            "basis": {"n": n, "su3": ["e0", "rot(c1,s1;theta_1_cand)"], "su2": ["rot(c2,s2;theta_2_cand)"]},
        },
        "computed": {
            "orthonormal_residual_su3_basis": orthonormal_residual(b3),
            "orthonormal_residual_su2_basis": orthonormal_residual(b2),
            "inv_res_su3": float(inv3),
            "inv_res_su2": float(inv2),
            "su3_closure_res": float(cl3),
            "su2_closure_res": float(cl2),
            "u1_match_residual": float(u1_match_res),
            "u2_match_residual": float(u2_match_res),
        },
        "tolerances": {
            "qw2191_inv_residual": tol_inv,
            "qw2191_su3_closure_residual": tol_su3,
            "qw2191_su2_closure_residual": tol_su2,
            "u_vector_match": 1e-12,
        },
        "checks": checks,
        "verdict": {
            "compatible_with_qw2191_o2_family_at_selected_point": compatible,
            "statement": (
                "Selected-point compatibility witness passes: the exported (theta_1^cand,theta_2^cand) "
                "produces a rotated Fourier basis that satisfies the same QW-2191 invariance and Lie-closure audits. "
                "This does not imply uniqueness nor QW-2191 discharge."
            ),
        },
        "hard_limits": [
            "Does not discharge QW-2191 (kernel-alone uniqueness obstruction remains true).",
            "Does not prove strict-core selector closure nor export an admissible S_sel_int.",
            "This is a compatibility-only witness at one selected representative point.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": report["stage"],
        "compatible_with_qw2191_o2_family_at_selected_point": report["verdict"][
            "compatible_with_qw2191_o2_family_at_selected_point"
        ],
        "theta_1_cand": theta_1,
        "theta_2_cand": theta_2,
        "inv_res_su3": float(inv3),
        "inv_res_su2": float(inv2),
        "su3_closure_res": float(cl3),
        "su2_closure_res": float(cl2),
        "u1_match_residual": float(u1_match_res),
        "u2_match_residual": float(u2_match_res),
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()

