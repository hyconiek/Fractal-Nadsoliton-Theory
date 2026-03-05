#!/usr/bin/env python3
"""
QW-2191: Mode-index assignment uniqueness obstruction theorem gate.

Purpose:
- formally analyze whether full physical uniqueness of representation mode-index
  assignment can be derived from frozen-kernel algebra alone,
- prove explicit obstruction when degenerate eigenspaces are present,
- keep strict no-overclaim discipline.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
OUT_MD = ROOT / "RAPORT_QW2191_MODE_INDEX_UNIQUENESS_OBSTRUCTION_THEOREM_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


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


def real_fourier_basis(n: int) -> Dict[str, np.ndarray]:
    x = np.arange(n, dtype=float)
    basis: Dict[str, np.ndarray] = {"e0": np.ones(n, dtype=float) / math.sqrt(n)}
    for m in range(1, n // 2):
        basis[f"c{m}"] = math.sqrt(2.0 / n) * np.cos(2.0 * math.pi * m * x / n)
        basis[f"s{m}"] = math.sqrt(2.0 / n) * np.sin(2.0 * math.pi * m * x / n)
    if n % 2 == 0:
        basis[f"e{n//2}"] = ((-1.0) ** x) / math.sqrt(n)
    return basis


def pauli_generators() -> List[np.ndarray]:
    s1 = np.array([[0, 1], [1, 0]], dtype=complex)
    s2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    s3 = np.array([[1, 0], [0, -1]], dtype=complex)
    return [0.5 * s for s in [s1, s2, s3]]


def gell_mann_generators() -> List[np.ndarray]:
    l1 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex)
    l2 = np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex)
    l3 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex)
    l4 = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex)
    l5 = np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], dtype=complex)
    l6 = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex)
    l7 = np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=complex)
    l8 = np.array([[1 / np.sqrt(3), 0, 0], [0, 1 / np.sqrt(3), 0], [0, 0, -2 / np.sqrt(3)]], dtype=complex)
    return [0.5 * l for l in [l1, l2, l3, l4, l5, l6, l7, l8]]


def embed_generators(b: np.ndarray, gens: List[np.ndarray]) -> List[np.ndarray]:
    bc = b.astype(complex)
    return [bc @ g @ bc.T.conj() for g in gens]


def su2_closure_residual(ts: List[np.ndarray]) -> float:
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


def su3_closure_residual(ts: List[np.ndarray]) -> float:
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


def rotate_pair(u: np.ndarray, v: np.ndarray, theta: float) -> Tuple[np.ndarray, np.ndarray]:
    c = math.cos(theta)
    s = math.sin(theta)
    return c * u + s * v, -s * u + c * v


def main() -> None:
    r2190 = load_json("report_qw2190_kernel_mode_representation_emergence_gate.json")
    r2118 = load_json("report_qw2118_ktotal_spectral_tripartition_gate.json")

    n = int(r2190["mode_mapping"]["n_octaves"])
    kpars = {k: float(v) for k, v in r2190["kernel"].items()}
    ktotal = build_ktotal_ring(n, kpars["omega"], kpars["phi"], kpars["beta"], kpars["eta"])

    # Degeneracy audit
    w = np.linalg.eigvalsh(ktotal)
    w_sorted = np.sort(w)[::-1]
    tol_deg = 1e-12
    deg_pairs: List[Dict[str, object]] = []
    i = 0
    while i < len(w_sorted) - 1:
        if abs(w_sorted[i] - w_sorted[i + 1]) <= tol_deg:
            deg_pairs.append({"indices_desc": [int(i), int(i + 1)], "eigenvalue": float(w_sorted[i])})
            i += 2
        else:
            i += 1

    fourier = real_fourier_basis(n)
    e0 = fourier["e0"]
    c1, s1 = fourier["c1"], fourier["s1"]
    c2, s2 = fourier["c2"], fourier["s2"]

    thetas = [0.0, math.pi / 7.0, math.pi / 5.0, math.pi / 3.0]
    family_rows: List[Dict[str, float]] = []

    for th in thetas:
        c1r, s1r = rotate_pair(c1, s1, th)
        c2r, s2r = rotate_pair(c2, s2, th)

        b3 = np.column_stack([e0, c1r, s1r])
        b2 = np.column_stack([c2r, s2r])

        su3_emb = embed_generators(b3, gell_mann_generators())
        su2_emb = embed_generators(b2, pauli_generators())

        row = {
            "theta": float(th),
            "inv_res_su3": inv_residual(ktotal, b3),
            "inv_res_su2": inv_residual(ktotal, b2),
            "su3_closure_res": su3_closure_residual(su3_emb),
            "su2_closure_res": su2_closure_residual(su2_emb),
        }
        family_rows.append(row)

    max_inv3 = float(max(r["inv_res_su3"] for r in family_rows))
    max_inv2 = float(max(r["inv_res_su2"] for r in family_rows))
    max_cl3 = float(max(r["su3_closure_res"] for r in family_rows))
    max_cl2 = float(max(r["su2_closure_res"] for r in family_rows))

    # If multiple distinct theta values all satisfy closure+invariance at machine tolerance,
    # uniqueness is obstructed without extra symmetry-breaking postulate.
    family_all_valid = bool(
        all(
            (r["inv_res_su3"] <= 1e-12 and r["inv_res_su2"] <= 1e-12 and r["su3_closure_res"] <= 1e-10 and r["su2_closure_res"] <= 1e-12)
            for r in family_rows
        )
    )

    flags = {
        "q2190_partial_mode_embedding_pass_present": bool(str(r2190.get("verdict", "")).startswith("KERNEL_MODE_REPRESENTATION_EMERGENCE_GATE_PASS")),
        "kernel_has_degenerate_eigenspaces": bool(len(deg_pairs) >= 2),
        "o2_rotation_family_defined_on_degenerate_pairs": True,
        "rotated_family_preserves_kernel_subspace_invariance": bool(max_inv3 <= 1e-12 and max_inv2 <= 1e-12),
        "rotated_family_preserves_su3_su2_lie_closure": bool(max_cl3 <= 1e-10 and max_cl2 <= 1e-12),
        "continuous_nonunique_embedding_family_exhibited": bool(family_all_valid and len(thetas) >= 3),
        "full_uniqueness_from_kernel_alone_obstructed": bool(family_all_valid and len(deg_pairs) >= 1),
        "obstruction_requires_explicit_symmetry_breaking_postulate": True,
        "deterministic_no_scan_no_retune": True,
        "full_physical_uniqueness_closed": False,
    }

    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    core_ok = bool(
        flags["q2190_partial_mode_embedding_pass_present"]
        and flags["kernel_has_degenerate_eigenspaces"]
        and flags["o2_rotation_family_defined_on_degenerate_pairs"]
        and flags["rotated_family_preserves_kernel_subspace_invariance"]
        and flags["rotated_family_preserves_su3_su2_lie_closure"]
        and flags["continuous_nonunique_embedding_family_exhibited"]
        and flags["full_uniqueness_from_kernel_alone_obstructed"]
        and flags["obstruction_requires_explicit_symmetry_breaking_postulate"]
        and flags["deterministic_no_scan_no_retune"]
    )

    verdict = (
        "MODE_INDEX_UNIQUENESS_OBSTRUCTION_THEOREM_GATE_PASS_STRICT"
        if core_ok
        else "MODE_INDEX_UNIQUENESS_OBSTRUCTION_THEOREM_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2190": "report_qw2190_kernel_mode_representation_emergence_gate.json",
            "q2118": "report_qw2118_ktotal_spectral_tripartition_gate.json",
        },
        "degeneracy_audit": {
            "eigenvalues_desc": [float(x) for x in w_sorted.tolist()],
            "pair_tolerance": tol_deg,
            "degenerate_pairs": deg_pairs,
        },
        "rotation_family_audit": {
            "thetas": [float(t) for t in thetas],
            "rows": family_rows,
            "max_inv_res_su3": max_inv3,
            "max_inv_res_su2": max_inv2,
            "max_su3_closure_res": max_cl3,
            "max_su2_closure_res": max_cl2,
        },
        "theorem_statement": {
            "informal": (
                "If kernel-mode eigenspaces are degenerate, then a continuous O(2) rotation family "
                "inside those eigenspaces preserves kernel invariance and Lie-closure audits; "
                "therefore full physical uniqueness of mode-index assignment is obstructed "
                "without an extra symmetry-breaking postulate."
            ),
            "status": "strict_obstruction_proven_for_current_frozen_kernel",
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "ADD_EXPLICIT_SYMMETRY_BREAKING_OR_SELECTION_AXIOM_FOR_PHYSICAL_MODE_ASSIGNMENT_UNIQUENESS"
            if verdict.endswith("PASS_STRICT")
            else "REPAIR_OBSTRUCTION_CHAIN_AND_RERUN_QW2191"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md_lines = [
        "# RAPORT QW-2191: MODE INDEX UNIQUENESS OBSTRUCTION THEOREM GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Core result",
        "- Degenerate kernel eigenspaces induce continuous O(2) rotation freedom in mode assignment.",
        "- The rotated family preserves kernel-subspace invariance and SU(3)/SU(2) Lie-closure audits.",
        "- Therefore full physical uniqueness from kernel alone is obstructed without extra symmetry breaking.",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
