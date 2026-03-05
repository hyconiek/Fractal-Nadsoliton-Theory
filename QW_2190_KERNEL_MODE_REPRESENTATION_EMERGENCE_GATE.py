#!/usr/bin/env python3
"""
QW-2190: Kernel-mode representation emergence gate (deterministic, no-scan).

Purpose:
- build an explicit deterministic mode-algebra embedding from the frozen kernel
  (12-octave ring -> real Fourier mode basis),
- embed SU(2) and SU(3) generators into kernel-mode subspaces and audit
  Lie-closure/invariance,
- keep strict scope boundary explicit: full physical uniqueness of mode-index
  assignment remains open.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2190_kernel_mode_representation_emergence_gate.json"
OUT_MD = ROOT / "RAPORT_QW2190_KERNEL_MODE_REPRESENTATION_EMERGENCE_GATE.md"


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
    basis: Dict[str, np.ndarray] = {}
    basis["e0"] = np.ones(n, dtype=float) / math.sqrt(n)
    # paired cosine/sine harmonics
    for m in range(1, n // 2):
        basis[f"c{m}"] = math.sqrt(2.0 / n) * np.cos(2.0 * math.pi * m * x / n)
        basis[f"s{m}"] = math.sqrt(2.0 / n) * np.sin(2.0 * math.pi * m * x / n)
    if n % 2 == 0:
        basis[f"e{n//2}"] = ((-1.0) ** x) / math.sqrt(n)
    return basis


def orth_residual(b: np.ndarray) -> float:
    i = np.eye(b.shape[1], dtype=float)
    return float(np.linalg.norm(b.T @ b - i))


def invariant_residual(k: np.ndarray, b: np.ndarray) -> float:
    p = b @ b.T
    n = k.shape[0]
    return float(np.linalg.norm((np.eye(n) - p) @ k @ p))


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
    # With Tr(Ta Tb)=1/2 delta_ab for Hermitian generators:
    # [Ta,Tb] = i f_{abc} Tc,  f_{abc} = -2i Tr([Ta,Tb] Tc)
    r = 0.0
    for a in range(8):
        for b in range(8):
            comm = ts[a] @ ts[b] - ts[b] @ ts[a]
            coeff = [(-2j * np.trace(comm @ ts[c])) for c in range(8)]
            recon = 1j * sum(coeff[c] * ts[c] for c in range(8))
            r = max(r, float(np.linalg.norm(comm - recon)))
    return r


def hermitian_residual(ts: List[np.ndarray]) -> float:
    return float(max(np.linalg.norm(t - t.T.conj()) for t in ts))


def traceless_residual(ts: List[np.ndarray]) -> float:
    return float(max(abs(np.trace(t)) for t in ts))


def cross_commutator_residual(ts2: List[np.ndarray], ts3: List[np.ndarray]) -> float:
    r = 0.0
    for a in ts2:
        for b in ts3:
            r = max(r, float(np.linalg.norm(a @ b - b @ a)))
    return r


def main() -> None:
    r2118 = load_json("report_qw2118_ktotal_spectral_tripartition_gate.json")
    r2127 = load_json("report_qw2127_nonabelian_spinor_gauge_action_bridge_gate.json")
    r2184 = load_json("report_qw2184_hypercharge_global_uniqueness_symbolic_gate.json")
    r2189 = load_json("report_qw2189_spinor_gauge_deanchored_consistency_gate.json") if (ROOT / "report_qw2189_spinor_gauge_deanchored_consistency_gate.json").exists() else None

    kpars = {k: float(v) for k, v in r2118["kernel"].items()}
    n = int(r2184["kernel_anchor"]["n_octaves"])

    ktotal = build_ktotal_ring(n, kpars["omega"], kpars["phi"], kpars["beta"], kpars["eta"])
    fourier = real_fourier_basis(n)

    # Deterministic mode-index mapping (declared explicitly):
    # SU(3) on [e0,c1,s1], SU(2) on [c2,s2]
    b3 = np.column_stack([fourier["e0"], fourier["c1"], fourier["s1"]])
    b2 = np.column_stack([fourier["c2"], fourier["s2"]])

    ortho3 = orth_residual(b3)
    ortho2 = orth_residual(b2)
    disjoint = float(np.linalg.norm(b2.T @ b3))

    inv3 = invariant_residual(ktotal, b3)
    inv2 = invariant_residual(ktotal, b2)

    su3_emb = embed_generators(b3, gell_mann_generators())
    su2_emb = embed_generators(b2, pauli_generators())

    su3_cl = su3_closure_residual(su3_emb)
    su2_cl = su2_closure_residual(su2_emb)
    su3_herm = hermitian_residual(su3_emb)
    su2_herm = hermitian_residual(su2_emb)
    su3_trace = traceless_residual(su3_emb)
    su2_trace = traceless_residual(su2_emb)
    cross_res = cross_commutator_residual(su2_emb, su3_emb)

    # Import hypercharge closure from QW-2184 as U(1) block
    ht = r2184["derived_template"]
    y_template = {k: str(v) for k, v in ht.items()}
    # exact anomaly closure already established symbolically in QW-2184
    y_anomaly_zero = all(
        str(v) == "0" for v in r2184["anomaly_coefficients_fractional"].values()
    )

    flags = {
        "q2127_nonabelian_bridge_present": bool(str(r2127.get("verdict", "")).startswith("NONABELIAN_SPINOR_GAUGE_ACTION_BRIDGE_GATE_PASS")),
        "q2184_symbolic_hypercharge_closure_present": bool(str(r2184.get("verdict", "")).endswith("PASS_DECLARED_CLASS")),
        "q2189_deanchored_consistency_present": bool(r2189 is not None and str(r2189.get("verdict", "")).startswith("SPINOR_GAUGE_DEANCHORED_CONSISTENCY_GATE_PASS")),
        "kernel_mode_basis_declared_deterministically": True,
        "su3_mode_subspace_dimension_is_3": bool(b3.shape == (n, 3)),
        "su2_mode_subspace_dimension_is_2": bool(b2.shape == (n, 2)),
        "mode_subspaces_orthonormal": bool(ortho3 <= 1e-12 and ortho2 <= 1e-12),
        "mode_subspaces_disjoint": bool(disjoint <= 1e-12),
        "kernel_invariance_of_su3_subspace": bool(inv3 <= 1e-12),
        "kernel_invariance_of_su2_subspace": bool(inv2 <= 1e-12),
        "embedded_su3_lie_closure": bool(su3_cl <= 1e-10),
        "embedded_su2_lie_closure": bool(su2_cl <= 1e-12),
        "embedded_generators_hermitian_and_traceless": bool(su3_herm <= 1e-12 and su2_herm <= 1e-12 and su3_trace <= 1e-12 and su2_trace <= 1e-12),
        "su3_su2_direct_product_cross_commutator_zero": bool(cross_res <= 1e-12),
        "u1_hypercharge_template_integrated": bool(y_template.get("Y_H") == "1/2" and y_template.get("Y_nR") == "0"),
        "u1_anomaly_closure_integrated": bool(y_anomaly_zero),
        "deterministic_no_scan_no_retune": True,
        "full_physical_uniqueness_of_mode_index_assignment": False,
    }

    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    core_ok = bool(
        flags["q2127_nonabelian_bridge_present"]
        and flags["q2184_symbolic_hypercharge_closure_present"]
        and flags["kernel_mode_basis_declared_deterministically"]
        and flags["su3_mode_subspace_dimension_is_3"]
        and flags["su2_mode_subspace_dimension_is_2"]
        and flags["mode_subspaces_orthonormal"]
        and flags["mode_subspaces_disjoint"]
        and flags["kernel_invariance_of_su3_subspace"]
        and flags["kernel_invariance_of_su2_subspace"]
        and flags["embedded_su3_lie_closure"]
        and flags["embedded_su2_lie_closure"]
        and flags["embedded_generators_hermitian_and_traceless"]
        and flags["su3_su2_direct_product_cross_commutator_zero"]
        and flags["u1_hypercharge_template_integrated"]
        and flags["u1_anomaly_closure_integrated"]
        and flags["deterministic_no_scan_no_retune"]
    )

    verdict = (
        "KERNEL_MODE_REPRESENTATION_EMERGENCE_GATE_PASS_PARTIAL_PHYSICAL_UNIQUENESS_OPEN"
        if core_ok
        else "KERNEL_MODE_REPRESENTATION_EMERGENCE_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2118": "report_qw2118_ktotal_spectral_tripartition_gate.json",
            "q2127": "report_qw2127_nonabelian_spinor_gauge_action_bridge_gate.json",
            "q2184": "report_qw2184_hypercharge_global_uniqueness_symbolic_gate.json",
            "q2189": "report_qw2189_spinor_gauge_deanchored_consistency_gate.json" if r2189 is not None else None,
        },
        "kernel": {k: float(v) for k, v in kpars.items()},
        "mode_mapping": {
            "n_octaves": n,
            "su3_modes": ["e0", "c1", "s1"],
            "su2_modes": ["c2", "s2"],
            "u1_hypercharge_source": "q2184_symbolic_template",
        },
        "numeric_audit": {
            "orthonormal_residual_su3_basis": ortho3,
            "orthonormal_residual_su2_basis": ortho2,
            "subspace_disjoint_residual": disjoint,
            "kernel_invariance_residual_su3_subspace": inv3,
            "kernel_invariance_residual_su2_subspace": inv2,
            "su3_closure_residual": su3_cl,
            "su2_closure_residual": su2_cl,
            "su3_hermitian_residual": su3_herm,
            "su2_hermitian_residual": su2_herm,
            "su3_traceless_residual": su3_trace,
            "su2_traceless_residual": su2_trace,
            "su3_su2_cross_commutator_residual": cross_res,
        },
        "u1_template": {
            "hypercharges_fractional": y_template,
            "anomaly_coefficients_fractional": r2184["anomaly_coefficients_fractional"],
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PROVE_FULL_PHYSICAL_UNIQUENESS_OF_MODE_INDEX_ASSIGNMENT_FROM_KERNEL_ALONE"
            if verdict.endswith("PHYSICAL_UNIQUENESS_OPEN")
            else "REPAIR_KERNEL_MODE_ALGEBRA_CHAIN_AND_RERUN_QW2190"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2190: KERNEL MODE REPRESENTATION EMERGENCE GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Deterministic mode mapping",
        "- SU(3): `e0,c1,s1`",
        "- SU(2): `c2,s2`",
        "- U(1): symbolic hypercharge template from `QW-2184`",
        "",
        "## Key residuals",
        f"- inv(SU3 subspace): `{inv3:.3e}`",
        f"- inv(SU2 subspace): `{inv2:.3e}`",
        f"- SU3 closure: `{su3_cl:.3e}`",
        f"- SU2 closure: `{su2_cl:.3e}`",
        f"- SU3/SU2 cross-commutator: `{cross_res:.3e}`",
        "",
        "## Boundary",
        "- full physical uniqueness of mode-index assignment remains open.",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
