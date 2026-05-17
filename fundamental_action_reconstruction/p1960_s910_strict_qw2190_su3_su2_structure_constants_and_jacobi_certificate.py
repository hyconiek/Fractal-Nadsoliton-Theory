#!/usr/bin/env python3
"""P1960 S910 strict QW-2190 SU(3)/SU(2) structure constants.

P1959 correctly refused to promote the local Abelian ghost seed to a full
non-Abelian BRST theorem, but it checked only the FAR-local path for the
QW-2190 source report.  The current QW-2190 report used by P452/P454 is stored
at the repository root.  This executor ingests that report and exports the
smallest theorem-grade algebra data that can now be honestly exported:

* exact SU(3), SU(2), and U(1) structure constants in the QW-2190 convention,
* exact Jacobi certificates for the exported constants,
* a numerical embedding audit on the strict QW-2190 mode-assignment carrier.

It does not export a full BV/BRST operator map, a global Q^2 certificate, a
ghost-cancelled Cutkosky trace, or a global selector closure.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p1960_s910_strict_qw2190_su3_su2_structure_constants_and_jacobi_certificate.json"
QW2190_ROOT = REPO / "report_qw2190_kernel_mode_representation_emergence_gate.json"
QW2190_FAR_LOCAL = ROOT / "report_qw2190_kernel_mode_representation_emergence_gate.json"
ASSIGNMENT = GEN / "mode_index_assignment_canonical_local_diagonal_strict_derived_v1.json"
P1959 = GEN / "p1959_s909_strict_nonabelian_gauge_group_data_availability_audit.json"
P452_SUMMARY = GEN / "p452_current_strict_qw2191_residual_z2_sign_flip_gauge_equivalence_audit_probe_summary.json"
P454_SUMMARY = GEN / "p454_current_strict_qw2191_o2_rotation_gauge_equivalence_audit_probe_summary.json"

TEXT_SUFFIXES = {".py", ".md", ".json", ".txt", ".yaml", ".yml"}
SKIP_NAMES = {
    "TOE_FINAL_DOCUMENTATION_RELEASE_8_STRICT_FULL.pdf",
    Path(__file__).name,
    OUT.name,
}


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def digest_obj(obj: object) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def digest_path(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(REPO))
    except ValueError:
        return str(path)


def repo_text_files() -> list[Path]:
    paths: list[Path] = []
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        if path.name in SKIP_NAMES:
            continue
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        paths.append(path)
    return sorted(paths)


def scan_terms(terms: list[str]) -> dict[str, dict[str, Any]]:
    files = repo_text_files()
    result: dict[str, dict[str, Any]] = {}
    for term in terms:
        count = 0
        sample_paths: list[str] = []
        low = term.lower()
        for path in files:
            text = path.read_text(encoding="utf-8", errors="ignore").lower()
            hits = text.count(low)
            if hits:
                count += hits
                if len(sample_paths) < 8:
                    sample_paths.append(rel(path))
        result[term] = {"count": count, "sample_paths": sample_paths}
    return result


def gell_mann_generators_exact() -> list[sp.Matrix]:
    i = sp.I
    rt3 = sp.sqrt(3)
    lam = [
        sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]),
        sp.Matrix([[0, -i, 0], [i, 0, 0], [0, 0, 0]]),
        sp.Matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]]),
        sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
        sp.Matrix([[0, 0, -i], [0, 0, 0], [i, 0, 0]]),
        sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),
        sp.Matrix([[0, 0, 0], [0, 0, -i], [0, i, 0]]),
        sp.Matrix([[1 / rt3, 0, 0], [0, 1 / rt3, 0], [0, 0, -2 / rt3]]),
    ]
    return [m / 2 for m in lam]


def pauli_generators_exact() -> list[sp.Matrix]:
    i = sp.I
    sig = [
        sp.Matrix([[0, 1], [1, 0]]),
        sp.Matrix([[0, -i], [i, 0]]),
        sp.Matrix([[1, 0], [0, -1]]),
    ]
    return [m / 2 for m in sig]


def comm(a: sp.Matrix, b: sp.Matrix) -> sp.Matrix:
    return a * b - b * a


def is_zero_matrix(m: sp.Matrix) -> bool:
    return all(sp.simplify(x) == 0 for x in list(m))


def structure_constants(gens: list[sp.Matrix]) -> list[list[list[sp.Expr]]]:
    n = len(gens)
    out: list[list[list[sp.Expr]]] = []
    for a in range(n):
        layer_b: list[list[sp.Expr]] = []
        for b in range(n):
            layer_c: list[sp.Expr] = []
            cab = comm(gens[a], gens[b])
            for c in range(n):
                # With Tr(T_a T_b) = delta_ab/2 and [T_a,T_b]=i f_abc T_c.
                value = sp.simplify(-2 * sp.I * sp.trace(cab * gens[c]))
                layer_c.append(value)
            layer_b.append(layer_c)
        out.append(layer_b)
    return out


def sparse_rows(f: list[list[list[sp.Expr]]], factor: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for a, b, c in product(range(len(f)), repeat=3):
        value = sp.simplify(f[a][b][c])
        if value != 0:
            rows.append(
                {
                    "factor": factor,
                    "a": str(a + 1),
                    "b": str(b + 1),
                    "c": str(c + 1),
                    "value": str(value),
                }
            )
    return rows


def generator_certificate(gens: list[sp.Matrix], f: list[list[list[sp.Expr]]]) -> dict[str, Any]:
    n = len(gens)
    normalization_bad: list[str] = []
    hermitian_bad: list[str] = []
    traceless_bad: list[str] = []
    closure_bad: list[str] = []
    antisym_bad: list[str] = []

    for a, b in product(range(n), repeat=2):
        expected = sp.Integer(1) if a == b else sp.Integer(0)
        residual = sp.simplify(2 * sp.trace(gens[a] * gens[b]) - expected)
        if residual != 0:
            normalization_bad.append(f"{a + 1},{b + 1}:{residual}")

    for a in range(n):
        if not is_zero_matrix(gens[a].conjugate().transpose() - gens[a]):
            hermitian_bad.append(str(a + 1))
        if sp.simplify(sp.trace(gens[a])) != 0:
            traceless_bad.append(str(a + 1))

    for a, b in product(range(n), repeat=2):
        reconstructed = sp.zeros(gens[0].rows, gens[0].cols)
        for c in range(n):
            reconstructed += sp.I * f[a][b][c] * gens[c]
        if not is_zero_matrix(comm(gens[a], gens[b]) - reconstructed):
            closure_bad.append(f"{a + 1},{b + 1}")

    for a, b, c in product(range(n), repeat=3):
        if sp.simplify(f[a][b][c] + f[b][a][c]) != 0:
            antisym_bad.append(f"ab:{a + 1},{b + 1},{c + 1}")
        if sp.simplify(f[a][b][c] + f[a][c][b]) != 0:
            antisym_bad.append(f"bc:{a + 1},{b + 1},{c + 1}")

    return {
        "normalization_2Tr_TaTb_delta": len(normalization_bad) == 0,
        "hermitian_generators": len(hermitian_bad) == 0,
        "traceless_generators": len(traceless_bad) == 0,
        "closure_reconstruction_exact_zero": len(closure_bad) == 0,
        "f_abc_total_antisymmetry_exact_zero": len(antisym_bad) == 0,
        "bad_components": {
            "normalization": normalization_bad[:12],
            "hermitian": hermitian_bad[:12],
            "traceless": traceless_bad[:12],
            "closure": closure_bad[:12],
            "antisymmetry": antisym_bad[:12],
        },
    }


def jacobi_certificate(f: list[list[list[sp.Expr]]]) -> dict[str, Any]:
    n = len(f)
    bad: list[dict[str, str]] = []
    for a, b, c, d in product(range(n), repeat=4):
        value = sp.Integer(0)
        for e in range(n):
            value += (
                f[a][b][e] * f[e][c][d]
                + f[b][c][e] * f[e][a][d]
                + f[c][a][e] * f[e][b][d]
            )
        value = sp.simplify(value)
        if value != 0:
            bad.append(
                {
                    "a": str(a + 1),
                    "b": str(b + 1),
                    "c": str(c + 1),
                    "d": str(d + 1),
                    "value": str(value),
                }
            )
    return {
        "checked_components": n**4,
        "nonzero_components": len(bad),
        "all_zero": len(bad) == 0,
        "bad_components_sample": bad[:12],
    }


def kernel_value(d: float, omega: float, phi: float, beta: float, eta: float) -> float:
    return float(math.cos(omega * d + phi) / (1.0 + beta * (d**eta)))


def build_ktotal_ring(n: int, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    m = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d = min(abs(i - j), n - abs(i - j))
            m[i, j] = kernel_value(float(d), omega, phi, beta, eta)
    return 0.5 * (m + m.T)


def exact_mats_to_numpy(gens: list[sp.Matrix]) -> list[np.ndarray]:
    arrays: list[np.ndarray] = []
    for gen in gens:
        arrays.append(np.array([[complex(sp.N(x, 30)) for x in row] for row in gen.tolist()], dtype=complex))
    return arrays


def embed_generators(basis: np.ndarray, gens: list[np.ndarray]) -> list[np.ndarray]:
    bc = basis.astype(complex)
    return [bc @ g @ bc.T.conj() for g in gens]


def inv_residual(k: np.ndarray, basis: np.ndarray) -> float:
    p = basis @ basis.T
    n = k.shape[0]
    return float(np.linalg.norm((np.eye(n) - p) @ k @ p))


def orthonormal_residual(basis: np.ndarray) -> float:
    gram = basis.T @ basis
    return float(np.linalg.norm(gram - np.eye(gram.shape[0])))


def closure_residual(gens: list[np.ndarray], f: list[list[list[sp.Expr]]]) -> float:
    r = 0.0
    n = len(gens)
    f_np = np.zeros((n, n, n), dtype=complex)
    for a, b, c in product(range(n), repeat=3):
        f_np[a, b, c] = complex(sp.N(f[a][b][c], 30))
    for a, b in product(range(n), repeat=2):
        lhs = gens[a] @ gens[b] - gens[b] @ gens[a]
        rhs = 1j * sum(f_np[a, b, c] * gens[c] for c in range(n))
        r = max(r, float(np.linalg.norm(lhs - rhs)))
    return r


def cross_commutator_residual(gens_a: list[np.ndarray], gens_b: list[np.ndarray]) -> float:
    r = 0.0
    for ga in gens_a:
        for gb in gens_b:
            r = max(r, float(np.linalg.norm(ga @ gb - gb @ ga)))
    return r


def embedding_audit(
    q2190: dict[str, Any],
    assignment: dict[str, Any],
    su3_exact: list[sp.Matrix],
    su2_exact: list[sp.Matrix],
    f_su3: list[list[list[sp.Expr]]],
    f_su2: list[list[list[sp.Expr]]],
) -> dict[str, Any]:
    if "_missing" in q2190 or "_missing" in assignment:
        return {
            "status": "OPEN_MISSING_INPUT",
            "q2190_present": "_missing" not in q2190,
            "assignment_present": "_missing" not in assignment,
        }

    n = int(q2190["mode_mapping"]["n_octaves"])
    kpars = {key: float(value) for key, value in q2190["kernel"].items()}
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
    su3_emb = embed_generators(b3, exact_mats_to_numpy(su3_exact))
    su2_emb = embed_generators(b2, exact_mats_to_numpy(su2_exact))

    tolerance = 1e-12
    audits = {
        "b3_orthonormal_residual": orthonormal_residual(b3),
        "b2_orthonormal_residual": orthonormal_residual(b2),
        "kernel_invariance_residual_su3": inv_residual(ktotal, b3),
        "kernel_invariance_residual_su2": inv_residual(ktotal, b2),
        "embedded_su3_closure_residual": closure_residual(su3_emb, f_su3),
        "embedded_su2_closure_residual": closure_residual(su2_emb, f_su2),
        "embedded_su3_su2_cross_commutator_residual": cross_commutator_residual(su3_emb, su2_emb),
    }
    return {
        "status": "PASS_NUMERIC_QW2190_EMBEDDING_AUDIT" if all(v <= tolerance for v in audits.values()) else "FAIL",
        "tolerance": tolerance,
        "audits": audits,
        "all_audits_within_tolerance": all(v <= tolerance for v in audits.values()),
        "scope_guard": (
            "This is a QW-2190 embedding-carrier audit for SU(3)/SU(2) algebra data. "
            "It is not a global selector theorem and not a BV/BRST theorem."
        ),
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)

    q2190 = load_json(QW2190_ROOT)
    assignment = load_json(ASSIGNMENT)
    p1959 = load_json(P1959)
    p452_summary = load_json(P452_SUMMARY)
    p454_summary = load_json(P454_SUMMARY)

    su3_gens = gell_mann_generators_exact()
    su2_gens = pauli_generators_exact()
    f_su3 = structure_constants(su3_gens)
    f_su2 = structure_constants(su2_gens)
    su3_rows = sparse_rows(f_su3, "SU3")
    su2_rows = sparse_rows(f_su2, "SU2")

    su3_generator_cert = generator_certificate(su3_gens, f_su3)
    su2_generator_cert = generator_certificate(su2_gens, f_su2)
    su3_jacobi = jacobi_certificate(f_su3)
    su2_jacobi = jacobi_certificate(f_su2)

    embedding = embedding_audit(q2190, assignment, su3_gens, su2_gens, f_su3, f_su2)

    exact_algebra_pass = all(
        [
            su3_generator_cert["normalization_2Tr_TaTb_delta"],
            su3_generator_cert["hermitian_generators"],
            su3_generator_cert["traceless_generators"],
            su3_generator_cert["closure_reconstruction_exact_zero"],
            su3_generator_cert["f_abc_total_antisymmetry_exact_zero"],
            su3_jacobi["all_zero"],
            su2_generator_cert["normalization_2Tr_TaTb_delta"],
            su2_generator_cert["hermitian_generators"],
            su2_generator_cert["traceless_generators"],
            su2_generator_cert["closure_reconstruction_exact_zero"],
            su2_generator_cert["f_abc_total_antisymmetry_exact_zero"],
            su2_jacobi["all_zero"],
        ]
    )

    q2190_root_present = QW2190_ROOT.exists()
    p452_pass = bool(p452_summary.get("audits_pass") or p452_summary.get("status", "").startswith("PASS"))
    p454_pass = bool(p454_summary.get("audits_pass") or p454_summary.get("status", "").startswith("PASS"))
    embedding_pass = bool(embedding.get("all_audits_within_tolerance"))

    terms = [
        "P1960",
        "StructureConstants_SU3_SU2_U1_strict_v1",
        "JacobiCertificate_SU3_SU2_U1_strict_v1",
        "structure constants",
        "Jacobi identity",
        "QW-2190 source report",
        "stałe struktury",
        "stale struktury",
        "tozsamosc Jacobiego",
        "nieabelowy BRST",
        "non-Abelian BRST",
    ]

    F_TABLE = exact_algebra_pass
    JACOBI = su3_jacobi["all_zero"] and su2_jacobi["all_zero"]
    QW2190_REPORT = q2190_root_present
    EMBEDDING = embedding_pass and p452_pass and p454_pass
    CONNECTION_RULES = False
    BRST_RULES = False
    GHOST_SELF = False
    BV_MAP = False
    GHOST_CUT = False

    brst_ready = all(
        [
            F_TABLE,
            JACOBI,
            QW2190_REPORT,
            EMBEDDING,
            CONNECTION_RULES,
            BRST_RULES,
            GHOST_SELF,
            BV_MAP,
            GHOST_CUT,
        ]
    )

    structure_constants_table = {
        "normalization": "Hermitian generators T_a=lambda_a/2 and t_i=sigma_i/2 with 2*Tr(T_a*T_b)=delta_ab",
        "commutator_convention": "[T_a,T_b] = I*f_abc*T_c",
        "global_factor_index_sets": {
            "SU3": "a,b,c in 1..8",
            "SU2": "i,j,k in 1..3",
            "U1": "single hypercharge generator Y with all f=0",
            "cross_factor_constants": "zero for SU3-SU2-U1 direct-product factors",
        },
        "sparse_nonzero_full_antisymmetric_rows": su3_rows + su2_rows,
        "u1_zero_structure_constants": True,
        "direct_product_cross_structure_constants_zero": True,
    }

    jacobi_export = {
        "SU3": su3_jacobi,
        "SU2": su2_jacobi,
        "U1": {"checked_components": 1, "nonzero_components": 0, "all_zero": True},
        "cross_factor_jacobi_zero_by_direct_product": True,
        "all_declared_factors_pass": bool(JACOBI),
    }

    out = {
        "packet_id": "P1960",
        "stage_id": "S910",
        "status": "PARTIAL_STRICT_ALGEBRA_DATA_EXPORT__BRST_STILL_OPEN",
        "local_verdict": (
            "QW2190_ROOT_REPORT_INGESTED_AND_SU3_SU2_U1_STRUCTURE_CONSTANTS_JACOBI_EXPORTED__"
            "NONABELIAN_BRST_OPERATOR_MAP_STILL_OPEN"
        ),
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "route": "strict_only",
        "legacy_bridge_used": False,
        "higher_reasoning_used": True,
        "ignored_files": sorted(SKIP_NAMES),
        "pre_execution_grep_summary": {
            "english_terms": [
                "P1960",
                "StructureConstants_SU3_SU2_U1_strict_v1",
                "JacobiCertificate_SU3_SU2_U1_strict_v1",
                "structure constants",
                "Jacobi identity",
                "QW-2190 source report",
                "non-Abelian BRST",
            ],
            "polish_terms": [
                "stale struktury",
                "tozsamosc Jacobiego",
                "nieabelowy BRST",
            ],
            "term_scan_counts": scan_terms(terms),
            "key_existing_sources_found": [
                "report_qw2190_kernel_mode_representation_emergence_gate.json exists at repository root and is used by P452/P454.",
                "N493/N494/N495 scope QW-2191 rotations/signs as gauge/conjugation for QW-2190 embedding audits only.",
                "P1959 correctly left full non-Abelian BRST open, but its QW-2190 source-report check was FAR-local.",
            ],
        },
        "input_paths": {
            "qw2190_root_report": rel(QW2190_ROOT),
            "qw2190_root_report_present": q2190_root_present,
            "qw2190_far_local_report": rel(QW2190_FAR_LOCAL),
            "qw2190_far_local_report_present": QW2190_FAR_LOCAL.exists(),
            "mode_index_assignment": rel(ASSIGNMENT),
            "p452_summary": rel(P452_SUMMARY),
            "p454_summary": rel(P454_SUMMARY),
            "p1959_json": rel(P1959),
        },
        "input_hashes": {
            "qw2190_root_report_sha256": digest_path(QW2190_ROOT),
            "mode_index_assignment_sha256": digest_path(ASSIGNMENT),
            "p452_summary_sha256": digest_path(P452_SUMMARY),
            "p454_summary_sha256": digest_path(P454_SUMMARY),
            "p1959_json_sha256": digest_path(P1959),
        },
        "QW2190SourceReportIngestion_strict_v1": {
            "present": q2190_root_present,
            "path_correction": {
                "p1959_checked_far_local_path": rel(QW2190_FAR_LOCAL),
                "p1960_checked_repo_root_path": rel(QW2190_ROOT),
                "repo_root_path_is_the_path_used_by_P452_P454": True,
            },
            "q2190_verdict": q2190.get("verdict"),
            "q2190_flags_sample": {
                key: q2190.get("flags", {}).get(key)
                for key in [
                    "kernel_mode_basis_declared_deterministically",
                    "embedded_su3_lie_closure",
                    "embedded_su2_lie_closure",
                    "su3_su2_direct_product_cross_commutator_zero",
                    "full_physical_uniqueness_of_mode_index_assignment",
                ]
            },
            "scope_guard": "This ingests QW-2190 for embedding audits, not global physical uniqueness.",
        },
        "StructureConstants_SU3_SU2_U1_strict_v1": structure_constants_table,
        "structure_constants_table": structure_constants_table,
        "StructureConstantGeneratorCertificates_strict_v1": {
            "SU3": su3_generator_cert,
            "SU2": su2_generator_cert,
            "U1": {
                "normalization": "Abelian U(1) factor has zero structure constants in this algebra export.",
                "closure_reconstruction_exact_zero": True,
                "jacobi_exact_zero": True,
            },
            "exact_algebra_pass": exact_algebra_pass,
        },
        "JacobiCertificate_SU3_SU2_U1_strict_v1": jacobi_export,
        "jacobi_identity_certificate": jacobi_export,
        "QW2190EmbeddingNumericalAudit_strict_v1": {
            "p452_summary_status": p452_summary.get("status"),
            "p452_audits_pass": p452_pass,
            "p454_summary_status": p454_summary.get("status"),
            "p454_audits_pass": p454_pass,
            "p1960_embedding_audit": embedding,
        },
        "p1959_delta": {
            "corrected": [
                "N1/QW2190_REPORT is AVAILABLE at repository-root scope used by P452/P454.",
                "N2/F_TABLE is now AVAILABLE for SU(3), SU(2), U(1) direct-product algebra constants.",
                "N3/JACOBI is now AVAILABLE for the exported algebra constants.",
            ],
            "still_open": [
                "N4/non-Abelian gauge-connection and ghost transformation rules as repo field-level export",
                "BV/BRST operator map",
                "BRST charge Q and full Q^2=0 on the strict field bundle",
                "ghost self-interaction action tied to the full L_total/EOM bundle",
                "ghost-cancelled Cutkosky equality",
                "global QW-2191 selector closure",
            ],
            "p1959_truth_assignment_seen": p1959.get("nonabelian_extension_acceptance_formula", {}).get("truth_assignment"),
        },
        "BRSTReadinessAfterP1960_strict_v1": {
            "formula": (
                "F_TABLE & JACOBI & QW2190_REPORT & EMBEDDING & CONNECTION_RULES "
                "& BRST_RULES & GHOST_SELF & BV_MAP & GHOST_CUT"
            ),
            "truth_assignment": {
                "F_TABLE": bool(F_TABLE),
                "JACOBI": bool(JACOBI),
                "QW2190_REPORT": bool(QW2190_REPORT),
                "EMBEDDING": bool(EMBEDDING),
                "CONNECTION_RULES": CONNECTION_RULES,
                "BRST_RULES": BRST_RULES,
                "GHOST_SELF": GHOST_SELF,
                "BV_MAP": BV_MAP,
                "GHOST_CUT": GHOST_CUT,
            },
            "evaluated_brst_ready": brst_ready,
            "safe_conclusion": (
                "P1960 discharges the algebra-data/Jacobi part of P1959 only. "
                "It does not discharge the field-level BV/BRST ghost-sector obstruction."
            ),
        },
        "false_pass_guard": {
            "does_not_claim": [
                "global QW-2191 discharge",
                "strict-core selector closure",
                "full SU(3)xSU(2)xU(1) BV/BRST operator map",
                "Q^2=0 on the complete strict field bundle",
                "ghost-cancelled Cutkosky theorem",
                "TG2_BRST_GLOBAL_NILPOTENCY PASS",
                "TG3_CUTKOSKY_GLOBAL_UNITARITY PASS",
                "ToE closure",
            ],
            "may_claim": [
                "repo-root QW-2190 source report is present and ingested",
                "exact SU(3), SU(2), and U(1) direct-product structure constants are exported",
                "exact Jacobi certificates for those constants are exported",
                "QW-2190 embedding numerical residuals pass within declared tolerance",
            ],
        },
        "output_digest_sha256": digest_obj(
            {
                "structure_constants": structure_constants_table,
                "jacobi": jacobi_export,
                "embedding": embedding,
            }
        ),
        "next_honest_step": (
            "Build P1961 with high reasoning: construct the field-level non-Abelian BRST differential "
            "s A_mu^a, s c^a, s cbar^a, s B^a from the P1960 constants and P1958 gauge-fixing seed, "
            "then execute an explicit graded nilpotency check. If the strict L_total field registry cannot "
            "support that map, export the precise field-registry obstruction."
        ),
        "lay_explanation": (
            "Znaleziono brakujacy raport QW-2190 w katalogu glownym repo i zbudowano dokladna "
            "tabele reguł komutatorow SU(3), SU(2) oraz zerowa czesc U(1). Sprawdzono tez "
            "tozsamosc Jacobiego, czyli warunek, bez ktorego pelne duchy BRST nie moga dzialac. "
            "To nadal nie jest pelny BRST: mamy teraz algebre, ale nie cala maszyne dzialania na polach."
        ),
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
