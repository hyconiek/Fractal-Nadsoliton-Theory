#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

OUT_JSON = GENERATED / "p424_current_strict_qw2118_host_pair1_block_isotropy_audit_probe.json"
OUT_SUMMARY = GENERATED / "p424_current_strict_qw2118_host_pair1_block_isotropy_audit_probe_summary.json"


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def build_host_kernel_matrix(distance_profile: dict[str, float], n_octaves: int) -> list[list[float]]:
    rows: list[list[float]] = []
    for i in range(n_octaves):
        row: list[float] = []
        for j in range(n_octaves):
            if i == j:
                row.append(0.0)
                continue
            d = min(abs(i - j), n_octaves - abs(i - j))
            row.append(float(distance_profile[str(d)]))
        rows.append(row)
    return rows


def build_host_operator_matrix(kernel_rows: list[list[float]], broken_floor: float) -> list[list[float]]:
    rows: list[list[float]] = []
    for i, row in enumerate(kernel_rows):
        rows.append([float(entry + broken_floor) if i == j else float(entry) for j, entry in enumerate(row)])
    return rows


def pair1_basis(n: int) -> tuple[list[float], list[float]]:
    c1: list[float] = []
    s1: list[float] = []
    scale = math.sqrt(2.0 / n)
    for i in range(n):
        ang = 2.0 * math.pi * i / n
        c1.append(scale * math.cos(ang))
        s1.append(scale * math.sin(ang))
    return c1, s1


def dot(u: list[float], v: list[float]) -> float:
    return float(sum(ui * vi for ui, vi in zip(u, v)))


def matvec(a: list[list[float]], x: list[float]) -> list[float]:
    return [float(sum(aij * xj for aij, xj in zip(row, x))) for row in a]


def quad_block(a: list[list[float]], u: list[float], v: list[float]) -> list[list[float]]:
    au = matvec(a, u)
    av = matvec(a, v)
    a11 = dot(u, au)
    a12 = dot(u, av)
    a22 = dot(v, av)
    return [[float(a11), float(a12)], [float(a12), float(a22)]]


def delta_1(mat2: list[list[float]]) -> tuple[float, float]:
    a_1 = float(mat2[0][0])
    b_1 = float(mat2[0][1])
    d_1 = float(mat2[1][1])
    return float(a_1 - d_1), float(b_1)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    q2118_path = REPO / "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2118_ktotal_spectral_tripartition_gate.json"
    q2124_path = REPO / "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2124_scalar_vacuum_closure_branch_resolved_gate.json"
    p10_path = REPO / "fundamental_action_reconstruction/generated/p10_existing_kernel_feedback_to_kobs_rerun_after_explicit_current_pair_chain.json"

    q2118 = load_json(q2118_path)
    q2124 = load_json(q2124_path)

    n = int(q2118["matrix_spec"]["n_octaves"])
    distance_profile = {str(k): float(v) for k, v in q2118["distance_profile"].items()}
    broken_floor = float(q2124["inputs"]["broken_floor"])

    k_rows = build_host_kernel_matrix(distance_profile, n)
    a_rows = build_host_operator_matrix(k_rows, broken_floor)

    c1, s1 = pair1_basis(n)
    host_pair1 = quad_block(a_rows, c1, s1)
    host_delta = delta_1(host_pair1)

    tol = 1e-12
    host_isotropic = abs(host_delta[0]) <= tol and abs(host_delta[1]) <= tol

    comparison: dict[str, Any] = {"p10_loaded": False}
    if p10_path.exists():
        p10 = load_json(p10_path)
        p10_mat = p10["computed"]["full_current_pair_h3_projected_block"]["matrix"]
        p10_delta = delta_1(p10_mat)
        comparison = {
            "p10_loaded": True,
            "p10_current_pair_h3_block": {
                "matrix": p10_mat,
                "Delta_1": [p10_delta[0], p10_delta[1]],
            },
            "host_pair1_can_equal_p10_block": bool(host_isotropic and abs(p10_delta[0]) <= tol and abs(p10_delta[1]) <= tol),
            "reason_if_not_equal": (
                "host_pair1_is_isotropic_but_p10_block_is_anisotropic"
                if host_isotropic and (abs(p10_delta[0]) > tol or abs(p10_delta[1]) > tol)
                else "unknown_or_missing"
            ),
        }

    artifact = {
        "probe_id": "P424_CURRENT_STRICT_QW2118_HOST_PAIR1_BLOCK_ISOTROPY_AUDIT_PROBE",
        "as_of": "2026-03-12",
        "no_false_pass": True,
        "sources": {
            "q2118": str(q2118_path.relative_to(REPO)),
            "q2124": str(q2124_path.relative_to(REPO)),
            "p10": str(p10_path.relative_to(REPO)),
        },
        "host_operator": {
            "symbol": "A_host = K_total + m0^2 I",
            "n_octaves": n,
            "broken_floor_m0_sq": broken_floor,
            "pair1_basis_order": ["c1", "s1"],
            "pair1_block": host_pair1,
            "Delta_1": [host_delta[0], host_delta[1]],
            "tolerance": tol,
            "is_isotropic_within_tolerance": host_isotropic,
        },
        "comparison": comparison,
        "theorem_pointer": "N465",
    }

    summary = {
        "probe_id": artifact["probe_id"],
        "as_of": artifact["as_of"],
        "status": "EXECUTED_HOST_PAIR1_BLOCK_ISOTROPY_AUDIT_NO_FALSE_PASS",
        "host_pair1_isotropic": host_isotropic,
        "host_pair1_Delta_1": artifact["host_operator"]["Delta_1"],
        "p10_loaded": comparison.get("p10_loaded", False),
        "host_pair1_can_equal_p10_block": comparison.get("host_pair1_can_equal_p10_block"),
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

