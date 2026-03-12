#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

OUT_JSON = GENERATED / "p425_current_strict_pair1_diagonal_profile_o2_cut_audit_probe.json"
OUT_SUMMARY = GENERATED / "p425_current_strict_pair1_diagonal_profile_o2_cut_audit_probe_summary.json"


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


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


def diag_matvec(d: list[float], x: list[float]) -> list[float]:
    return [float(di * xi) for di, xi in zip(d, x)]


def quad_block_diag(d: list[float], u: list[float], v: list[float]) -> list[list[float]]:
    du = diag_matvec(d, u)
    dv = diag_matvec(d, v)
    a1 = dot(u, du)
    b1 = dot(u, dv)
    d1 = dot(v, dv)
    return [[float(a1), float(b1)], [float(b1), float(d1)]]


def delta_1(mat2: list[list[float]]) -> tuple[float, float]:
    a1 = float(mat2[0][0])
    b1 = float(mat2[0][1])
    d1 = float(mat2[1][1])
    return float(a1 - d1), float(b1)


def delta_1_closed_form(d: list[float], n: int) -> tuple[float, float]:
    s_cos = 0.0
    s_sin = 0.0
    for i, di in enumerate(d):
        ang = 4.0 * math.pi * i / n
        s_cos += float(di) * math.cos(ang)
        s_sin += float(di) * math.sin(ang)
    return float((2.0 / n) * s_cos), float((1.0 / n) * s_sin)


def profile_constant(n: int, value: float) -> list[float]:
    return [float(value) for _ in range(n)]


def profile_mode2_cos(n: int, amplitude: float) -> list[float]:
    return [float(amplitude * math.cos(4.0 * math.pi * i / n)) for i in range(n)]


def profile_mode2_sin(n: int, amplitude: float) -> list[float]:
    return [float(amplitude * math.sin(4.0 * math.pi * i / n)) for i in range(n)]


def profile_mode1_cos(n: int, amplitude: float) -> list[float]:
    return [float(amplitude * math.cos(2.0 * math.pi * i / n)) for i in range(n)]


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    q2118_path = REPO / "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2118_ktotal_spectral_tripartition_gate.json"
    q2118 = load_json(q2118_path)
    n = int(q2118["matrix_spec"]["n_octaves"])

    c1, s1 = pair1_basis(n)
    tol = 1e-12

    profiles = [
        ("constant_1", profile_constant(n, 1.0)),
        ("mode1_cos_amp_1", profile_mode1_cos(n, 1.0)),
        ("mode2_cos_amp_1", profile_mode2_cos(n, 1.0)),
        ("mode2_sin_amp_1", profile_mode2_sin(n, 1.0)),
        ("mixed", [1.0 + 0.3 * profile_mode2_cos(n, 1.0)[i] + 0.4 * profile_mode2_sin(n, 1.0)[i] for i in range(n)]),
    ]

    rows: list[dict[str, Any]] = []
    for name, d in profiles:
        mat = quad_block_diag(d, c1, s1)
        delta = delta_1(mat)
        delta_cf = delta_1_closed_form(d, n)
        rows.append(
            {
                "name": name,
                "diag_profile": [float(x) for x in d],
                "pair1_block": mat,
                "Delta_1_computed": [delta[0], delta[1]],
                "Delta_1_closed_form": [delta_cf[0], delta_cf[1]],
                "abs_error": [abs(delta[0] - delta_cf[0]), abs(delta[1] - delta_cf[1])],
                "breaks_O2_on_pair1": bool(abs(delta[0]) > tol or abs(delta[1]) > tol),
            }
        )

    artifact = {
        "probe_id": "P425_CURRENT_STRICT_PAIR1_DIAGONAL_PROFILE_O2_CUT_AUDIT_PROBE",
        "as_of": "2026-03-12",
        "no_false_pass": True,
        "sources": {"q2118": str(q2118_path.relative_to(REPO))},
        "n_octaves": n,
        "pair1_basis_order": ["c1", "s1"],
        "tolerance": tol,
        "profiles": rows,
        "theorem_pointer": "N466",
    }

    summary = {
        "probe_id": artifact["probe_id"],
        "as_of": artifact["as_of"],
        "status": "EXECUTED_PAIR1_DIAGONAL_PROFILE_O2_CUT_AUDIT_NO_FALSE_PASS",
        "n_octaves": n,
        "tolerance": tol,
        "profile_results": [
            {
                "name": row["name"],
                "breaks_O2_on_pair1": row["breaks_O2_on_pair1"],
                "Delta_1": row["Delta_1_computed"],
                "abs_error": row["abs_error"],
            }
            for row in rows
        ],
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

