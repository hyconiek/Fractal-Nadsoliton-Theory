#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

OUT_JSON = GENERATED / "p427_current_strict_pair1_diagonal_o2_cut_canonical_eigenbasis_audit_probe.json"
OUT_SUMMARY = GENERATED / "p427_current_strict_pair1_diagonal_o2_cut_canonical_eigenbasis_audit_probe_summary.json"


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def pair1_basis(n: int) -> tuple[list[float], list[float]]:
    c1: list[float] = []
    s1: list[float] = []
    scale = math.sqrt(2.0 / n)
    for k in range(n):
        ang = 2.0 * math.pi * k / n
        c1.append(scale * math.cos(ang))
        s1.append(scale * math.sin(ang))
    return c1, s1


def rotate_pair1_basis(c1: list[float], s1: list[float], theta: float) -> tuple[list[float], list[float]]:
    ct = math.cos(theta)
    st = math.sin(theta)
    c1p = [float(ct * c + st * s) for c, s in zip(c1, s1)]
    s1p = [float(-st * c + ct * s) for c, s in zip(c1, s1)]
    return c1p, s1p


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


def f2(d: list[float], n: int) -> tuple[float, float, float]:
    re = 0.0
    im = 0.0
    for k, dk in enumerate(d):
        ang = 4.0 * math.pi * k / n
        re += float(dk) * math.cos(ang)
        im += float(dk) * math.sin(ang)
    mag = math.sqrt(re * re + im * im)
    return float(re), float(im), float(mag)


def eigenvalues_2x2_sym(mat: list[list[float]]) -> tuple[float, float]:
    a = float(mat[0][0])
    b = float(mat[0][1])
    d = float(mat[1][1])
    tr = a + d
    disc = (a - d) * (a - d) + 4.0 * b * b
    root = math.sqrt(max(0.0, disc))
    lam_plus = 0.5 * (tr + root)
    lam_minus = 0.5 * (tr - root)
    if lam_plus < lam_minus:
        lam_plus, lam_minus = lam_minus, lam_plus
    return float(lam_plus), float(lam_minus)


def profile_constant(n: int, value: float) -> list[float]:
    return [float(value) for _ in range(n)]


def profile_mode1_cos(n: int, amplitude: float) -> list[float]:
    return [float(amplitude * math.cos(2.0 * math.pi * k / n)) for k in range(n)]


def profile_mode2_cos(n: int, amplitude: float) -> list[float]:
    return [float(amplitude * math.cos(4.0 * math.pi * k / n)) for k in range(n)]


def profile_mode2_sin(n: int, amplitude: float) -> list[float]:
    return [float(amplitude * math.sin(4.0 * math.pi * k / n)) for k in range(n)]


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
        ("mixed", [1.0 + 0.3 * profile_mode2_cos(n, 1.0)[k] + 0.4 * profile_mode2_sin(n, 1.0)[k] for k in range(n)]),
    ]

    rows: list[dict[str, Any]] = []
    for name, d in profiles:
        re_f2, im_f2, mag_f2 = f2(d, n)
        breaks_o2 = bool(mag_f2 > tol)

        mat = quad_block_diag(d, c1, s1)
        lam_plus_mat, lam_minus_mat = eigenvalues_2x2_sym(mat)

        mean_d = float(sum(d) / n)
        lam_plus_pred = float(mean_d + (mag_f2 / n))
        lam_minus_pred = float(mean_d - (mag_f2 / n))

        theta_star: float | None
        if not breaks_o2:
            theta_star = None
            rot_block = mat
            offdiag = abs(float(rot_block[0][1]))
        else:
            theta_star = 0.5 * math.atan2(im_f2, re_f2)
            c1p, s1p = rotate_pair1_basis(c1, s1, theta_star)
            rot_block = quad_block_diag(d, c1p, s1p)

            # Canonicalize by ordering eigenvalues (swap via θ -> θ + π/2 if needed).
            if float(rot_block[0][0]) < float(rot_block[1][1]):
                theta_star = float(theta_star + 0.5 * math.pi)
                c1p, s1p = rotate_pair1_basis(c1, s1, theta_star)
                rot_block = quad_block_diag(d, c1p, s1p)

            offdiag = abs(float(rot_block[0][1]))

        rows.append(
            {
                "name": name,
                "diag_profile": [float(x) for x in d],
                "F2": {"Re": re_f2, "Im": im_f2, "abs": mag_f2},
                "breaks_O2_on_pair1": breaks_o2,
                "pair1_block_c1s1": mat,
                "theta_star": theta_star,
                "pair1_block_rotated_by_theta_star": rot_block,
                "offdiag_abs_after_rotation": offdiag,
                "eigenvalues_from_block": [lam_plus_mat, lam_minus_mat],
                "eigenvalues_predicted_from_F2": [lam_plus_pred, lam_minus_pred],
                "abs_error_eigs": [abs(lam_plus_mat - lam_plus_pred), abs(lam_minus_mat - lam_minus_pred)],
            }
        )

    artifact = {
        "probe_id": "P427_CURRENT_STRICT_PAIR1_DIAGONAL_O2_CUT_CANONICAL_EIGENBASIS_AUDIT_PROBE",
        "as_of": "2026-03-12",
        "no_false_pass": True,
        "sources": {"q2118": str(q2118_path.relative_to(REPO))},
        "n_octaves": n,
        "pair1_basis_order": ["c1", "s1"],
        "tolerance": tol,
        "profiles": rows,
        "theorem_pointer": "N468",
    }

    summary = {
        "probe_id": artifact["probe_id"],
        "as_of": artifact["as_of"],
        "status": "EXECUTED_PAIR1_DIAGONAL_CANONICAL_EIGENBASIS_AUDIT_NO_FALSE_PASS",
        "n_octaves": n,
        "tolerance": tol,
        "profile_results": [
            {
                "name": row["name"],
                "breaks_O2_on_pair1": row["breaks_O2_on_pair1"],
                "F2_abs": row["F2"]["abs"],
                "theta_star": row["theta_star"],
                "offdiag_abs_after_rotation": row["offdiag_abs_after_rotation"],
                "abs_error_eigs": row["abs_error_eigs"],
            }
            for row in rows
        ],
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

