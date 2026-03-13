#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"

OUT_JSON = GENERATED / "p430_current_strict_aut_z12_invariant_diagonal_profile_mode2_defect_orbit_reduction_audit_probe.json"
OUT_SUMMARY = (
    GENERATED / "p430_current_strict_aut_z12_invariant_diagonal_profile_mode2_defect_orbit_reduction_audit_probe_summary.json"
)


def pair1_basis(n: int) -> tuple[list[float], list[float]]:
    c1: list[float] = []
    s1: list[float] = []
    scale = math.sqrt(2.0 / n)
    for k in range(n):
        ang = 2.0 * math.pi * k / n
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


def f2(d: list[float], n: int) -> tuple[float, float, float]:
    re = 0.0
    im = 0.0
    for k, dk in enumerate(d):
        ang = 4.0 * math.pi * k / n
        re += float(dk) * math.cos(ang)
        im += float(dk) * math.sin(ang)
    mag = math.sqrt(re * re + im * im)
    return float(re), float(im), float(mag)


def orbit_profile_aut_z12_invariant(
    n: int,
    *,
    t0: float,
    t6: float,
    t1: float,
    t2: float,
    t3: float,
    t4: float,
) -> list[float]:
    if n != 12:
        raise ValueError("This probe is specialized to n=12.")
    d = [0.0] * 12
    orbit_map: dict[int, float] = {
        0: float(t0),
        6: float(t6),
        1: float(t1),
        5: float(t1),
        7: float(t1),
        11: float(t1),
        2: float(t2),
        10: float(t2),
        3: float(t3),
        9: float(t3),
        4: float(t4),
        8: float(t4),
    }
    for k in range(12):
        d[k] = float(orbit_map[k])
    return d


def f2_orbit_formula(*, t0: float, t6: float, t1: float, t2: float, t3: float, t4: float) -> float:
    # N471: F2(d) = t0 + t6 + 2*t1 - t2 - 2*t3 - t4, and Im(F2)=0.
    return float(t0 + t6 + 2.0 * t1 - t2 - 2.0 * t3 - t4)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n = 12
    tol = 1e-12

    c1, s1 = pair1_basis(n)

    examples = [
        {
            "name": "example_nonzero_real_defect",
            "orbit_values": {"t0": 1.0, "t6": 2.0, "t1": 0.5, "t2": -0.7, "t3": 0.25, "t4": 3.0},
        },
        {
            "name": "example_parity_only_profile",
            "orbit_values": {"t0": 1.2, "t6": 1.2, "t1": -0.3, "t2": 1.2, "t3": -0.3, "t4": 1.2},
        },
        {
            "name": "example_balanced_orbit_profile",
            "orbit_values": {"t0": 0.0, "t6": 0.0, "t1": 1.0, "t2": 2.0, "t3": 0.0, "t4": 0.0},
        },
    ]

    rows: list[dict[str, Any]] = []
    for ex in examples:
        vals = ex["orbit_values"]
        d = orbit_profile_aut_z12_invariant(n, **vals)
        re_f2, im_f2, mag_f2 = f2(d, n)
        predicted_re = f2_orbit_formula(**vals)
        predicted_im = 0.0

        block = quad_block_diag(d, c1, s1)
        a1, b1, d1 = float(block[0][0]), float(block[0][1]), float(block[1][1])

        rows.append(
            {
                "name": ex["name"],
                "orbit_values": vals,
                "diag_profile": d,
                "F2_computed": {"Re": re_f2, "Im": im_f2, "abs": mag_f2},
                "F2_predicted_by_N471": {"Re": predicted_re, "Im": predicted_im, "abs": abs(predicted_re)},
                "errors": {
                    "abs_Re_error": abs(re_f2 - predicted_re),
                    "abs_Im_error": abs(im_f2 - predicted_im),
                    "abs_abs_error": abs(mag_f2 - abs(predicted_re)),
                },
                "passes_tolerance": {
                    "Im_below_tol": bool(abs(im_f2) <= tol),
                    "Re_matches_formula": bool(abs(re_f2 - predicted_re) <= tol),
                    "abs_matches_formula": bool(abs(mag_f2 - abs(predicted_re)) <= tol),
                },
                "cuts_O2_on_pair1_by_N466": bool(mag_f2 > tol),
                "pair1_block_c1s1": block,
                "pair1_anisotropy_signature": {"a1_minus_d1": float(a1 - d1), "b1": b1},
            }
        )

    artifact = {
        "probe_id": "P430_CURRENT_STRICT_AUT_Z12_INVARIANT_DIAGONAL_PROFILE_MODE2_DEFECT_ORBIT_REDUCTION_AUDIT_PROBE",
        "as_of": "2026-03-12",
        "no_false_pass": True,
        "n_octaves": n,
        "tolerance": tol,
        "profiles": rows,
        "theorem_pointer": "N471",
        "notes": [
            "This probe audits the pure-math orbit-reduction identity for F2(d) under Aut(Z_12)-invariance on n=12.",
            "It does not claim that any physically exported diagonal/local sector is Aut-invariant.",
        ],
    }

    summary = {
        "probe_id": artifact["probe_id"],
        "as_of": artifact["as_of"],
        "status": "EXECUTED_AUT_Z12_INVARIANT_DIAGONAL_MODE2_DEFECT_ORBIT_REDUCTION_AUDIT_NO_FALSE_PASS",
        "n_octaves": n,
        "tolerance": tol,
        "profile_results": [
            {
                "name": row["name"],
                "F2_abs": row["F2_computed"]["abs"],
                "F2_Re": row["F2_computed"]["Re"],
                "F2_Im": row["F2_computed"]["Im"],
                "passes_tolerance": row["passes_tolerance"],
                "cuts_O2_on_pair1_by_N466": row["cuts_O2_on_pair1_by_N466"],
            }
            for row in rows
        ],
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

