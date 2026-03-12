#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

OUT_JSON = GENERATED / "p429_current_strict_plus3_shift_invariant_diagonal_profile_mode2_defect_audit_probe.json"
OUT_SUMMARY = GENERATED / "p429_current_strict_plus3_shift_invariant_diagonal_profile_mode2_defect_audit_probe_summary.json"


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


def profile_plus3_invariant(n: int, d0: float, d1: float, d2: float) -> list[float]:
    base = [float(d0), float(d1), float(d2)]
    return [base[k % 3] for k in range(n)]


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n = 12
    tol = 1e-12

    c1, s1 = pair1_basis(n)

    examples = [
        ("example_1", profile_plus3_invariant(n, 1.0, 2.3, -0.7)),
        ("example_2", profile_plus3_invariant(n, 0.0, 1.0, 0.0)),
        ("example_3", profile_plus3_invariant(n, -3.0, 0.5, 4.25)),
    ]

    rows: list[dict[str, Any]] = []
    for name, d in examples:
        re_f2, im_f2, mag_f2 = f2(d, n)
        block = quad_block_diag(d, c1, s1)
        a1, b1, d1 = float(block[0][0]), float(block[0][1]), float(block[1][1])
        rows.append(
            {
                "name": name,
                "diag_profile": d,
                "F2": {"Re": re_f2, "Im": im_f2, "abs": mag_f2},
                "cuts_O2_on_pair1_by_N466": bool(mag_f2 > tol),
                "pair1_block_c1s1": block,
                "pair1_anisotropy_signature": {"a1_minus_d1": float(a1 - d1), "b1": b1},
            }
        )

    artifact = {
        "probe_id": "P429_CURRENT_STRICT_PLUS3_SHIFT_INVARIANT_DIAGONAL_PROFILE_MODE2_DEFECT_AUDIT_PROBE",
        "as_of": "2026-03-12",
        "no_false_pass": True,
        "n_octaves": n,
        "tolerance": tol,
        "profiles": rows,
        "theorem_pointer": "N470",
        "notes": [
            "This probe audits the pure-math fact that +3-shift invariance on n=12 forces the mode-2 defect F2(d) to vanish.",
            "It does not claim that any physically exported diagonal/local sector is +3-shift invariant.",
        ],
    }

    summary = {
        "probe_id": artifact["probe_id"],
        "as_of": artifact["as_of"],
        "status": "EXECUTED_PLUS3_SHIFT_INVARIANT_DIAGONAL_MODE2_DEFECT_AUDIT_NO_FALSE_PASS",
        "n_octaves": n,
        "tolerance": tol,
        "profile_results": [
            {
                "name": row["name"],
                "F2_abs": row["F2"]["abs"],
                "cuts_O2_on_pair1_by_N466": row["cuts_O2_on_pair1_by_N466"],
                "pair1_anisotropy_signature": row["pair1_anisotropy_signature"],
            }
            for row in rows
        ],
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

