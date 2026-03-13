#!/usr/bin/env python3
from __future__ import annotations

import cmath
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"

OUT_JSON = GENERATED / "p431_current_strict_canonical_local_diagonal_mode2_defect_underdetermination_audit_probe.json"
OUT_SUMMARY = GENERATED / "p431_current_strict_canonical_local_diagonal_mode2_defect_underdetermination_audit_probe_summary.json"


def f2(d: list[float], n: int) -> complex:
    return sum(complex(dk) * cmath.exp(1j * (4.0 * math.pi * k / n)) for k, dk in enumerate(d))


def pair1_anisotropy_from_f2(f2_val: complex, n: int) -> dict[str, float]:
    # N466: Δ1(D) = (a1-d1, b1) = (2/n Re(F2), 1/n Im(F2)).
    # For n=12 this matches P426: (1/6 Re(F2), 1/12 Im(F2)).
    return {"a1_minus_d1": float((2.0 / n) * f2_val.real), "b1": float((1.0 / n) * f2_val.imag)}


def constant_profile(n: int, c: float) -> list[float]:
    return [float(c)] * n


def mode2_cos_profile(n: int) -> list[float]:
    # d_k := cos(4πk/n)  (pure mode-2 cosine profile for the n-ring)
    return [float(math.cos(4.0 * math.pi * k / n)) for k in range(n)]


def realize_via_r15_free_m2_only(d: list[float], *, m0_sq: float) -> dict[str, Any]:
    # R15 class: d_k = (3*g4_k*vpsi_k^2 + 5*g6_k*vpsi_k^4 + 2*gY_k*vphi^2 + m2_k) - m0^2.
    # Witness specialization: set g4=g6=gY=0 and vpsi=vphi=0, so d_k = m2_k - m0^2.
    m2 = [float(m0_sq + dk) for dk in d]
    return {
        "specialization": {
            "g4_psi_k": 0.0,
            "g6_psi_k": 0.0,
            "gY_k": 0.0,
            "vpsi_k": 0.0,
            "vphi": 0.0,
        },
        "implied_assignment": {"m2_psi_k": m2},
        "checks": {"d_k_equals_m2_minus_m0_sq": True},
    }


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n = 12
    tol = 1e-12

    # Exported by R15/QW-2124 as the certified host scalar floor.
    m0_sq = 1.013551972358388

    profiles = [
        {
            "name": "witness_F2_zero_constant_profile",
            "d": constant_profile(n, 0.0),
        },
        {
            "name": "witness_F2_nonzero_mode2_cos_profile",
            "d": mode2_cos_profile(n),
        },
    ]

    rows: list[dict[str, Any]] = []
    for p in profiles:
        d = p["d"]
        f2_val = f2(d, n)
        rows.append(
            {
                "name": p["name"],
                "diag_profile_d": d,
                "F2": {"Re": float(f2_val.real), "Im": float(f2_val.imag), "abs": float(abs(f2_val))},
                "cuts_O2_on_pair1_by_N466": bool(abs(f2_val) > tol),
                "pair1_anisotropy_signature_from_N466": pair1_anisotropy_from_f2(f2_val, n),
                "r15_realization_witness": realize_via_r15_free_m2_only(d, m0_sq=m0_sq),
            }
        )

    artifact = {
        "probe_id": "P431_CURRENT_STRICT_CANONICAL_LOCAL_DIAGONAL_MODE2_DEFECT_UNDERDETERMINATION_AUDIT_PROBE",
        "as_of": "2026-03-12",
        "no_false_pass": True,
        "n_octaves": n,
        "tolerance": tol,
        "exported_host_floor_m0_sq": m0_sq,
        "witness_profiles": rows,
        "notes": [
            "This probe demonstrates underdetermination of F2(d) at the level of the exported R15 coefficient class.",
            "It does not claim any physical admissibility of the witness assignments (no vacuum/EOM solving).",
        ],
    }

    summary = {
        "probe_id": artifact["probe_id"],
        "as_of": artifact["as_of"],
        "status": "EXECUTED_CANONICAL_LOCAL_DIAGONAL_MODE2_DEFECT_UNDERDETERMINATION_AUDIT_NO_FALSE_PASS",
        "n_octaves": n,
        "tolerance": tol,
        "witness_results": [
            {
                "name": row["name"],
                "F2_abs": row["F2"]["abs"],
                "F2_Re": row["F2"]["Re"],
                "F2_Im": row["F2"]["Im"],
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

