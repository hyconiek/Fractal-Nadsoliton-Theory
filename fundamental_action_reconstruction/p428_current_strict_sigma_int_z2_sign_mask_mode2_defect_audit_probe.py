#!/usr/bin/env python3
from __future__ import annotations

import cmath
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_MASK = GENERATED / "b_sigma_int_E_pair_sign_mask_strict_provenance_v1.json"

OUT_JSON = GENERATED / "p428_current_strict_sigma_int_z2_sign_mask_mode2_defect_audit_probe.json"
OUT_SUMMARY = GENERATED / "p428_current_strict_sigma_int_z2_sign_mask_mode2_defect_audit_probe_summary.json"


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def fm(d: list[float], m: int) -> complex:
    n = len(d)
    total = 0.0 + 0.0j
    for k, dk in enumerate(d):
        ang = 2.0 * math.pi * m * k / n
        total += float(dk) * complex(math.cos(ang), math.sin(ang))
    return complex(total)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    mask = load_json(IN_MASK)
    tol = 1e-12

    rows: list[dict[str, Any]] = []
    for pair_key in ("pair1", "pair2"):
        d_raw = mask["value"][pair_key]
        d = [float(x) for x in d_raw]
        n = len(d)

        f2 = fm(d, 2)
        f6 = fm(d, 6)

        mean = float(sum(d) / n)
        # From N466: (a1-d1, b1) = ((2/n)Re F2, (1/n)Im F2)
        a1_minus_d1 = float((2.0 / n) * f2.real)
        b1 = float((1.0 / n) * f2.imag)

        rows.append(
            {
                "pair": pair_key,
                "n": n,
                "diag_profile": d_raw,
                "mean": mean,
                "F2": {"Re": float(f2.real), "Im": float(f2.imag), "abs": float(abs(f2))},
                "F6": {"Re": float(f6.real), "Im": float(f6.imag), "abs": float(abs(f6))},
                "pair1_anisotropy_signature_from_N466": {"a1_minus_d1": a1_minus_d1, "b1": b1},
                "cuts_O2_on_pair1_by_N466": bool(abs(f2) > tol),
            }
        )

    artifact = {
        "probe_id": "P428_CURRENT_STRICT_SIGMA_INT_Z2_SIGN_MASK_MODE2_DEFECT_AUDIT_PROBE",
        "as_of": "2026-03-12",
        "no_false_pass": True,
        "input": str(IN_MASK.relative_to(REPO)),
        "tolerance": tol,
        "results": rows,
        "theorem_pointer": "N469",
        "notes": [
            "This probe only audits the strict-exported FR-derived Z2 parity sign mask b_{i,k}=(-1)^k / (-1)^{k+1}.",
            "It does not assert that diag(weights) is a physically exported selector operator; it only audits the mode-2 defect condition relevant for any diagonal pair1 O(2)-cut attempt (N466).",
        ],
    }

    summary = {
        "probe_id": artifact["probe_id"],
        "as_of": artifact["as_of"],
        "status": "EXECUTED_SIGMA_INT_Z2_MASK_MODE2_DEFECT_AUDIT_NO_FALSE_PASS",
        "tolerance": tol,
        "pair_results": [
            {
                "pair": row["pair"],
                "F2_abs": row["F2"]["abs"],
                "cuts_O2_on_pair1_by_N466": row["cuts_O2_on_pair1_by_N466"],
                "F6_abs": row["F6"]["abs"],
            }
            for row in rows
        ],
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

