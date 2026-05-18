#!/usr/bin/env python3
"""P1967 S917 strict Shannon axis source to Delta_sel tensor map.

Repository grep after P1966 found existing strict, lane-scoped axis-only
symmetry-breaking sources: the Shannon element-order reference mode-index
assignment (F454/N480/N488/N496/N514) and its diagonal/local alignment audit
(P455/N497).  This script performs the next honest step: it maps that strict
source into the exact P1966 minimal traceless tensor shape for pair1..pair5 and
verifies the nonzero O(2)->Z2 spectral gaps.

It deliberately does not claim residual-Z2 sign-sensitive selector closure,
admissible S_sel_int, global QW-2191 discharge, or ToE closure.
"""

from __future__ import annotations

import cmath
import hashlib
import json
import math
import platform
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1967_s917_strict_shannon_axis_source_to_delta_sel_tensor_map.json"

SOURCE_FILES = [
    "N480_CURRENT_FIRST_STRICT_SHANNON_ELEMENT_ORDER_REFERENCE_CROSS_ENTROPY_CUTS_PAIR1_O2_TO_Z2_UNIQUENESS_THEOREM.md",
    "N488_CURRENT_FIRST_STRICT_SHANNON_ELEMENT_ORDER_REFERENCE_CROSS_ENTROPY_CUTS_PAIR2_O2_TO_Z2_UNIQUENESS_THEOREM.md",
    "N496_CURRENT_FIRST_STRICT_SHANNON_ELEMENT_ORDER_REFERENCE_CROSS_ENTROPY_CUTS_PAIR3_TO_PAIR5_O2_TO_Z2_UNIQUENESS_THEOREM.md",
    "N514_CURRENT_FIRST_STRICT_ZN_ELEMENT_ORDER_REFERENCE_FOURIER_DEFECT_NONVANISHING_THEOREM.md",
    "F454_CURRENT_STRICT_SHANNON_ELEMENT_ORDER_REFERENCE_MODE_INDEX_ASSIGNMENT_PACKET.md",
    "P455_CURRENT_STRICT_MODE_INDEX_ASSIGNMENT_SHANNON_VS_DIAGONAL_ALIGNMENT_AUDIT_PROBE.md",
    "B8_SELECTOR_TRACK_ANTI_OVERCLAIM_AUDIT.md",
]

KEYWORDS = [
    "O(2)",
    "Z2",
    "symmetry-breaking",
    "axis",
    "direction",
    "Shannon",
    "element-order",
    "mode-index assignment",
    "strict-core selector closure",
    "QW-2191",
]


def load_generated(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def file_digest(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repo_grep_evidence() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for name in SOURCE_FILES:
        path = ROOT / name
        if not path.exists():
            rows.append({"path": name, "present": False, "matches": []})
            continue
        matches = []
        for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
            lower = line.lower()
            if any(keyword.lower() in lower for keyword in KEYWORDS):
                matches.append({"line": lineno, "text": line.strip()[:240]})
            if len(matches) >= 8:
                break
        rows.append({"path": name, "present": True, "sha256": file_digest(path), "matches": matches})
    return rows


def ord_z12(x: int) -> int:
    return 1 if x == 0 else 12 // math.gcd(x, 12)


def exact_fourier_defect_parts(m: int, n: int = 12) -> tuple[sp.Expr, sp.Expr]:
    re = sp.simplify(
        sum(sp.Integer(ord_z12(x)) * sp.cos(sp.Rational(4 * m * x, n) * sp.pi) for x in range(n))
    )
    im = sp.simplify(
        sum(sp.Integer(ord_z12(x)) * sp.sin(sp.Rational(4 * m * x, n) * sp.pi) for x in range(n))
    )
    return re, im


def numeric_fourier_defect(m: int, n: int = 12) -> complex:
    return sum(ord_z12(x) * cmath.exp(1j * 4 * math.pi * m * x / n) for x in range(n))


def build_pair_rows(shannon_assignment: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    exported_pairs = shannon_assignment.get("outputs", {}).get("pairs", {})
    for m in range(1, 6):
        re_exact, im_exact = exact_fourier_defect_parts(m)
        f_exact = sp.sstr(re_exact) if im_exact == 0 else f"{sp.sstr(re_exact)} + I*({sp.sstr(im_exact)})"
        f_num = numeric_fourier_defect(m)
        re = float(sp.N(re_exact, 30))
        im = float(sp.N(im_exact, 30))
        mat = np.array([[re, im], [im, -re]], dtype=float)
        vals = np.linalg.eigvalsh(mat)
        gap = float(vals[1] - vals[0])
        exported = exported_pairs.get(f"pair{m}", {})
        exported_f = exported.get("F_2m_ord", {})
        exported_abs = float(exported_f.get("abs", float("nan")))
        exported_re = float(exported_f.get("Re", float("nan")))
        exported_im = float(exported_f.get("Im", float("nan")))
        rows.append(
            {
                "pair": f"pair{m}",
                "m": m,
                "ord_Z12_profile": [ord_z12(x) for x in range(12)],
                "F_2m_ord_exact": f_exact,
                "F_2m_ord_re_exact": sp.sstr(re_exact),
                "F_2m_ord_im_exact": sp.sstr(im_exact),
                "F_2m_ord_numeric_direct": {"Re": f_num.real, "Im": f_num.imag, "abs": abs(f_num)},
                "Delta_sel_pair_m": [[re, im], [im, -re]],
                "Delta_sel_trace": float(np.trace(mat)),
                "Delta_sel_eigenvalues": [float(v) for v in vals],
                "Delta_sel_gap": gap,
                "nonzero_gap_pass": bool(gap > 1e-12),
                "theta_star_from_delta": 0.5 * math.atan2(im, re),
                "exported_F454_F_2m_ord": {"Re": exported_re, "Im": exported_im, "abs": exported_abs},
                "matches_exported_F454_abs": bool(abs(abs(f_num) - exported_abs) < 1e-9),
                "matches_exported_F454_re_im": bool(
                    abs(f_num.real - exported_re) < 1e-9 and abs(f_num.imag - exported_im) < 1e-9
                ),
            }
        )
    return rows


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    shannon_assignment = load_generated("mode_index_assignment_shannon_element_order_reference_strict_core_v1.json")
    p455_summary = load_generated("p455_current_strict_mode_index_assignment_shannon_vs_diagonal_alignment_audit_probe_summary.json")
    p1966 = load_generated("p1966_s916_strict_qw2191_selector_premise_obstruction_and_minimal_breaking.json")

    pair_rows = build_pair_rows(shannon_assignment)
    all_gaps_nonzero = all(row["nonzero_gap_pass"] for row in pair_rows)
    all_match_export = all(row["matches_exported_F454_abs"] and row["matches_exported_F454_re_im"] for row in pair_rows)
    p455_aligned = p455_summary.get("aligned_all_pairs_up_to_residual_z2") is True

    out = {
        "packet_id": "P1967",
        "stage_id": "S917",
        "status": "STRICT_SHANNON_AXIS_SOURCE_INSTANTIATES_P1966_DELTA_SEL_ON_PAIR1_TO_PAIR5__GLOBAL_SELECTOR_CLOSURE_STILL_OPEN",
        "route": "strict_only_shannon_element_order_reference_lane",
        "legacy_bridge_used": False,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "sympy": sp.__version__,
        "repo_grep_evidence": repo_grep_evidence(),
        "input_generated_hashes": {
            "mode_index_assignment_shannon_element_order_reference_strict_core_v1_sha256": hashlib.sha256(
                json.dumps(shannon_assignment, sort_keys=True).encode("utf-8")
            ).hexdigest(),
            "p455_summary_sha256": hashlib.sha256(json.dumps(p455_summary, sort_keys=True).encode("utf-8")).hexdigest(),
            "p1966_sha256": hashlib.sha256(json.dumps(p1966, sort_keys=True).encode("utf-8")).hexdigest(),
        },
        "construction_formula": {
            "ord_Z12": "ord_Z12(0)=1; ord_Z12(x)=12/gcd(x,12) for x != 0",
            "F_2m_ord": "sum_{x=0}^{11} ord_Z12(x) * exp(i*4*pi*m*x/12)",
            "Delta_sel_pair_m": "[[Re(F_2m_ord), Im(F_2m_ord)], [Im(F_2m_ord), -Re(F_2m_ord)]]",
            "gap": "2*abs(F_2m_ord)",
        },
        "pair_rows": pair_rows,
        "machine_checks": {
            "all_pair_gaps_nonzero": all_gaps_nonzero,
            "all_pair_defects_match_exported_F454": all_match_export,
            "p455_shannon_vs_diagonal_axis_alignment_up_to_z2": p455_aligned,
            "p1966_minimal_tensor_shape_present": p1966.get("minimal_symmetry_breaking_premise", {}).get("all_numeric_samples_pass") is True,
            "strict_axis_only_source_pass": bool(all_gaps_nonzero and all_match_export and p455_aligned),
        },
        "theorem_export": {
            "axis_only_positive_statement": "For the exported n=12 strict Shannon element-order reference lane, F_2m(ord_Z12) is nonzero for m=1..5, so the induced Delta_sel_pair_m has nonzero spectral gap and cuts each pair_m O(2) family to residual Z2.",
            "scope": "n=12 QW-2190 Fourier scaffold, pair1..pair5, Shannon element-order reference lane; aligned with diagonal/local lane up to residual Z2 by P455.",
            "not_promoted_to": [
                "sign-sensitive physical orientation datum",
                "admissible S_sel_int",
                "global QW-2191 discharge",
                "F(nadsoliton) => L_SM + L_GR mapping witness",
                "ToE closure",
            ],
        },
        "false_pass_guard": "This packet upgrades P1966 from abstract missing-shape identification to a lane-scoped strict axis-only instantiation, but it does not solve residual Z2 sign or global selector closure.",
        "next_honest_step": "Test whether the downstream observables needed for strict ToE closure are Z2 gauge-irrelevant under the P1967 axes; if not, construct or refute a strict sign-sensitive datum.",
        "lay_explanation": "Repo ma już ścisłe źródło, które wybiera osie w parach trybów: rozkład Shannona oparty o rząd elementu w Z12. P1967 przelicza to na macierze wybierające kierunek osi i pokazuje niezerową szczelinę. Nadal nie wybiera absolutnego zwrotu strzałki, więc pełny globalny selektor pozostaje otwarty.",
    }
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
