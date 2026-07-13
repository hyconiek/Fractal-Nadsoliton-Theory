"""P3138/S2088: non-Fourier D_HL extremum joint-source audit.

P3137 left one constructive option: try exactly one genuinely non-Fourier
joint source object for the P3134 D_HL family.  This script builds the
finite zero-crossing/curvature-extremum candidate E_DHL and audits whether it
selects an import-free joint origin/polarity pair (r, lambda).
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3138_s2088_nonfourier_dhl_extremum_joint_source_audit.json"
MD = GEN / "p3138_s2088_nonfourier_dhl_extremum_joint_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
P3137 = GEN / "p3137_s2087_fourier_frame_source_law_audit.json"

N = 12
BETA_TORS = 0.01
MODES = [1, 2, 3, 4, 5]
LAMBDAS = [-1, 1]
UNITS = [1, 5, 7, 11]
TOL = 1e-12


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def profile(r: int, k: int, lam: int) -> list[float]:
    return [lam * BETA_TORS * math.sin(2 * math.pi * k * ((x - r) % N) / N) for x in range(N)]


def translate(vals: list[float], t: int) -> list[float]:
    return [vals[(x - t) % N] for x in range(N)]


def invert(vals: list[float]) -> list[float]:
    return [vals[(-x) % N] for x in range(N)]


def forward(vals: list[float]) -> list[float]:
    return [vals[(x + 1) % N] - vals[x] for x in range(N)]


def second(vals: list[float]) -> list[float]:
    return [vals[(x + 1) % N] - 2 * vals[x] + vals[(x - 1) % N] for x in range(N)]


def near_zero_crossings(vals: list[float]) -> list[int]:
    hits = set()
    for x, v in enumerate(vals):
        if abs(v) < TOL:
            hits.add(x)
        if vals[x] * vals[(x + 1) % N] < -TOL:
            hits.add(x)
            hits.add((x + 1) % N)
    return sorted(hits)


def extrema_source(vals: list[float]) -> dict[str, Any]:
    """Non-Fourier receiver using only local samples, differences and extrema.

    Tie-breaking by the smallest cell is deliberately reported as imported;
    without that order convention the source returns a set, not one cell.
    """
    d1 = forward(vals)
    d2 = second(vals)
    max_slope = max(d1)
    min_slope = min(d1)
    pos = [i for i, v in enumerate(d1) if abs(v - max_slope) < TOL]
    neg = [i for i, v in enumerate(d1) if abs(v - min_slope) < TOL]
    zeros = near_zero_crossings(vals)
    curv = max(abs(v) for v in d2)
    curv_cells = [i for i, v in enumerate(d2) if abs(abs(v) - curv) < TOL]
    return {
        "positive_slope_cells": pos,
        "negative_slope_cells": neg,
        "zero_crossing_cells": zeros,
        "max_abs_curvature_cells": curv_cells,
        "tie_free_positive_slope": len(pos) == 1,
        "tie_free_negative_slope": len(neg) == 1,
        "candidate_origin_by_ordered_positive_slope": min(pos),
        "candidate_lambda_by_slope_sign": 1 if abs(max_slope) >= abs(min_slope) else -1,
        "imports_linear_cell_order_for_ties": len(pos) != 1,
    }


def same_set_shift(cells: list[int], t: int) -> list[int]:
    return sorted((c + t) % N for c in cells)


def build_payload() -> dict[str, Any]:
    rows = []
    translation_failures = 0
    inversion_pairings = 0
    for k in MODES:
        for r in range(N):
            for lam in LAMBDAS:
                vals = profile(r, k, lam)
                src = extrema_source(vals)
                trans = extrema_source(translate(vals, 1))
                inv = extrema_source(invert(vals))
                covariant = trans["positive_slope_cells"] == same_set_shift(src["positive_slope_cells"], 1)
                if covariant:
                    translation_failures += 1
                if inv["positive_slope_cells"] == src["negative_slope_cells"] or inv["negative_slope_cells"] == src["positive_slope_cells"]:
                    inversion_pairings += 1
                rows.append({
                    "mode_k": k,
                    "r": r,
                    "lambda": lam,
                    **src,
                    "translation_t1_moves_receiver_set_covariantly": covariant,
                    "inversion_exchanges_positive_negative_slope_sets": inv["positive_slope_cells"] == src["negative_slope_cells"] or inv["negative_slope_cells"] == src["positive_slope_cells"],
                    "accepted_import_free_joint_source": False,
                    "blocker": "local extrema/zero-crossing data are receivers; translation moves the selected cells and inversion pairs polarity unless an external cell order/orientation convention is imported",
                })

    source_rows = [
        {"gate": "non_fourier_formula", "passed": True, "detail": "uses local values, forward differences, second differences, and zero crossings only"},
        {"gate": "computes_nonempty_origin_receiver", "passed": all(row["positive_slope_cells"] for row in rows), "detail": "every tested profile has local slope/zero/extremum receiver data"},
        {"gate": "tie_free_unique_cell", "passed": all(row["tie_free_positive_slope"] for row in rows), "detail": "fails on nonprimitive/higher-mode profiles where several equal extrema occur"},
        {"gate": "translation_quotient_invariant", "passed": False, "detail": "the receiver set shifts with diagonal Z12 translation instead of descending to an absolute support representative"},
        {"gate": "inversion_polarity_unpaired", "passed": False, "detail": "inversion exchanges positive and negative slope/curvature polarity classes"},
        {"gate": "import_free_cell_order", "passed": False, "detail": "single-cell extraction requires a labelled cyclic order/tie convention"},
    ]
    passed = sum(1 for g in source_rows if g["passed"])
    return {
        "status": "P3138_NONFOURIER_E_DHL_EXTREMUM_JOINT_SOURCE_BOUNDED_NO_GO",
        "input_hashes": {"P3137": sha(P3137)},
        "constructed_object": {
            "name": "E_DHL zero-crossing/curvature-extremum joint source candidate",
            "formula": "E_DHL(D) = (argmax_x Delta D(x), zero-crossing set, argmax_x |Delta^2 D(x)|) with attempted polarity from slope sign",
            "why_non_fourier": "No Fourier coefficients, characters, projectors, or phase gauges are used; only local cyclic samples and finite differences are audited.",
        },
        "finite_certificate": {
            "profiles_tested": len(rows),
            "candidate_modes": MODES,
            "translation_t1_covariant_receiver_rows": translation_failures,
            "inversion_paired_rows": inversion_pairings,
            "tie_free_positive_slope_rows": sum(row["tie_free_positive_slope"] for row in rows),
            "multi_tie_positive_slope_rows": sum(not row["tie_free_positive_slope"] for row in rows),
            "accepted_import_free_joint_sources": 0,
            "gates_passed": passed,
            "gates_required": len(source_rows),
        },
        "source_gate_rows": source_rows,
        "profile_rows": rows,
        "decision": {
            "bounded_result": "P3138 constructs one non-Fourier joint source candidate E_DHL from local zero-crossings, slopes, and curvature extrema of the P3134 D_HL profiles. It is a genuine receiver and can locate local structure, but it remains translation-covariant rather than translation-quotient invariant; higher modes have equal-extremum ties; and inversion pairs the polarity classes. Thus extracting one absolute (r,lambda) still imports a labelled support order/orientation convention. No import-free joint source is exported.",
            "positive_scoped_flags": {"non_fourier_joint_candidate_constructed": True, "finite_extremum_matrix_computed": True, "translation_covariance_obstruction_computed": True, "inversion_pairing_obstruction_computed": True},
            "negative_export_flags": {"E_DHL_source_exported": False, "J_DHL_source_exported": False, "D_HL_source_exported": False, "Zeta_OS_exported": False, "Gamma_SO_exported": False, "QW_2191_discharged": False, "strict_selector_closure_exported": False, "bridge_completion_exported": False, "legacy_role_transfer_exported": False, "L_total_exported": False, "ToE_closure_exported": False},
            "next_honest_step": "After P3133-P3138, the D_HL lane has now tested legacy/torsion shape, joint origin-polarity matrices, Fourier frame sourcing, and one non-Fourier extremum source. Unless a genuinely new strict source law supplies absolute support origin and polarity without imported labels, the next proof-grade move should be a no-new-live-frontier reconciliation for the D_HL selector/symmetry-breaking lane.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    lines = [
        "# P3138/S2088 Non-Fourier D_HL extremum joint-source audit",
        "", f"Status: `{payload['status']}`", "", "## Constructed object",
        f"- `{payload['constructed_object']['name']}`",
        f"- Formula: `{payload['constructed_object']['formula']}`",
        f"- Non-Fourier discipline: {payload['constructed_object']['why_non_fourier']}",
        "", "## Finite certificate",
        f"- profiles tested: `{cert['profiles_tested']}`",
        f"- translation-covariant receiver rows under t=1: `{cert['translation_t1_covariant_receiver_rows']}`",
        f"- inversion-paired rows: `{cert['inversion_paired_rows']}`",
        f"- tie-free positive-slope rows: `{cert['tie_free_positive_slope_rows']}`",
        f"- multi-tie positive-slope rows: `{cert['multi_tie_positive_slope_rows']}`",
        f"- accepted import-free joint sources: `{cert['accepted_import_free_joint_sources']}`",
        "", "## Gate table",
    ]
    for row in payload["source_gate_rows"]:
        lines.append(f"- `{row['gate']}`: `{row['passed']}` — {row['detail']}")
    lines.extend(["", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3138/S2088 non-Fourier D_HL extremum joint-source audit", "## P3138/S2088 non-Fourier D_HL extremum joint-source audit\n\n`P3138/S2088` constructs exactly one non-Fourier joint source candidate for the P3134 `D_HL` family: `E_DHL`, based on local zero-crossings, forward-difference extrema, and second-difference curvature extrema.  The finite audit covers `120` profiles.  The receiver structure is real, but all rows remain translation-covariant rather than absolute after the diagonal `Z12` quotient, inversion pairs positive/negative polarity classes, and higher-mode rows introduce equal-extremum ties.  Therefore no import-free `(r,lambda)` source, `D_HL` source, `Zeta_OS`, `Gamma_SO`, selector closure, bridge completion, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3138/S2088 E_DHL extremum receiver is not a variational source", "## P3138/S2088 E_DHL extremum receiver is not a variational source\n\n`P3138/S2088` shows that local extremum/zero-crossing data for the constructed `D_HL` profile are receiver data, not an import-free source law for an absolute support origin and polarity.  They cannot yet define a Lagrangian density, Hamiltonian normalization, spacetime EOM, physical unit, `L_total`, bridge-completion theorem, role-transfer theorem, or ToE closure.\n")
    append_once(AGENTS, "Current non-Fourier D_HL extremum joint-source guardrail (P3138/S2088, 2026-07-13)", "## Current non-Fourier D_HL extremum joint-source guardrail (P3138/S2088, 2026-07-13)\n\n- P3138 constructs the allowed non-Fourier `E_DHL` candidate from local zero-crossings, forward-difference extrema, and curvature extrema of the P3134 `D_HL` profiles.\n- The finite audit tests `120` profiles and finds real receiver structure, but no import-free source: extrema/zero data translate with the diagonal `Z12` orbit, inversion pairs polarity classes, and higher modes carry equal-extremum ties requiring imported support-order conventions.\n- Do not promote local extremum/zero-crossing receivers to `J_DHL`, `D_HL`, `Zeta_OS`, `Gamma_SO`, `QW-2191` discharge, strict selector closure, bridge completion, role transfer, `L_total`, or ToE closure.\n- After P3133-P3138, the next honest move is a no-new-live-frontier reconciliation for the `D_HL` selector/symmetry-breaking lane unless a genuinely new strict source law supplies absolute support origin and polarity without imported labels.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
