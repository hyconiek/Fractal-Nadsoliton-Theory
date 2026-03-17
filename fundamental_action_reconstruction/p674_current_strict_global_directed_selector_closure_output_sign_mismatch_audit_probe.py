#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_DIRECTED_STATE = GENERATED / "selector_state_global_c_v1_directed_strict_v1.json"
IN_GLOBAL_OUTPUT_OPERATOR = (
    GENERATED / "selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json"
)

OUT_JSON = (
    GENERATED
    / "p674_current_strict_global_directed_selector_closure_output_sign_mismatch_audit_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p674_current_strict_global_directed_selector_closure_output_sign_mismatch_audit_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def dot(a: list[float], b: list[float]) -> float:
    return float(sum(float(x) * float(y) for x, y in zip(a, b, strict=True)))


def fourier_c_s_basis_vectors_on_z12(k: int) -> tuple[list[float], list[float]]:
    if k <= 0 or k >= 6:
        raise ValueError("This probe expects k in {1,2,3,4,5} for the Z_12 Fourier-degenerate pairs.")
    n = 12
    norm = math.sqrt(2.0 / n)  # k in {1..5} uses the standard sqrt(2/n) normalization
    c = [norm * math.cos(2.0 * math.pi * k * x / n) for x in range(n)]
    s = [norm * math.sin(2.0 * math.pi * k * x / n) for x in range(n)]
    return c, s


def sign_of(x: float, tol: float) -> int:
    if abs(x) <= tol:
        return 0
    return 1 if x > 0.0 else -1


def is_2x2_matrix(v: Any) -> bool:
    return (
        isinstance(v, list)
        and len(v) == 2
        and all(isinstance(row, list) and len(row) == 2 for row in v)
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_DIRECTED_STATE.exists() or not IN_GLOBAL_OUTPUT_OPERATOR.exists():
        status = "P674_NOT_COMPUTABLE_MISSING_REQUIRED_INPUTS"
        missing = [
            str(p.relative_to(REPO))
            for p in (IN_DIRECTED_STATE, IN_GLOBAL_OUTPUT_OPERATOR)
            if not p.exists()
        ]
        artifact = {
            "stage": "P674",
            "lane": "current_strict_global_directed_selector_closure_output_sign_mismatch_audit_probe_only",
            "status": status,
            "missing": missing,
            "no_false_pass": True,
        }
        summary = {
            "stage": "P674",
            "lane": artifact["lane"],
            "status": status,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    state = load_json(IN_DIRECTED_STATE)
    outop = load_json(IN_GLOBAL_OUTPUT_OPERATOR)

    lane = "current_strict_global_directed_selector_closure_output_sign_mismatch_audit_probe_only"
    tol = 1e-12

    u_vectors = (((state.get("outputs") or {}).get("u_vectors_directed")) or {}) if isinstance(state, dict) else {}
    charts = (outop.get("charts") or {}) if isinstance(outop, dict) else {}

    per_pair: dict[str, Any] = {}
    signs: dict[str, int] = {}
    blocking_mismatches: list[str] = []

    expected_pairs = [f"pair{i}" for i in range(1, 6)]
    for i, pair in enumerate(expected_pairs, start=1):
        u_key = f"u_{i}"
        u = u_vectors.get(u_key)
        if not (isinstance(u, list) and len(u) == 12 and all(isinstance(x, (int, float)) for x in u)):
            blocking_mismatches.append(f"missing_or_invalid_{u_key}")
            continue

        chart = charts.get(pair) or {}
        y = chart.get("Y_sel_matrix_in_c_s_to_o")
        if not is_2x2_matrix(y):
            blocking_mismatches.append(f"missing_or_invalid_Y_sel_{pair}")
            continue

        c_vec, s_vec = fourier_c_s_basis_vectors_on_z12(i)
        c_coord = dot(c_vec, [float(x) for x in u])
        s_coord = dot(s_vec, [float(x) for x in u])

        o_plus = float(y[0][0]) * c_coord + float(y[0][1]) * s_coord
        o_minus = float(y[1][0]) * c_coord + float(y[1][1]) * s_coord
        sgn = sign_of(o_plus, tol=tol)

        per_pair[pair] = {
            "pair": pair,
            "k": i,
            "u_key": u_key,
            "coords_in_c_s": [c_coord, s_coord],
            "Y_sel_matrix_in_c_s_to_o": y,
            "output_coords_in_o_plus_o_minus": [o_plus, o_minus],
            "o_plus_sign": sgn,
        }
        signs[pair] = sgn

    nonzero_signs = [v for v in signs.values() if v != 0]
    unique_signs = sorted(set(nonzero_signs))
    all_same_sign = len(unique_signs) == 1 and len(nonzero_signs) == 5

    # A single global sign flip multiplies all charts by the same factor; it cannot reconcile mixed signs.
    global_sign_fix_possible = all_same_sign

    if blocking_mismatches:
        status = "P674_REQUIRES_REVIEW_INSUFFICIENT_OR_INVALID_INPUT_EXPORTS"
    else:
        status = (
            "PASS_AUDIT_DIRECTED_CLOSURE_OUTPUT_SIGN_MISMATCH_ACROSS_CHARTS"
            if not all_same_sign
            else "PASS_AUDIT_DIRECTED_CLOSURE_OUTPUT_SIGN_CONSISTENT_ACROSS_CHARTS"
        )

    meaning = (
        "Under the current exported directed selector state (premise-based) and the exported fixed global output-channel bases, "
        "the induced directed output vector is not chart-independent (o_+ sign differs across charts). "
        "Therefore no global directed closure outcome can be promoted from these exports without adding an explicit sign-lift/section choice; "
        "projective (ray) closure remains well-defined."
        if (not blocking_mismatches and not all_same_sign)
        else (
            "Under the current exported directed selector state and exported fixed output-channel bases, the induced directed output vector is chart-consistent. "
            "This probe does not claim uniqueness beyond this sign check and does not imply strict-core selector closure."
        )
    )

    artifact = {
        "stage": "P674",
        "lane": lane,
        "status": status,
        "inputs": {
            "directed_state": str(IN_DIRECTED_STATE.relative_to(REPO)),
            "global_output_operator": str(IN_GLOBAL_OUTPUT_OPERATOR.relative_to(REPO)),
        },
        "tolerance": tol,
        "per_pair": per_pair,
        "o_plus_signs_by_pair": signs,
        "sign_consistency": {
            "nonzero_signs": nonzero_signs,
            "unique_nonzero_signs": unique_signs,
            "all_same_sign": all_same_sign,
            "global_sign_fix_possible": global_sign_fix_possible,
        },
        "meaning": meaning,
        "blocking_mismatches": blocking_mismatches,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P674",
        "lane": lane,
        "status": status,
        "all_same_sign": all_same_sign,
        "global_sign_fix_possible": global_sign_fix_possible,
        "unique_nonzero_signs": unique_signs,
        "blocking_mismatches": blocking_mismatches,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

