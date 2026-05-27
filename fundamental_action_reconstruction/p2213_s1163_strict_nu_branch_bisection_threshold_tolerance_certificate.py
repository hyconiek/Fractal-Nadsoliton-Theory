#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2209 = GEN / "p2209_s1159_strict_nu_branch_transport_nonvanishing_threshold_map.json"
IN_2212 = GEN / "p2212_s1162_strict_nu_branch_monotonic_threshold_bracket_certificate.json"
OUT = GEN / "p2213_s1163_strict_nu_branch_bisection_threshold_tolerance_certificate.json"
MD = GEN / "p2213_s1163_strict_nu_branch_bisection_threshold_tolerance_certificate.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def effective_bound(dm: float, c1: float, c2: float) -> float:
    return max(c1 * dm, c2 * dm * dm)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2209 = load(IN_2209)
    p2212 = load(IN_2212)

    inputs = (p2209.get("strict_nu_branch_transport_nonvanishing_threshold_map", {}) or {}).get("inputs", {}) or {}
    c1 = float(inputs.get("linear_lower_bound_coeff_c1", 0.0) or 0.0)
    c2 = float(inputs.get("quadratic_lower_bound_coeff_c2", 0.0) or 0.0)
    min_margin = float(inputs.get("min_majorant_margin", 0.0) or 0.0)

    cert = (p2212.get("strict_nu_branch_monotonic_threshold_bracket_certificate", {}) or {})
    bracket = cert.get("input_bracket", {}) or {}
    lo = float(bracket.get("abs_dm_lo", 0.0) or 0.0)
    hi = float(bracket.get("abs_dm_hi", 0.0) or 0.0)

    if not (c1 > 0.0 and c2 > 0.0 and min_margin > 0.0 and 0.0 <= lo < hi):
        raise RuntimeError("Invalid upstream inputs for P2213 bisection tolerance certificate")

    def f(dm: float) -> float:
        return effective_bound(dm, c1, c2) - min_margin

    f_lo = f(lo)
    f_hi = f(hi)
    if not (f_lo <= 0.0 and f_hi > 0.0):
        raise RuntimeError("Input bracket does not straddle threshold crossing")

    tol = 1e-12
    max_iter = 200
    trace = []
    a, b = lo, hi
    fa, fb = f_lo, f_hi

    for i in range(max_iter):
        m = 0.5 * (a + b)
        fm = f(m)
        trace.append({"iter": i, "a": a, "b": b, "m": m, "f_m": fm, "width": b - a})
        if (b - a) <= tol:
            break
        if fm <= 0.0:
            a, fa = m, fm
        else:
            b, fb = m, fm

    dm_star_est = 0.5 * (a + b)
    width = b - a

    payload = {
        "schema_version": "p2213_s1163_v1",
        "packet_id": "P2213",
        "stage_id": "S1163",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_BISECTION_THRESHOLD_TOLERANCE_CERTIFICATE",
        "strict_nu_branch_bisection_threshold_tolerance_certificate": {
            "certificate_id": "STRICT_NU_BRANCH_BISECTION_THRESHOLD_TOLERANCE_CERTIFICATE_V1",
            "source_packets": [str(IN_2209.relative_to(ROOT)), str(IN_2212.relative_to(ROOT))],
            "input_bracket": {"abs_dm_lo": lo, "abs_dm_hi": hi, "f_lo": f_lo, "f_hi": f_hi},
            "bisection_result": {
                "abs_dm_star_estimate": dm_star_est,
                "final_bracket": {"abs_dm_lo": a, "abs_dm_hi": b, "width": width},
                "tolerance_target": tol,
                "iterations_used": len(trace),
                "trace_tail": trace[-8:],
            },
            "certified_statement": "Threshold crossing dm* is bracketed within final interval; effective lower bound crosses min_margin within certified tolerance on modeled lane.",
            "theorem_scope_limit": "modeled lower-bound lane bisection certificate only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2214_candidate",
            "goal": "derive closed-form dm* candidate per active branch and compare with bisection bracket for algebraic consistency",
        },
        "gatekeeper_checks": {
            "bisection_certificate_exported": True,
            "input_bracket_straddles_crossing": (f_lo <= 0.0 and f_hi > 0.0),
            "final_width_within_tolerance": width <= tol,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2213 S1163: strict nu-branch bisection threshold tolerance certificate",
            "",
            f"- input bracket: `[{lo:.12e}, {hi:.12e}]`",
            f"- dm* estimate: `{dm_star_est:.12e}`",
            f"- final width: `{width:.12e}`",
            f"- tolerance target met: `{width <= tol}`",
            "",
            "Modeled-lane bisection certificate only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
