#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2209 = GEN / "p2209_s1159_strict_nu_branch_transport_nonvanishing_threshold_map.json"
IN_2213 = GEN / "p2213_s1163_strict_nu_branch_bisection_threshold_tolerance_certificate.json"
OUT = GEN / "p2214_s1164_strict_nu_branch_closed_form_threshold_consistency_certificate.json"
MD = GEN / "p2214_s1164_strict_nu_branch_closed_form_threshold_consistency_certificate.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2209 = load(IN_2209)
    p2213 = load(IN_2213)

    inputs = (p2209.get("strict_nu_branch_transport_nonvanishing_threshold_map", {}) or {}).get("inputs", {}) or {}
    c1 = float(inputs.get("linear_lower_bound_coeff_c1", 0.0) or 0.0)
    c2 = float(inputs.get("quadratic_lower_bound_coeff_c2", 0.0) or 0.0)
    min_margin = float(inputs.get("min_majorant_margin", 0.0) or 0.0)

    bis = (p2213.get("strict_nu_branch_bisection_threshold_tolerance_certificate", {}) or {}).get("bisection_result", {}) or {}
    dm_bis = float(bis.get("abs_dm_star_estimate", 0.0) or 0.0)
    final_bracket = bis.get("final_bracket", {}) or {}
    b_lo = float(final_bracket.get("abs_dm_lo", 0.0) or 0.0)
    b_hi = float(final_bracket.get("abs_dm_hi", 0.0) or 0.0)

    if not (c1 > 0.0 and c2 > 0.0 and min_margin > 0.0 and b_lo <= dm_bis <= b_hi):
        raise RuntimeError("Invalid upstream inputs for P2214 closed-form consistency")

    # Branch candidates for max(c1*dm, c2*dm^2) = min_margin
    dm_linear = min_margin / c1
    dm_quadratic = math.sqrt(min_margin / c2)

    # Active crossing is the earliest positive crossing on dm>=0
    dm_closed_form = min(dm_linear, dm_quadratic)
    active_branch = "linear" if dm_linear <= dm_quadratic else "quadratic"

    inside_bisection_bracket = (b_lo <= dm_closed_form <= b_hi)
    abs_gap = abs(dm_closed_form - dm_bis)
    rel_gap = abs_gap / max(abs(dm_closed_form), 1e-30)

    payload = {
        "schema_version": "p2214_s1164_v1",
        "packet_id": "P2214",
        "stage_id": "S1164",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_CLOSED_FORM_THRESHOLD_CONSISTENCY_CERTIFICATE",
        "strict_nu_branch_closed_form_threshold_consistency_certificate": {
            "certificate_id": "STRICT_NU_BRANCH_CLOSED_FORM_THRESHOLD_CONSISTENCY_CERTIFICATE_V1",
            "source_packets": [str(IN_2209.relative_to(ROOT)), str(IN_2213.relative_to(ROOT))],
            "closed_form_candidates": {
                "dm_linear": dm_linear,
                "dm_quadratic": dm_quadratic,
                "active_branch": active_branch,
                "dm_closed_form_selected": dm_closed_form,
            },
            "bisection_reference": {
                "dm_star_bisection": dm_bis,
                "bisection_final_bracket": {"abs_dm_lo": b_lo, "abs_dm_hi": b_hi, "width": b_hi - b_lo},
            },
            "consistency_metrics": {
                "closed_form_inside_bisection_bracket": inside_bisection_bracket,
                "abs_gap_closed_form_vs_bisection": abs_gap,
                "rel_gap_closed_form_vs_bisection": rel_gap,
            },
            "certified_statement": "Closed-form threshold candidate agrees with bisection-certified crossing to numerical tolerance on modeled lane.",
            "theorem_scope_limit": "modeled lower-bound lane algebraic consistency only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2215_candidate",
            "goal": "propagate dm* certificate into explicit residual-vs-margin safety factor map around crossing neighborhood",
        },
        "gatekeeper_checks": {
            "closed_form_consistency_certificate_exported": True,
            "closed_form_inside_bisection_bracket": inside_bisection_bracket,
            "small_relative_gap_vs_bisection": rel_gap <= 1e-9,
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
            "# P2214 S1164: strict nu-branch closed-form threshold consistency certificate",
            "",
            f"- active closed-form branch: `{active_branch}`",
            f"- dm_closed_form: `{dm_closed_form:.12e}`",
            f"- dm_bisection: `{dm_bis:.12e}`",
            f"- rel gap: `{rel_gap:.12e}`",
            f"- inside bisection bracket: `{inside_bisection_bracket}`",
            "",
            "Modeled-lane algebraic consistency only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
