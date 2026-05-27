#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2203 = GEN / "p2203_s1153_strict_frw_bianchi_transport_residual_map_under_shared_majorant.json"
IN_2208 = GEN / "p2208_s1158_strict_nu_branch_separation_lower_bound_certificate.json"
OUT = GEN / "p2209_s1159_strict_nu_branch_transport_nonvanishing_threshold_map.json"
MD = GEN / "p2209_s1159_strict_nu_branch_transport_nonvanishing_threshold_map.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2203 = load(IN_2203)
    p2208 = load(IN_2208)

    residual_block = p2203.get("strict_frw_bianchi_transport_residual_map_under_shared_majorant", {}) or {}
    summary = residual_block.get("summary", {}) or {}
    min_margin = float(summary.get("min_majorant_margin_from_p2201", 0.0) or 0.0)

    cert_block = p2208.get("strict_nu_branch_separation_lower_bound_certificate", {}) or {}
    c1 = float(cert_block.get("min_linear_coefficient_c1", 0.0) or 0.0)
    c2 = float(cert_block.get("min_quadratic_coefficient_c2", 0.0) or 0.0)
    eps = float(cert_block.get("scan_exclusion_radius", 0.0) or 0.0)

    dm_grid = [0.0, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2]
    rows = []
    for dm in dm_grid:
        lower_linear = c1 * dm
        lower_quad = c2 * dm * dm
        lower_bound = max(lower_linear, lower_quad)
        exceeds_margin = (min_margin > 0.0) and (lower_bound > min_margin)
        rows.append(
            {
                "abs_dm": dm,
                "lower_bound_linear": lower_linear,
                "lower_bound_quadratic": lower_quad,
                "effective_lower_bound": lower_bound,
                "effective_over_min_margin": (lower_bound / min_margin) if min_margin > 0.0 else float("inf"),
                "certifies_margin_exceedance": exceeds_margin,
            }
        )

    if min_margin > 0.0 and c1 > 0.0:
        dm_threshold_linear = min_margin / c1
    else:
        dm_threshold_linear = float("inf")

    if min_margin > 0.0 and c2 > 0.0:
        dm_threshold_quadratic = math.sqrt(min_margin / c2)
    else:
        dm_threshold_quadratic = float("inf")

    first_certified = next((r for r in rows if r["certifies_margin_exceedance"]), None)

    payload = {
        "schema_version": "p2209_s1159_v1",
        "packet_id": "P2209",
        "stage_id": "S1159",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_TRANSPORT_NONVANISHING_THRESHOLD_MAP",
        "strict_nu_branch_transport_nonvanishing_threshold_map": {
            "map_id": "STRICT_NU_BRANCH_TRANSPORT_NONVANISHING_THRESHOLD_MAP_V1",
            "source_packets": [
                str(IN_2203.relative_to(ROOT)),
                str(IN_2208.relative_to(ROOT)),
            ],
            "inputs": {
                "min_majorant_margin": min_margin,
                "linear_lower_bound_coeff_c1": c1,
                "quadratic_lower_bound_coeff_c2": c2,
                "scan_exclusion_radius_eps": eps,
            },
            "threshold_estimates": {
                "dm_threshold_from_linear_bound": dm_threshold_linear,
                "dm_threshold_from_quadratic_bound": dm_threshold_quadratic,
                "interpretation": "if |m-1| is above threshold, certified lower bound exceeds min_margin on scanned compact lane",
            },
            "threshold_map_rows": rows,
            "first_certified_row": first_certified,
            "theorem_scope_limit": "certificate links compact-scan lower bounds and shared-majorant margin; it is not a global all-background Task-3 closure theorem",
        },
        "recommended_next_honest_step": {
            "id": "P2210_candidate",
            "goal": "refine near-threshold dm grid adaptively and export a narrow interval certificate where margin exceedance flips from unproven to certified",
        },
        "gatekeeper_checks": {
            "threshold_map_exported": True,
            "positive_linear_coeff": c1 > 0.0,
            "positive_quadratic_coeff": c2 > 0.0,
            "first_certified_row_present": first_certified is not None,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2209 S1159: strict nu-branch transport non-vanishing threshold map",
                "",
                f"- c1: `{c1:.12e}`",
                f"- c2: `{c2:.12e}`",
                f"- min margin (from P2201 via P2203): `{min_margin:.12e}`",
                f"- linear threshold |m-1|: `{dm_threshold_linear:.12e}`",
                f"- quadratic threshold |m-1|: `{dm_threshold_quadratic:.12e}`",
                f"- first certified row present: `{first_certified is not None}`",
                "",
                "Threshold map is compact-scan scoped; no global Task-3 closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
