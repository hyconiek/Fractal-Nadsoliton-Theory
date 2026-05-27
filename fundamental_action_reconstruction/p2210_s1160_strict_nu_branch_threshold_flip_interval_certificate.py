#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2209 = GEN / "p2209_s1159_strict_nu_branch_transport_nonvanishing_threshold_map.json"
OUT = GEN / "p2210_s1160_strict_nu_branch_threshold_flip_interval_certificate.json"
MD = GEN / "p2210_s1160_strict_nu_branch_threshold_flip_interval_certificate.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2209 = load(IN_2209)
    block = p2209.get("strict_nu_branch_transport_nonvanishing_threshold_map", {}) or {}
    inputs = block.get("inputs", {}) or {}

    min_margin = float(inputs.get("min_majorant_margin", 0.0) or 0.0)
    c1 = float(inputs.get("linear_lower_bound_coeff_c1", 0.0) or 0.0)
    c2 = float(inputs.get("quadratic_lower_bound_coeff_c2", 0.0) or 0.0)

    if not (min_margin > 0.0 and c1 > 0.0 and c2 > 0.0):
        raise RuntimeError("P2209 inputs invalid for threshold flip interval certification")

    base_rows = block.get("threshold_map_rows", []) or []
    certified = [r for r in base_rows if bool(r.get("certifies_margin_exceedance", False))]
    uncertified = [r for r in base_rows if not bool(r.get("certifies_margin_exceedance", False))]

    has_flip_in_coarse_window = bool(certified) and bool(uncertified)
    max_scanned_dm = max(float(r["abs_dm"]) for r in base_rows) if base_rows else 0.0
    max_scanned_ratio = max(float(r.get("effective_over_min_margin", 0.0)) for r in base_rows) if base_rows else 0.0

    linear_needed_dm = min_margin / c1
    quadratic_needed_dm = (min_margin / c2) ** 0.5
    conservative_needed_dm = min(linear_needed_dm, quadratic_needed_dm)

    if has_flip_in_coarse_window:
        first_cert = min(certified, key=lambda r: float(r["abs_dm"]))
        last_uncert = max(uncertified, key=lambda r: float(r["abs_dm"]))
        flip_interval = {
            "abs_dm_lo_uncertified": float(last_uncert["abs_dm"]),
            "abs_dm_hi_certified": float(first_cert["abs_dm"]),
        }
    else:
        flip_interval = None

    payload = {
        "schema_version": "p2210_s1160_v1",
        "packet_id": "P2210",
        "stage_id": "S1160",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_THRESHOLD_FLIP_INTERVAL_CERTIFICATE",
        "strict_nu_branch_threshold_flip_interval_certificate": {
            "certificate_id": "STRICT_NU_BRANCH_THRESHOLD_FLIP_INTERVAL_CERTIFICATE_V1",
            "source_packet": str(IN_2209.relative_to(ROOT)),
            "has_flip_in_coarse_window": has_flip_in_coarse_window,
            "coarse_flip_interval": flip_interval,
            "coarse_window_summary": {
                "max_scanned_abs_dm": max_scanned_dm,
                "max_scanned_effective_over_min_margin": max_scanned_ratio,
                "all_rows_uncertified": len(certified) == 0,
            },
            "implied_needed_abs_dm_estimates": {
                "linear_needed_abs_dm": linear_needed_dm,
                "quadratic_needed_abs_dm": quadratic_needed_dm,
                "conservative_needed_abs_dm": conservative_needed_dm,
                "interpretation": "minimum modeled |m-1| scale where lower bound can reach min_margin",
            },
            "theorem_scope_limit": "finite-grid compact-lane certificate: either flip interval if observed or explicit no-flip-in-window evidence; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2211_candidate",
            "goal": "extend |m-1| scan up to conservative_needed_abs_dm and certify first exceedance interval",
        },
        "gatekeeper_checks": {
            "flip_interval_certificate_exported": True,
            "coarse_rows_present": len(base_rows) > 0,
            "no_flip_detected_in_coarse_window": len(certified) == 0,
            "needed_abs_dm_exceeds_current_window": conservative_needed_dm > max_scanned_dm,
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
                "# P2210 S1160: strict nu-branch threshold flip interval certificate",
                "",
                f"- has flip in coarse P2209 window: `{has_flip_in_coarse_window}`",
                f"- max scanned |m-1|: `{max_scanned_dm:.12e}`",
                f"- max effective_over_min_margin on scan: `{max_scanned_ratio:.12e}`",
                f"- conservative needed |m-1| estimate: `{conservative_needed_dm:.12e}`",
                "",
                "No global closure claim; compact-lane finite-window certificate only.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
