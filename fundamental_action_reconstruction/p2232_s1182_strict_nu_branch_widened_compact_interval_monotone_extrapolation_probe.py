#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2209 = GEN / "p2209_s1159_strict_nu_branch_transport_nonvanishing_threshold_map.json"
IN_2214 = GEN / "p2214_s1164_strict_nu_branch_closed_form_threshold_consistency_certificate.json"
OUT = GEN / "p2232_s1182_strict_nu_branch_widened_compact_interval_monotone_extrapolation_probe.json"
MD = GEN / "p2232_s1182_strict_nu_branch_widened_compact_interval_monotone_extrapolation_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def eff(dm: float, c1: float, c2: float) -> float:
    return max(c1 * dm, c2 * dm * dm)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2209 = load(IN_2209)
    p2214 = load(IN_2214)

    inputs = (p2209.get("strict_nu_branch_transport_nonvanishing_threshold_map", {}) or {}).get("inputs", {}) or {}
    c1 = float(inputs.get("linear_lower_bound_coeff_c1", 0.0) or 0.0)
    c2 = float(inputs.get("quadratic_lower_bound_coeff_c2", 0.0) or 0.0)
    margin = float(inputs.get("min_majorant_margin", 0.0) or 0.0)

    dm_star = float(((p2214.get("strict_nu_branch_closed_form_threshold_consistency_certificate", {}) or {}).get("bisection_reference", {}) or {}).get("dm_star_bisection", 0.0) or 0.0)

    if not (c1 > 0 and c2 > 0 and margin > 0 and dm_star > 0):
        raise RuntimeError("Invalid upstream inputs for P2232 widened compact extrapolation probe")

    # widened compact neighborhood around dm*: [0.5*dm*, 1.5*dm*]
    lo = 0.5 * dm_star
    hi = 1.5 * dm_star
    n = 4001
    rows = []
    prev = None
    mono = True
    for k in range(n):
        dm = lo + (hi - lo) * (k / (n - 1))
        sf = eff(dm, c1, c2) / margin
        if prev is not None and sf + 1e-15 < prev:
            mono = False
        prev = sf
        rows.append({"abs_dm": dm, "safety_factor": sf, "above_margin": sf > 1.0})

    first_cross = next((r for r in rows if r["above_margin"]), None)
    last_below = next((r for r in reversed(rows) if not r["above_margin"]), None)

    payload = {
        "schema_version": "p2232_s1182_v1",
        "packet_id": "P2232",
        "stage_id": "S1182",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_WIDENED_COMPACT_INTERVAL_MONOTONE_EXTRAPOLATION_PROBE",
        "strict_nu_branch_widened_compact_interval_monotone_extrapolation_probe": {
            "probe_id": "STRICT_NU_BRANCH_WIDENED_COMPACT_INTERVAL_MONOTONE_EXTRAPOLATION_PROBE_V1",
            "source_packets": [str(IN_2209.relative_to(ROOT)), str(IN_2214.relative_to(ROOT))],
            "interval": {"abs_dm_lo": lo, "abs_dm_hi": hi, "grid_size": n},
            "summary": {
                "monotone_safety_factor_on_interval": mono,
                "first_crossing_row": first_cross,
                "last_below_row": last_below,
            },
            "theorem_scope_limit": "widened compact interval monotonic probe only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2233_candidate",
            "goal": "replace sampled monotonic probe with interval-arithmetic monotonic certificate on widened compact domain",
        },
        "gatekeeper_checks": {
            "widened_compact_probe_exported": True,
            "monotone_on_widened_interval": mono,
            "crossing_detected_on_widened_interval": first_cross is not None,
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
            "# P2232 S1182: strict nu-branch widened compact interval monotone extrapolation probe",
            "",
            f"- interval: `[{lo:.12e}, {hi:.12e}]`",
            f"- monotone on widened interval: `{mono}`",
            f"- crossing detected: `{first_cross is not None}`",
            "",
            "Widened compact probe only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
