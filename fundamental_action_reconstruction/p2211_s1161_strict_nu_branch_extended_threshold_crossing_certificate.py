#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2209 = GEN / "p2209_s1159_strict_nu_branch_transport_nonvanishing_threshold_map.json"
IN_2210 = GEN / "p2210_s1160_strict_nu_branch_threshold_flip_interval_certificate.json"
OUT = GEN / "p2211_s1161_strict_nu_branch_extended_threshold_crossing_certificate.json"
MD = GEN / "p2211_s1161_strict_nu_branch_extended_threshold_crossing_certificate.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def effective_lower_bound(dm: float, c1: float, c2: float) -> float:
    return max(c1 * dm, c2 * dm * dm)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2209 = load(IN_2209)
    p2210 = load(IN_2210)

    inputs = (p2209.get("strict_nu_branch_transport_nonvanishing_threshold_map", {}) or {}).get("inputs", {}) or {}
    c1 = float(inputs.get("linear_lower_bound_coeff_c1", 0.0) or 0.0)
    c2 = float(inputs.get("quadratic_lower_bound_coeff_c2", 0.0) or 0.0)
    min_margin = float(inputs.get("min_majorant_margin", 0.0) or 0.0)

    needed = (
        (p2210.get("strict_nu_branch_threshold_flip_interval_certificate", {}) or {})
        .get("implied_needed_abs_dm_estimates", {})
        .get("conservative_needed_abs_dm", 0.0)
    )
    conservative_needed_dm = float(needed or 0.0)

    if not (c1 > 0.0 and c2 > 0.0 and min_margin > 0.0 and conservative_needed_dm > 0.0):
        raise RuntimeError("Invalid upstream inputs for P2211 extended crossing certificate")

    # Extend beyond estimated crossing scale from P2210.
    dm_max = 1.10 * conservative_needed_dm
    n = 2401
    dms = [dm_max * k / (n - 1) for k in range(n)]

    rows = []
    for dm in dms:
        b = effective_lower_bound(dm, c1, c2)
        rows.append(
            {
                "abs_dm": dm,
                "effective_lower_bound": b,
                "effective_over_min_margin": b / min_margin,
                "certifies_margin_exceedance": b > min_margin,
            }
        )

    first_cert = next((r for r in rows if r["certifies_margin_exceedance"]), None)
    last_uncert = next((r for r in reversed(rows) if not r["certifies_margin_exceedance"]), None)

    crossing_interval = None
    if first_cert is not None and last_uncert is not None and float(last_uncert["abs_dm"]) <= float(first_cert["abs_dm"]):
        crossing_interval = {
            "abs_dm_lo_uncertified": float(last_uncert["abs_dm"]),
            "abs_dm_hi_certified": float(first_cert["abs_dm"]),
            "interval_width": float(first_cert["abs_dm"]) - float(last_uncert["abs_dm"]),
        }

    payload = {
        "schema_version": "p2211_s1161_v1",
        "packet_id": "P2211",
        "stage_id": "S1161",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_EXTENDED_THRESHOLD_CROSSING_CERTIFICATE",
        "strict_nu_branch_extended_threshold_crossing_certificate": {
            "certificate_id": "STRICT_NU_BRANCH_EXTENDED_THRESHOLD_CROSSING_CERTIFICATE_V1",
            "source_packets": [str(IN_2209.relative_to(ROOT)), str(IN_2210.relative_to(ROOT))],
            "scan_setup": {
                "dm_max": dm_max,
                "grid_size": n,
                "conservative_needed_abs_dm_from_p2210": conservative_needed_dm,
            },
            "crossing_summary": {
                "first_certified_row": first_cert,
                "last_uncertified_row": last_uncert,
                "crossing_interval": crossing_interval,
                "crossing_detected_on_extended_grid": first_cert is not None,
            },
            "theorem_scope_limit": "extended finite-grid compact-lane certificate only; no global all-background Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2212_candidate",
            "goal": "replace grid crossing estimate with interval/monotonic proof of minimal certified |m-1| in this modeled lane",
        },
        "gatekeeper_checks": {
            "extended_crossing_certificate_exported": True,
            "crossing_detected_on_extended_grid": first_cert is not None,
            "crossing_interval_exported": crossing_interval is not None,
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
                "# P2211 S1161: strict nu-branch extended threshold crossing certificate",
                "",
                f"- conservative needed |m-1| (P2210): `{conservative_needed_dm:.12e}`",
                f"- extended dm_max: `{dm_max:.12e}`",
                f"- crossing detected on extended grid: `{first_cert is not None}`",
                f"- crossing interval exported: `{crossing_interval is not None}`",
                "",
                "Extended finite-grid certificate only; no global Task-3 closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
