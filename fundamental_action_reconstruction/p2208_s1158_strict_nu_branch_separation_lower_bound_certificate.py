#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2207 = GEN / "p2207_s1157_strict_nu_branch_interval_certified_mismatch_bound.json"
OUT = GEN / "p2208_s1158_strict_nu_branch_separation_lower_bound_certificate.json"
MD = GEN / "p2208_s1158_strict_nu_branch_separation_lower_bound_certificate.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2207 = load(IN_2207)
    cert = p2207.get("strict_nu_branch_interval_certified_mismatch_bound", {}) or {}
    table = cert.get("sup_norm_table", []) or []
    eps = float(cert.get("m_scan", {}).get("exclude_radius_around_1", 5e-4) or 5e-4)

    ratios = []
    for row in table:
        m = float(row["m"])
        s = float(row["sup_norm_diff_on_compact_d"])
        dm = abs(m - 1.0)
        if dm <= 0.0:
            continue
        ratios.append({
            "m": m,
            "abs_dm": dm,
            "sup_norm": s,
            "ratio_sup_over_dm": s / dm,
            "ratio_sup_over_dm2": s / (dm * dm),
        })

    c1 = min(r["ratio_sup_over_dm"] for r in ratios)
    c2 = min(r["ratio_sup_over_dm2"] for r in ratios)
    worst = min(ratios, key=lambda r: r["ratio_sup_over_dm"])

    # certified inequality on scanned compact set
    # for all scanned m with |m-1|>=eps: sup_norm >= c1*|m-1|
    # and also sup_norm >= c2*|m-1|^2
    check_linear = all(r["sup_norm"] + 1e-18 >= c1 * r["abs_dm"] for r in ratios)
    check_quad = all(r["sup_norm"] + 1e-18 >= c2 * r["abs_dm"] * r["abs_dm"] for r in ratios)

    payload = {
        "schema_version": "p2208_s1158_v1",
        "packet_id": "P2208",
        "stage_id": "S1158",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_SEPARATION_LOWER_BOUND_CERTIFICATE",
        "strict_nu_branch_separation_lower_bound_certificate": {
            "certificate_id": "STRICT_NU_BRANCH_SEPARATION_LOWER_BOUND_CERTIFICATE_V1",
            "source_packet": str(IN_2207.relative_to(ROOT)),
            "scan_exclusion_radius": eps,
            "min_linear_coefficient_c1": c1,
            "min_quadratic_coefficient_c2": c2,
            "worst_case_row_for_linear_bound": worst,
            "certified_bounds": {
                "linear_bound_holds_on_scan": check_linear,
                "quadratic_bound_holds_on_scan": check_quad,
                "statement": "For scanned m with |m-1|>=eps on compact d-domain: sup_norm >= c1*|m-1| and sup_norm >= c2*|m-1|^2",
            },
            "theorem_scope_limit": "numeric lower-bound certificate on scanned compact set only; not a global all-background theorem",
        },
        "recommended_next_honest_step": {
            "id": "P2209_candidate",
            "goal": "combine lower-bound certificate with transport residual integral to derive explicit non-vanishing threshold map vs |m-1|",
        },
        "gatekeeper_checks": {
            "lower_bound_certificate_exported": True,
            "linear_bound_holds_on_scan": check_linear,
            "quadratic_bound_holds_on_scan": check_quad,
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
            "# P2208 S1158: strict nu-branch separation lower-bound certificate",
            "",
            f"- c1 (linear lower-bound coeff): `{c1:.12e}`",
            f"- c2 (quadratic lower-bound coeff): `{c2:.12e}`",
            f"- linear bound holds on scan: `{check_linear}`",
            f"- quadratic bound holds on scan: `{check_quad}`",
            "",
            "Lower-bound certificate is scan/compact-set scoped; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
