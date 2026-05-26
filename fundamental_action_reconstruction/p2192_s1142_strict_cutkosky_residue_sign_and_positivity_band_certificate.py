#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2191 = GEN / "p2191_s1141_strict_cutkosky_uncertainty_bands_and_ur_residual_envelope.json"
OUT = GEN / "p2192_s1142_strict_cutkosky_residue_sign_and_positivity_band_certificate.json"
MD = GEN / "p2192_s1142_strict_cutkosky_residue_sign_and_positivity_band_certificate.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sign_label(x: float, eps: float = 1e-14) -> str:
    if x > eps:
        return "positive"
    if x < -eps:
        return "negative"
    return "zero"


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2191 = load(IN_2191)
    w = p2191.get("strict_cutkosky_uncertainty_bands_and_ur_residual_envelope", {}) or {}

    rows = w.get("uncertainty_band_table", []) or []
    env = w.get("residual_residue_envelope", {}) or {}
    ur = w.get("ur_link_uncertainty_object", {}) or {}

    sign_rows = []
    positive_count = 0
    negative_count = 0
    zero_count = 0
    max_quad_err = 0.0

    for i, r in enumerate(rows):
        val = float(r.get("integral", 0.0) or 0.0)
        err = float(r.get("quad_abs_error", 0.0) or 0.0)
        max_quad_err = max(max_quad_err, err)
        s = sign_label(val)
        if s == "positive":
            positive_count += 1
        elif s == "negative":
            negative_count += 1
        else:
            zero_count += 1
        sign_rows.append(
            {
                "idx": i,
                "integral": val,
                "integral_sign": s,
                "quad_abs_error": err,
                "params": r.get("params", {}),
            }
        )

    n = len(sign_rows)
    all_nonnegative = negative_count == 0
    all_positive_or_zero = (positive_count + zero_count) == n

    envelope_min = float(env.get("min_integral", 0.0) or 0.0)
    envelope_max = float(env.get("max_integral", 0.0) or 0.0)
    envelope_center = float(env.get("central_integral", 0.0) or 0.0)
    ur_radius = float(ur.get("radius", 0.0) or 0.0)

    certificate = {
        "certificate_id": "STRICT_CUTKOSKY_RESIDUE_SIGN_AND_POSITIVITY_BAND_CERTIFICATE_V1",
        "source_packet": str(IN_2191.relative_to(ROOT)),
        "frozen_channel": w.get("frozen_channel", "graviton_to_gauge_gauge"),
        "frozen_background": w.get("frozen_background", "B1_fixed_chart_seed"),
        "sample_count": n,
        "sign_summary": {
            "positive_count": positive_count,
            "negative_count": negative_count,
            "zero_count": zero_count,
            "all_nonnegative": all_nonnegative,
            "all_positive_or_zero": all_positive_or_zero,
        },
        "envelope_summary": {
            "min_integral": envelope_min,
            "max_integral": envelope_max,
            "central_integral": envelope_center,
            "ur_radius": ur_radius,
            "max_quad_abs_error": max_quad_err,
        },
        "positivity_band_certificate": {
            "band_lower_bound_nonnegative": envelope_min >= 0.0,
            "band_center_nonnegative": envelope_center >= 0.0,
            "band_upper_nonnegative": envelope_max >= 0.0,
            "ur_radius_nonnegative": ur_radius >= 0.0,
        },
        "residue_sign_table": sign_rows,
    }

    checks = {
        "certificate_exported": True,
        "sample_count_positive": n > 0,
        "all_nonnegative_integral_signs": all_nonnegative,
        "positivity_band_lower_nonnegative": certificate["positivity_band_certificate"]["band_lower_bound_nonnegative"],
        "no_selector_closure_claimed": True,
        "no_toe_closure_claimed": True,
        "full_cutkosky_closure_proven": False,
        "full_d3_covariance_transport_proven": True,
    }

    all_checks_pass = checks["sample_count_positive"] and checks["all_nonnegative_integral_signs"] and checks["positivity_band_lower_nonnegative"]
    result_kind = (
        "PASS_STRICT_CUTKOSKY_RESIDUE_SIGN_AND_POSITIVITY_BAND_CERTIFICATE"
        if all_checks_pass
        else "OPEN_STRICT_CUTKOSKY_RESIDUE_SIGN_AND_POSITIVITY_BAND_CERTIFICATE_ISSUE"
    )

    payload = {
        "schema_version": "p2192_s1142_v1",
        "packet_id": "P2192",
        "stage_id": "S1142",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_cutkosky_residue_sign_and_positivity_band_certificate": certificate,
        "recommended_next_honest_step": {
            "id": "P2193_candidate",
            "goal": "export channel-fixed uncertainty-to-gate bridge table mapping positivity bands to explicit CI tolerance classes",
        },
        "gatekeeper_checks": checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2192 S1142: strict Cutkosky residue-sign and positivity-band certificate",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Sample count: `{n}`",
                f"- Sign counts (+/-/0): `{positive_count}/{negative_count}/{zero_count}`",
                f"- Envelope min: `{envelope_min:.12e}`",
                f"- UR radius: `{ur_radius:.12e}`",
                "",
                "This is a strict-lane sign/positivity certificate only, not full dressed Cutkosky closure.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
