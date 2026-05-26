#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2192 = GEN / "p2192_s1142_strict_cutkosky_residue_sign_and_positivity_band_certificate.json"
OUT = GEN / "p2193_s1143_strict_cutkosky_uncertainty_to_ci_tolerance_bridge.json"
MD = GEN / "p2193_s1143_strict_cutkosky_uncertainty_to_ci_tolerance_bridge.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2192 = load(IN_2192)
    cert = p2192.get("strict_cutkosky_residue_sign_and_positivity_band_certificate", {}) or {}

    env = cert.get("envelope_summary", {}) or {}
    signs = cert.get("sign_summary", {}) or {}

    ur_radius = float(env.get("ur_radius", 0.0) or 0.0)
    lower_nonneg = bool((cert.get("positivity_band_certificate", {}) or {}).get("band_lower_bound_nonnegative", False))
    negative_count = int(signs.get("negative_count", 0) or 0)

    if lower_nonneg and negative_count == 0 and ur_radius <= 0.001:
        tolerance_class = "TOL_A_STRICT"
        ci_recommendation = "PASS"
    elif lower_nonneg and negative_count == 0 and ur_radius <= 0.005:
        tolerance_class = "TOL_B_MONITORED"
        ci_recommendation = "WARN"
    else:
        tolerance_class = "TOL_C_BLOCK"
        ci_recommendation = "FAIL"

    bridge = {
        "bridge_id": "STRICT_CUTKOSKY_UNCERTAINTY_TO_CI_TOLERANCE_BRIDGE_V1",
        "source_packet": str(IN_2192.relative_to(ROOT)),
        "input_metrics": {
            "ur_radius": ur_radius,
            "negative_count": negative_count,
            "band_lower_bound_nonnegative": lower_nonneg,
        },
        "tolerance_policy": {
            "TOL_A_STRICT": "ur_radius <= 1e-3 and no negative rows and nonnegative lower band",
            "TOL_B_MONITORED": "ur_radius <= 5e-3 and no negative rows and nonnegative lower band",
            "TOL_C_BLOCK": "otherwise",
        },
        "output": {
            "tolerance_class": tolerance_class,
            "ci_recommendation": ci_recommendation,
            "release_gate_hint": "OPEN" if ci_recommendation in {"PASS", "WARN"} else "BLOCK",
        },
    }

    payload = {
        "schema_version": "p2193_s1143_v1",
        "packet_id": "P2193",
        "stage_id": "S1143",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CUTKOSKY_UNCERTAINTY_TO_CI_TOLERANCE_BRIDGE",
        "strict_cutkosky_uncertainty_to_ci_tolerance_bridge": bridge,
        "recommended_next_honest_step": {
            "id": "P2194_candidate",
            "goal": "add replayed multi-background comparison table for the same tolerance policy",
        },
        "gatekeeper_checks": {
            "ci_tolerance_bridge_exported": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2193 S1143: strict Cutkosky uncertainty-to-CI tolerance bridge",
                "",
                f"- Tolerance class: `{tolerance_class}`",
                f"- CI recommendation: `{ci_recommendation}`",
                f"- UR radius: `{ur_radius:.12e}`",
                "",
                "This is a strict-lane tolerance bridge only, not full Cutkosky closure.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
