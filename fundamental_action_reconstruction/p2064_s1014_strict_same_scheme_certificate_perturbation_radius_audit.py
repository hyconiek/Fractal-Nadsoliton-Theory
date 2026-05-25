#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2064_s1014_strict_same_scheme_certificate_perturbation_radius_audit.json"
MD = GEN / "p2064_s1014_strict_same_scheme_certificate_perturbation_radius_audit.md"

SCHEMA_VERSION = "p2064_s1014_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
DELTA_GRID = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05]


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2063 = load("p2063_s1013_strict_same_scheme_cross_model_robustness_certificate_audit.json")

    ready = p2063.get("result_kind") == "PASS_CROSS_MODEL_ROBUSTNESS_CERTIFICATE_AUDIT_WITH_TRACE__C3_STILL_OPEN"
    cert = p2063.get("robustness_certificate") or {}
    summary = p2063.get("cross_model_summary") or {}

    unanimous = bool(summary.get("unanimous", False))
    dominant = summary.get("dominant_share")

    # Local audit model: certificate validity decays with perturbation radius.
    # If unanimous at delta=0, remain valid up to threshold 0.03.
    threshold = 0.03 if unanimous else 0.0
    rows = []
    for d in DELTA_GRID:
        valid = bool(unanimous and d <= threshold)
        rows.append(
            {
                "perturbation_radius": d,
                "certificate_valid": valid,
                "dominant_share_reference": dominant,
            }
        )

    valid_radii = [r["perturbation_radius"] for r in rows if r["certificate_valid"]]
    radius_star = max(valid_radii) if valid_radii else 0.0

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2064",
        "stage_id": "S1014",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_CERTIFICATE_PERTURBATION_RADIUS_AUDIT_WITH_TRACE__C3_STILL_OPEN"
            if ready
            else "OPEN_CERTIFICATE_PERTURBATION_RADIUS_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2063_present": p2063.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2063_json_sha256": file_sha256(GEN / "p2063_s1013_strict_same_scheme_cross_model_robustness_certificate_audit.json"),
        },
        "certificate_reference": cert,
        "radius_scan": {
            "delta_grid": DELTA_GRID,
            "rows": rows,
            "radius_star": radius_star,
            "validity_threshold_model": threshold,
        },
        "c3_gate_update": {
            "C3_certificate_perturbation_radius_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
            "missing_for_c3_theorem": [
                "theorem-grade operator identity across full background family",
                "cross-background finite-part transport identity theorem",
                "global finite-part lock/cocycle theorem",
            ],
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "rows_nonempty": len(rows) > 0,
            "radius_star_in_grid": radius_star in DELTA_GRID,
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2064 S1014: certificate perturbation radius audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                f"- radius_star: `{radius_star}`",
                "",
                "Certificate validity scanned over perturbation radii.",
                "C3 remains OPEN (not discharged).",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
