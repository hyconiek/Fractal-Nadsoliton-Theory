#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2126 = GEN / "p2126_s1076_strict_cmp2_backend_evidence_weighted_posterior_calibration_audit.json"
OUT = GEN / "p2132_s1082_strict_cmp2_real_support_acquisition_gate.json"
MD = GEN / "p2132_s1082_strict_cmp2_real_support_acquisition_gate.md"

SCHEMA_VERSION = "p2132_s1082_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"
MIN_NEW_ROWS_FOR_REPLAY = 8


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2126 = load(IN_2126)

    base_rows = ((p2126.get("backend_evidence_weighted_posterior_predictive_calibration_audit", {}) or {}).get("rows") or [])
    n_base = len(base_rows)

    acquisition_path = GEN / "cmp2_backend_rows_extension_v1.json"
    ext = load(acquisition_path) if acquisition_path.exists() else {"rows": []}
    new_rows = ext.get("rows") or []
    n_new = len(new_rows)

    # Minimal structural sanity for future non-synthetic replay.
    required_keys = {
        "cmp_bin_index",
        "backend_s",
        "posterior_weights_backend_evidence",
        "posterior_predictive_covered",
    }
    valid_rows = 0
    for r in new_rows:
        if isinstance(r, dict) and required_keys.issubset(r.keys()):
            valid_rows += 1

    ready = n_new >= MIN_NEW_ROWS_FOR_REPLAY and valid_rows == n_new and n_base > 0

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2132",
        "stage_id": "S1082",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CMP2_REAL_SUPPORT_ACQUISITION_GATE_WITH_TRACE" if ready else "OPEN_STRICT_CMP2_REAL_SUPPORT_ACQUISITION_GATE_BLOCKED",
        "real_support_acquisition_gate": {
            "base_rows_from_p2126": n_base,
            "extension_path": str(acquisition_path.relative_to(ROOT)),
            "new_rows_detected": n_new,
            "new_rows_valid": valid_rows,
            "min_new_rows_for_replay": MIN_NEW_ROWS_FOR_REPLAY,
            "required_row_keys": sorted(required_keys),
            "acquisition_ready_for_non_synthetic_replay": ready,
            "scope_limit": "acquisition/readiness gate only; no theorem-grade closure",
        },
        "recommended_next_honest_step": {
            "id": "P2133_candidate",
            "goal": "if acquisition gate passes, merge real rows and rerun P2127-P2131 with non-synthetic support",
        },
        "gatekeeper_checks": {
            "base_rows_present": n_base > 0,
            "new_rows_detected": n_new > 0,
            "new_rows_schema_valid": valid_rows == n_new,
            "acquisition_ready_for_replay": ready,
            "full_d3_covariance_transport_proven": False,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2132 S1082: strict CMP2 real-support acquisition gate",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- New rows detected: `{n_new}`",
            "",
            "This stage checks readiness for non-synthetic replay on truly extended backend rows.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
