#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2090_s1040_strict_delta_bg_yf_projection_theorem_object_decision.json"
MD = GEN / "p2090_s1040_strict_delta_bg_yf_projection_theorem_object_decision.md"

SCHEMA_VERSION = "p2090_s1040_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


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
    p2089 = load("p2089_s1039_strict_delta_bg_yf_lift_readiness_nonavailability_refresh.json")
    p1965 = load("p1965_s915_strict_delta_bg_yf_eom_normal_form_extraction_nonavailability.json")

    ready = p2089.get("result_kind") == "PASS_STRICT_DELTA_BG_YF_LIFT_READINESS_NONAVAILABILITY_REFRESH_WITH_NEXT_STEP__C3_STILL_OPEN"
    p1965_formal = p1965.get("formal_nonavailability_theorem", {})

    req = p2089.get("requirement_register", {})
    r2_open = req.get("R2", {}).get("status") == "OPEN_MISSING_EXPORT"
    r3_open = req.get("R3", {}).get("status") == "OPEN_MISSING_EXPORT"

    missing_projection_rows = [
        "explicit full-covariant termwise projection rows from EOM bundle to {delta_R, delta_RicUU, delta_gradchi2}",
        "linked metric/Higgs/Yukawa projection operator contract proving no Omega_unexported remainder",
    ]

    theorem_mode = "NONAVAILABILITY_CURRENT_EXPORT_STATE" if (ready and r2_open and r3_open) else "INCONCLUSIVE"

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2090",
        "stage_id": "S1040",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_DELTA_BG_YF_PROJECTION_THEOREM_OBJECT_DECISION__CURRENT_EXPORT_NONAVAILABILITY__C3_STILL_OPEN"
            if theorem_mode == "NONAVAILABILITY_CURRENT_EXPORT_STATE"
            else "OPEN_STRICT_DELTA_BG_YF_PROJECTION_THEOREM_OBJECT_DECISION_INCONCLUSIVE"
        ),
        "depends_on": {
            "p2089_present": p2089.get("_missing") is None,
            "p1965_present": p1965.get("_missing") is None,
            "preconditions_ready": ready,
            "r2_open": r2_open,
            "r3_open": r3_open,
        },
        "input_hashes": {
            "p2089_json_sha256": file_sha256(GEN / "p2089_s1039_strict_delta_bg_yf_lift_readiness_nonavailability_refresh.json"),
            "p1965_json_sha256": file_sha256(GEN / "p1965_s915_strict_delta_bg_yf_eom_normal_form_extraction_nonavailability.json"),
        },
        "projection_theorem_object_decision": {
            "mode": theorem_mode,
            "availability_object_exported": False,
            "nonavailability_object_exported": theorem_mode == "NONAVAILABILITY_CURRENT_EXPORT_STATE",
            "missing_projection_rows": missing_projection_rows,
            "inherits_prior_nonavailability_statement": p1965_formal.get("statement", "MISSING"),
            "licensed_statement": (
                "Current export state does not provide the required projection theorem object from full EOM to DELTA_BG_Yf tensor basis; "
                "this is a current-state nonavailability result, not a no-go theorem."
            ),
        },
        "recommended_next_honest_step": {
            "id": "P2091_candidate",
            "goal": "export explicit projection-row pack or certify exact nonavailability of each missing row/operator",
            "decision_target": "close R2/R3 one-by-one with machine-checkable row-level witnesses",
        },
        "c3_gate_update": {
            "C3_delta_bg_yf_projection_theorem_object_decision": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "decision_mode_is_nonavailability": theorem_mode == "NONAVAILABILITY_CURRENT_EXPORT_STATE",
            "missing_rows_registered": len(missing_projection_rows) >= 2,
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2090 S1040: strict DELTA_BG_Yf projection theorem-object decision",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Decision mode: `{payload['projection_theorem_object_decision']['mode']}`",
            f"- Missing projection rows registered: `{len(missing_projection_rows)}`",
            "",
            "This stage publishes the theorem-object decision in current export state as nonavailability.",
            "No theorem-grade closure claim is made; C3 remains OPEN.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
