"""P1321: strictification feasibility audit for P_sel_v1.

Checks whether premise P_sel_v1 can be promoted to strict-core using currently
exported O3 artifacts.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1321_strictification_feasibility_audit_report_v1.json"


def load_json(name: str) -> dict:
    return json.loads((GEN / name).read_text(encoding="utf-8"))


def main() -> None:
    p1318 = load_json("p1318_o35_independent_replay_adversarial_audit_report_v1.json")
    p1319 = load_json("p1319_o3_final_residual_slot_neutrality_probe_report_v1.json")
    p1320 = load_json("p1320_premise_augmented_selector_law_probe_report_v1.json")

    checks = {
        "o3_5_replay_pass": p1318.get("replay_status") == "PASS",
        "o3_5_adversarial_pass": p1318.get("adversarial_status") == "PASS",
        "strict_neutrality_proven": p1319.get("neutrality_holds") is True,
        "premise_augmented_closure_available": p1320.get("uniqueness_under_premise") is True,
        "strict_core_already_closed": p1320.get("qw2191_strict_status") == "CLOSED",
    }

    strictification_ready = (
        checks["o3_5_replay_pass"]
        and checks["o3_5_adversarial_pass"]
        and checks["strict_neutrality_proven"]
    )

    missing_obligations = []
    if not checks["strict_neutrality_proven"]:
        missing_obligations.append("RESIDUAL_SLOT_NEUTRALITY_OR_ELIMINATION")
    if not checks["strict_core_already_closed"]:
        missing_obligations.append("INTERNAL_STRICT_SELECTOR_SOURCE_EXPORT")

    payload = {
        "packet_id": "P1321_STRICTIFICATION_FEASIBILITY_AUDIT_REPORT_V1",
        "date_utc": "2026-05-12",
        "input_artifacts": [
            "p1318_o35_independent_replay_adversarial_audit_report_v1.json",
            "p1319_o3_final_residual_slot_neutrality_probe_report_v1.json",
            "p1320_premise_augmented_selector_law_probe_report_v1.json",
        ],
        "checks": checks,
        "strictification_ready": strictification_ready,
        "missing_obligations": missing_obligations,
        "verdict": "STRICTIFICATION_NOT_READY" if not strictification_ready else "STRICTIFICATION_READY",
        "qw2191_status": "OPEN_STRICT_CLOSED_NON_STRICT",
    }

    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
