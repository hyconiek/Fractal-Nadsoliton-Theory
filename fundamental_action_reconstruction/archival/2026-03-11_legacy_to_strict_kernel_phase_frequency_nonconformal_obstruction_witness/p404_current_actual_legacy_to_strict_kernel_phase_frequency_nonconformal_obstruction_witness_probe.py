#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"

IN_F326_SUMMARY = GENERATED / "f326_phase_frequency_nonconformal_obstruction_witness_packet_summary.json"

OUT_JSON = GENERATED / "p404_phase_frequency_nonconformal_obstruction_witness_probe.json"
OUT_SUMMARY = GENERATED / "p404_phase_frequency_nonconformal_obstruction_witness_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    exists = IN_F326_SUMMARY.exists()
    f326_summary: dict[str, Any] | None = load_json(IN_F326_SUMMARY) if exists else None

    explicit_phase_frequency_bridge_present = None
    obstructed = None
    if f326_summary is not None:
        explicit_phase_frequency_bridge_present = bool(
            f326_summary.get("explicit_phase_frequency_bridge_present", False)
        )
        obstructed = bool(
            f326_summary.get("phase_frequency_layer_obstructed_on_current_exports", False)
        )

    passed = bool(exists and (explicit_phase_frequency_bridge_present is False) and (obstructed is True))

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P404_ARCHIVAL",
        "archival": True,
        "status": "ARCHIVAL_SCRATCH_DRAFT_DO_NOT_CITE_AS_ACTIVE_EXECUTION_LANE",
        "goal": "Probe (archival) that the draft F326 archival computation reports phase/frequency bridge absence on current exports, without active-lane promotion.",
        "inputs": {
            "f326_archival_summary": str(IN_F326_SUMMARY),
        },
        "checks": {
            "f326_summary_exists": exists,
            "explicit_phase_frequency_bridge_present": explicit_phase_frequency_bridge_present,
            "phase_frequency_layer_obstructed_on_current_exports": obstructed,
        },
        "verdict": {
            "pass": passed,
            "note": "ARCHIVAL ONLY: this is not an active-lane discharge.",
        },
        "hard_limits": [
            "ARCHIVAL ONLY: route frozen by S2 (legacy kernel retirement decree).",
            "Does not claim strict-core discharge.",
            "Does not claim selector closure or QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": report["stage"],
        "archival": True,
        "pass": bool(passed),
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()

