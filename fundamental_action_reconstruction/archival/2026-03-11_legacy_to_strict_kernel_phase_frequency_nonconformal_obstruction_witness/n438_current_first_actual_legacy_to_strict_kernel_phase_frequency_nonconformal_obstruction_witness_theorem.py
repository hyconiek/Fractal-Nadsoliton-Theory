#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"

IN_P404_SUMMARY = GENERATED / "p404_phase_frequency_nonconformal_obstruction_witness_probe_summary.json"

OUT_JSON = GENERATED / "n438_phase_frequency_nonconformal_obstruction_witness_theorem.json"
OUT_SUMMARY = GENERATED / "n438_phase_frequency_nonconformal_obstruction_witness_theorem_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    exists = IN_P404_SUMMARY.exists()
    p404 = load_json(IN_P404_SUMMARY) if exists else None
    passed = bool(p404.get("pass", False)) if isinstance(p404, dict) else False

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "N438_ARCHIVAL",
        "archival": True,
        "status": "ARCHIVAL_SCRATCH_DRAFT_DO_NOT_CITE_AS_ACTIVE_EXECUTION_LANE",
        "goal": "Package (archival) the draft phase/frequency nonconformal obstruction witness statement, conditional on the archival probe output.",
        "inputs": {
            "p404_archival_probe_summary": str(IN_P404_SUMMARY),
        },
        "verdict": {
            "conditional_on_p404_pass": True,
            "p404_pass": passed,
            "phase_frequency_obstruction_witness_materialized_in_archive": bool(passed),
            "note": "ARCHIVAL ONLY: not an active-lane theorem/discharge.",
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
        "p404_pass": bool(passed),
        "phase_frequency_obstruction_witness_materialized_in_archive": bool(
            report["verdict"]["phase_frequency_obstruction_witness_materialized_in_archive"]
        ),
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

