#!/usr/bin/env python3
"""P1796 S746 checkpoint generator."""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT_CHECKPOINT = ROOT / "generated" / "p1796_s746_strict_sv1_to_sv5_executable_ledger_template_and_state_sync_checkpoint.json"
OUT_TEMPLATE = ROOT / "generated" / "p1796_s746_sv1_sv5_executable_ledger_template.json"


def checkpoint() -> dict:
    return {
        "packet_id": "P1796_S746",
        "as_of": str(date(2026, 5, 15)),
        "mode": "strict_only",
        "allowed_verdicts": ["PASS_ZERO", "OPEN_OBSTRUCTION_WITH_TRACE"],
        "required_classification_axes": ["scope", "representation", "artifact_level"],
        "state_sync_scope": ["SV1", "SV2", "SV3", "SV4", "SV5"],
        "forbidden_state_sync_scope": ["SV6", "SV7", "SV8"],
        "gate_locks": {
            "G_BW": "CAN_CHANGE_ONLY_WITH_PASS_ZERO_WITNESS",
            "G_BRST": "LOCKED_IF_G_BW_NOT_PASS_ZERO",
            "G_CUT": "LOCKED_IF_G_BRST_NOT_PASS",
        },
    }


def template() -> dict:
    row = {
        "ledger_id": "SVX_RUN_YYYYMMDD_HHMM",
        "sv_key": "SVX",
        "freeze_id": "FREEZE_ID_REQUIRED",
        "classification": {
            "scope": "LOCAL",
            "representation": "NONPROXY",
            "artifact_level": "FULL_EXPORT",
        },
        "verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
        "residual_or_witness": {
            "residual_value": None,
            "witness_id": None,
            "check_command": "python3 run_componentwise_check.py --sv SVX",
            "check_output_digest": "REPLACE_WITH_HASH",
        },
        "obstruction_trace": "REQUIRED_IF_OPEN",
        "upstream_dependencies": [],
        "downstream_locks": ["G_BW", "G_BRST", "G_CUT"],
    }
    return {"sv_ledgers": [row], "notes": "Duplicate row for SV1..SV5 with common freeze_id."}


def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write("\n")


def main() -> None:
    write_json(OUT_CHECKPOINT, checkpoint())
    write_json(OUT_TEMPLATE, template())
    print(OUT_CHECKPOINT)
    print(OUT_TEMPLATE)


if __name__ == "__main__":
    main()
