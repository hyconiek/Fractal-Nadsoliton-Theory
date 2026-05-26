#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
EXT = GEN / "cmp2_backend_rows_extension_v1.json"
TEMPLATE = GEN / "cmp2_backend_rows_extension_v1.template.json"
OUT = GEN / "p2152_s1102_strict_cmp2_real_extension_delivery_readiness_gate.json"
MD = GEN / "p2152_s1102_strict_cmp2_real_extension_delivery_readiness_gate.md"


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def required_row_keys(template: dict[str, Any]) -> list[str]:
    rows = template.get("rows", [])
    if not rows or not isinstance(rows, list) or not isinstance(rows[0], dict):
        return []
    return sorted(rows[0].keys())


def row_missing_keys(row: dict[str, Any], req: list[str]) -> list[str]:
    return [k for k in req if k not in row]


def main() -> None:
    GEN.mkdir(exist_ok=True)
    template = load_json(TEMPLATE)
    extension = load_json(EXT)

    req_keys = required_row_keys(template)
    rows = extension.get("rows", []) if isinstance(extension, dict) else []
    if not isinstance(rows, list):
        rows = []

    missing_by_row: dict[str, list[str]] = {}
    valid_rows = 0
    for i, row in enumerate(rows):
        if not isinstance(row, dict):
            missing_by_row[str(i)] = req_keys[:]
            continue
        missing = row_missing_keys(row, req_keys)
        if missing:
            missing_by_row[str(i)] = missing
        else:
            valid_rows += 1

    schema_ok = extension.get("schema_version") == "cmp2_backend_rows_extension_v1"
    has_min_rows = valid_rows >= 8
    delivery_ready = schema_ok and has_min_rows and len(missing_by_row) == 0

    result_kind = (
        "PASS_STRICT_CMP2_REAL_EXTENSION_DELIVERY_READINESS_GATE"
        if delivery_ready
        else "OPEN_STRICT_CMP2_REAL_EXTENSION_DELIVERY_READINESS_GATE"
    )

    payload = {
        "schema_version": "p2152_s1102_v1",
        "packet_id": "P2152",
        "stage_id": "S1102",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "real_extension_delivery_readiness_gate": {
            "extension_file": str(EXT.relative_to(ROOT)),
            "template_file": str(TEMPLATE.relative_to(ROOT)),
            "required_row_keys": req_keys,
            "row_count": len(rows),
            "valid_rows": valid_rows,
            "min_valid_rows_required": 8,
            "schema_ok": schema_ok,
            "has_min_rows": has_min_rows,
            "missing_by_row": missing_by_row,
            "delivery_ready": delivery_ready,
        },
        "recommended_next_honest_step": {
            "id": "P2153_candidate",
            "goal": "if delivery_ready is true, rerun P2132/P2133/P2134/P2135 then refresh P2141/P2146/P2150",
        },
        "gatekeeper_checks": {
            "readiness_gate_exported": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": False,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2152 S1102: strict CMP2 real-extension delivery readiness gate",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Delivery ready: `{delivery_ready}`",
                f"- Valid rows: `{valid_rows}` / min `{8}`",
                "",
                "No theorem-grade closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
