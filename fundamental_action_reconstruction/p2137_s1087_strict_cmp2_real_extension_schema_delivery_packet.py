#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
TARGET = GEN / "cmp2_backend_rows_extension_v1.json"
TEMPLATE = GEN / "cmp2_backend_rows_extension_v1.template.json"
OUT = GEN / "p2137_s1087_strict_cmp2_real_extension_schema_delivery_packet.json"
MD = GEN / "p2137_s1087_strict_cmp2_real_extension_schema_delivery_packet.md"

SCHEMA_VERSION = "p2137_s1087_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"
MIN_ROWS = 8
REQ = [
    "cmp_bin_index",
    "backend_s",
    "posterior_weights_backend_evidence",
    "posterior_predictive_covered",
]


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "rows": []}
    return json.loads(path.read_text(encoding="utf-8"))


def validate_rows(rows: list[Any]) -> tuple[int, list[str]]:
    valid = 0
    errs: list[str] = []
    for i, r in enumerate(rows):
        if not isinstance(r, dict):
            errs.append(f"row[{i}] not object")
            continue
        miss = [k for k in REQ if k not in r]
        if miss:
            errs.append(f"row[{i}] missing keys: {miss}")
            continue
        if not isinstance(r.get("posterior_weights_backend_evidence"), dict):
            errs.append(f"row[{i}] posterior_weights_backend_evidence not dict")
            continue
        valid += 1
    return valid, errs


def main() -> None:
    GEN.mkdir(exist_ok=True)

    template = {
        "schema_version": "cmp2_backend_rows_extension_v1",
        "rows": [
            {
                "cmp_bin_index": 4,
                "backend_s": 0.625,
                "posterior_weights_backend_evidence": {
                    "M1_nn": 0.28,
                    "M2_monotone": 0.26,
                    "M3_nn_tiebreak": 0.24,
                    "M4_monotone_penalized": 0.22,
                },
                "posterior_predictive_covered": True,
            }
        ],
        "note": "replace mock row with real backend-exported rows only",
    }
    TEMPLATE.write_text(json.dumps(template, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    data = load(TARGET)
    rows = data.get("rows") or []
    valid, errs = validate_rows(rows)
    ready = valid == len(rows) and valid >= MIN_ROWS

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2137",
        "stage_id": "S1087",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CMP2_REAL_EXTENSION_SCHEMA_DELIVERY_PACKET_WITH_TRACE" if ready else "OPEN_STRICT_CMP2_REAL_EXTENSION_SCHEMA_DELIVERY_PACKET_BLOCKED",
        "real_extension_schema_delivery": {
            "target_file": str(TARGET.relative_to(ROOT)),
            "template_file": str(TEMPLATE.relative_to(ROOT)),
            "required_keys": REQ,
            "min_rows_required": MIN_ROWS,
            "rows_detected": len(rows),
            "rows_valid": valid,
            "validation_errors": errs[:20],
            "ready_for_p2132_p2133": ready,
            "scope_limit": "schema-delivery/readiness only; not theorem-grade closure",
        },
        "recommended_next_honest_step": {
            "id": "P2138_candidate",
            "goal": "after real rows are delivered and P2132/P2133 unlock, rerun P2134/P2135 and compare non-synthetic CI inflation stability",
        },
        "gatekeeper_checks": {
            "template_exported": True,
            "rows_detected": len(rows) > 0,
            "rows_schema_valid": valid == len(rows),
            "rows_count_sufficient": valid >= MIN_ROWS,
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
            "# P2137 S1087: strict CMP2 real extension schema-delivery packet",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Rows detected in target file: `{len(rows)}`",
            "",
            "This stage delivers schema/template for real extension rows and validates readiness for P2132/P2133 unlock.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
