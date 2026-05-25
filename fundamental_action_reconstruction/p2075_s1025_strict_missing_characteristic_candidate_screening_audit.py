#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2075_s1025_strict_missing_characteristic_candidate_screening_audit.json"
MD = GEN / "p2075_s1025_strict_missing_characteristic_candidate_screening_audit.md"

SCHEMA_VERSION = "p2075_s1025_v1"
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
    p2074 = load("p2074_s1024_strict_missing_characteristic_candidate_catalog_audit.json")
    ready = p2074.get("result_kind") == "PASS_MISSING_CHARACTERISTIC_CANDIDATE_CATALOG_AUDIT_WITH_TRACE__C3_STILL_OPEN"

    rows = ((p2074.get("candidate_catalog") or {}).get("rows") or []) if isinstance(p2074, dict) else []

    screened = []
    for r in rows:
        cid = r.get("candidate_id", "unknown")
        # conservative ranking heuristic for next numeric stage (not theorem-grade)
        if cid == "mc_phase_lock":
            priority = 1
            rationale = "direct phase-coupling candidate; highest bridge-signal priority"
        elif cid == "mc_grad_cross":
            priority = 2
            rationale = "kinetic-mixing candidate; medium bridge-signal priority"
        else:
            priority = 3
            rationale = "higher structure complexity; deferred until top candidates tested"

        screened.append(
            {
                "candidate_id": cid,
                "priority_rank": priority,
                "screening_outcome": "RETAIN_FOR_NUMERIC_FIT_STAGE",
                "screening_rationale": rationale,
                "selector_qw2191_screening": r.get("selector_qw2191_screening", "QW2191_UNSPECIFIED"),
                "legacy_strict_alignment_signal": r.get("legacy_strict_alignment_signal", "OPEN"),
            }
        )

    screened.sort(key=lambda x: x["priority_rank"])

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2075",
        "stage_id": "S1025",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_MISSING_CHARACTERISTIC_CANDIDATE_SCREENING_AUDIT_WITH_TRACE__C3_STILL_OPEN"
            if ready
            else "OPEN_MISSING_CHARACTERISTIC_CANDIDATE_SCREENING_AUDIT_BLOCKED"
        ),
        "depends_on": {"p2074_present": p2074.get("_missing") is None, "preconditions_ready": ready},
        "input_hashes": {
            "p2074_json_sha256": file_sha256(GEN / "p2074_s1024_strict_missing_characteristic_candidate_catalog_audit.json"),
        },
        "screening_table": {
            "candidate_count": len(screened),
            "rows": screened,
            "next_stage_contract": "run numeric fit on priority 1->2 candidates before any closure claims",
        },
        "c3_gate_update": {
            "C3_missing_characteristic_candidate_screening_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "screening_rows_nonempty": len(screened) > 0,
            "priority_ranking_monotone": all(screened[i]["priority_rank"] <= screened[i + 1]["priority_rank"] for i in range(len(screened) - 1)),
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
            "# P2075 S1025: missing-characteristic candidate screening audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Candidate count: `{len(screened)}`",
            "- Next step: numeric fit for top-priority candidates",
            "",
            "This stage ranks candidate operators for conservative next numeric testing.",
            "C3 remains OPEN (not discharged).",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
