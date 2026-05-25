#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2099 = GEN / "p2099_s1049_strict_u1_same_scheme_lock_witness.json"
IN_2100 = GEN / "p2100_s1050_strict_u2_phase_space_quadrature_witness.json"
IN_2101 = GEN / "p2101_s1051_strict_u3_residue_positivity_uncertainty_witness.json"
OUT = GEN / "p2102_s1052_strict_task2_entry_gate_summary.json"
MD = GEN / "p2102_s1052_strict_task2_entry_gate_summary.md"

SCHEMA_VERSION = "p2102_s1052_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2099 = load(IN_2099)
    p2100 = load(IN_2100)
    p2101 = load(IN_2101)

    u1 = p2099.get("gatekeeper_checks", {}).get("u1_computed") is True
    u2 = p2100.get("gatekeeper_checks", {}).get("u2_computed") is True
    u3 = p2101.get("gatekeeper_checks", {}).get("u3_computed") is True
    gate_ready = u1 and u2 and u3

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2102",
        "stage_id": "S1052",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_TASK2_ENTRY_GATE_SUMMARY_WITH_TRACE__U1_U2_U3_READY"
            if gate_ready
            else "OPEN_STRICT_TASK2_ENTRY_GATE_SUMMARY_BLOCKED"
        ),
        "depends_on": {
            "p2099_present": p2099.get("_missing") is None,
            "p2100_present": p2100.get("_missing") is None,
            "p2101_present": p2101.get("_missing") is None,
        },
        "task2_entry_gate": {
            "channel": "graviton -> gauge_gauge",
            "U1_shared_rg_scheme_lock": "COMPUTED" if u1 else "OPEN",
            "U2_exact_discontinuity_integration": "COMPUTED" if u2 else "OPEN",
            "U3_positive_residue_uncertainty_table": "COMPUTED" if u3 else "OPEN",
            "gate_ready": gate_ready,
            "scope_limit": "Task2 entry readiness summary only; not full dressed Cutkosky closure theorem.",
        },
        "recommended_next_honest_step": {
            "id": "P2103_candidate",
            "goal": "start first dressed discontinuity backend import and compare Disc_dressed vs CutSum on common basis in same scheme",
        },
        "c3_gate_update": {
            "C3_task2_entry_gate_summary": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "u1_ready": u1,
            "u2_ready": u2,
            "u3_ready": u3,
            "gate_ready": gate_ready,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2102 S1052: strict Task2 entry gate summary",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Gate ready: `{gate_ready}`",
            "",
            "This stage summarizes U1/U2/U3 readiness for Task2 entry.",
            "No full dressed Cutkosky closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
