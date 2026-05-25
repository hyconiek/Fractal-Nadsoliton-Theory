#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2102 = GEN / "p2102_s1052_strict_task2_entry_gate_summary.json"
IN_1992 = GEN / "p1992_s942_strict_cutkosky_graviton_gauge_exact_phase_space_witness.json"
IN_1954 = GEN / "p1954_s904_strict_dressed_amplitude_nonavailability_theorem.json"
OUT = GEN / "p2103_s1053_strict_dressed_discontinuity_backend_import_contract.json"
MD = GEN / "p2103_s1053_strict_dressed_discontinuity_backend_import_contract.md"

SCHEMA_VERSION = "p2103_s1053_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2102 = load(IN_2102)
    p1992 = load(IN_1992)
    p1954 = load(IN_1954)

    gate_ready = (p2102.get("gatekeeper_checks", {}) or {}).get("gate_ready") is True

    blockers = [
        {
            "id": "D1",
            "name": "dressed_pole_residue_backend_on_common_basis",
            "status": "OPEN",
            "evidence": "p1954 nonavailability: missing dressed amplitude source rows",
        },
        {
            "id": "D2",
            "name": "disc_dressed_to_cutsum_same_scheme_comparator",
            "status": "OPEN",
            "evidence": "current chain has phase-space witness/proxy but no full dressed comparator export",
        },
        {
            "id": "D3",
            "name": "uncertainty_propagation_from_dressed_backend_to_optical_defect",
            "status": "OPEN",
            "evidence": "U3 table exists, but not yet fed by dressed backend covariance",
        },
    ]

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2103",
        "stage_id": "S1053",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_DRESSED_DISCONTINUITY_BACKEND_IMPORT_CONTRACT_WITH_TRACE"
            if gate_ready
            else "OPEN_STRICT_DRESSED_DISCONTINUITY_BACKEND_IMPORT_CONTRACT_BLOCKED"
        ),
        "depends_on": {
            "p2102_present": p2102.get("_missing") is None,
            "p1992_present": p1992.get("_missing") is None,
            "p1954_present": p1954.get("_missing") is None,
            "task2_entry_gate_ready": gate_ready,
        },
        "import_contract": {
            "target": "first dressed discontinuity backend import on common basis",
            "channel": "graviton -> gauge_gauge",
            "available_seed": {
                "p1992_phase_space_witness_present": p1992.get("_missing") is None,
                "p1954_formal_nonavailability_present": "formal_nonavailability_theorem" in p1954,
            },
            "blocker_register": blockers,
            "entry_rule": "Proceed to backend import execution only with explicit row-level source mapping for dressed poles/residues and same-scheme comparator contract.",
        },
        "recommended_next_honest_step": {
            "id": "P2104_candidate",
            "goal": "export D1 minimal dressed-pole/residue source object on common basis (or object-specific nonavailability) with machine-checkable row map",
        },
        "c3_gate_update": {
            "C3_dressed_discontinuity_backend_import_contract": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "task2_entry_gate_ready": gate_ready,
            "all_blockers_registered_open": all(b["status"] == "OPEN" for b in blockers),
            "dressed_backend_import_executed": False,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2103 S1053: strict dressed discontinuity backend-import contract",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Task2 gate ready: `{gate_ready}`",
            f"- Registered import blockers: `{len(blockers)}`",
            "",
            "This stage formalizes the first dressed-backend import contract and blocker register.",
            "No full dressed Cutkosky closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
