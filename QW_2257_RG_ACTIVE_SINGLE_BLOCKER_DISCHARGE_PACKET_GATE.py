#!/usr/bin/env python3
"""QW-2257: RG active single-blocker discharge packet gate.

Builds a reduced obligation packet after QW-2255 showed that active RG path
contains exactly one core blocker.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
CORE_BLOCKER = "RGGlobalWellPosednessAllScales_DerivedOrPending"


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def main() -> None:
    q2255 = load("report_qw2255_rg_active_path_blocker_reduction_gate.json")
    q2251 = load("report_qw2251_rg_export_obligation_execution_status_gate.json")
    q2249 = load("report_qw2249_rg_export_axiom_free_obligation_packet_gate.json")

    active_blockers = q2255.get("active_blockers", [])
    is_singleton = active_blockers == [CORE_BLOCKER]

    obligations = {
        "RG_ACTIVE_CORE_O1": {
            "statement": (
                "construct theorem RG_CanonicalAction_to_WellPosedness_EXPORT : "
                "(FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales"
            ),
            "target_replacement_symbol": CORE_BLOCKER,
            "required": True,
            "acceptance": "symbol exists and is machine-checkable in *.lean",
        },
        "RG_ACTIVE_CORE_O2": {
            "statement": (
                "prove RG_CanonicalAction_to_WellPosedness_EXPORT without local/global axioms "
                "and without *_DerivedOrPending references, then rerun QW-2255 with active_blockers=[]"
            ),
            "target_replacement_symbol": CORE_BLOCKER,
            "required": True,
            "acceptance": "axiom-free proof refs + QW-2255 rerun active_blockers empty",
        },
    }

    spec = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "RG_ACTIVE_PATH_SINGLE_BLOCKER_DISCHARGE",
        "source_reports": {
            "q2255": "report_qw2255_rg_active_path_blocker_reduction_gate.json",
            "q2251": "report_qw2251_rg_export_obligation_execution_status_gate.json",
            "q2249": "report_qw2249_rg_export_axiom_free_obligation_packet_gate.json",
        },
        "core_blocker": CORE_BLOCKER,
        "active_blockers": active_blockers,
        "obligations": obligations,
        "scope_boundary": {
            "single_core_blocker_eliminated": False,
            "dax1_non_axiomatic_provider_completed": False,
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    spec_path = ROOT / "spec_qw2257_rg_active_single_blocker_discharge_packet.json"
    spec_path.write_text(json.dumps(spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "spec_file": spec_path.name,
        "spec_file_sha256": sha256_file(spec_path),
        "active_blockers": active_blockers,
        "derived_from": {
            "q2255_verdict": q2255.get("verdict"),
            "q2251_verdict": q2251.get("verdict"),
            "q2249_verdict": q2249.get("verdict"),
        },
        "scope_boundary": spec["scope_boundary"],
    }

    proof_path = ROOT / "proof_object_qw2257_rg_active_single_blocker_discharge_packet.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    flags = {
        "q2255_active_reduction_present": q2255.get("verdict")
        == "RG_ACTIVE_PATH_BLOCKER_REDUCTION_GATE_PASS_PARTIAL_SINGLE_CORE_BLOCKER",
        "q2249_obligation_packet_present": q2249.get("verdict")
        == "RG_EXPORT_AXIOM_FREE_OBLIGATION_PACKET_GATE_PASS_PACKET_READY_EXPORT_PENDING",
        "q2251_execution_status_present": q2251.get("verdict")
        == "RG_EXPORT_OBLIGATION_EXECUTION_STATUS_GATE_PASS_PARTIAL_EXPORT_PENDING",
        "active_single_core_blocker_confirmed": is_singleton,
        "core_blocker_matches_expected_symbol": CORE_BLOCKER in active_blockers,
        "reduced_obligation_packet_written": spec_path.exists(),
        "proof_object_generated": proof_path.exists(),
        "single_core_blocker_eliminated": False,
        "dax1_non_axiomatic_provider_completed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2255_active_reduction_present"]
        and flags["q2249_obligation_packet_present"]
        and flags["q2251_execution_status_present"]
        and flags["active_single_core_blocker_confirmed"]
        and flags["core_blocker_matches_expected_symbol"]
        and flags["reduced_obligation_packet_written"]
        and flags["proof_object_generated"]
    )

    verdict = (
        "RG_ACTIVE_SINGLE_BLOCKER_DISCHARGE_PACKET_GATE_PASS_PACKET_READY_CORE_BLOCKER_PENDING"
        if core_ok
        else "RG_ACTIVE_SINGLE_BLOCKER_DISCHARGE_PACKET_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": spec["source_reports"],
        "core_blocker": CORE_BLOCKER,
        "active_blockers": active_blockers,
        "spec_file": spec_path.name,
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "obligation_ids": list(obligations.keys()),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "DISCHARGE_RG_ACTIVE_CORE_O1_O2_AND_RERUN_QW2255_WITH_EMPTY_ACTIVE_BLOCKERS",
    }

    out_json = ROOT / "report_qw2257_rg_active_single_blocker_discharge_packet_gate.json"
    out_md = ROOT / "RAPORT_QW2257_RG_ACTIVE_SINGLE_BLOCKER_DISCHARGE_PACKET_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2257: RG ACTIVE SINGLE-BLOCKER DISCHARGE PACKET GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- core_blocker: `{CORE_BLOCKER}`",
                f"- active_blockers: `{active_blockers}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "verdict": verdict,
                "pass_count": pass_count,
                "total_flags": total_flags,
                "core_blocker": CORE_BLOCKER,
                "active_blockers": active_blockers,
            }
        )
    )


if __name__ == "__main__":
    main()
