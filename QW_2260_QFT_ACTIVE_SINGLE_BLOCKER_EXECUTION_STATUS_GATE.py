#!/usr/bin/env python3
"""QW-2260: QFT active single-blocker execution-status gate.

Executes status check for reduced obligations QFT_ACTIVE_CORE_O1..O2.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def main() -> None:
    q2258 = load("report_qw2258_qft_active_single_blocker_discharge_packet_gate.json")
    q2256 = load("report_qw2256_qft_active_path_blocker_reduction_gate.json")
    q2244 = load("report_qw2244_qft_dax1_non_axiomatic_provider_attempt_gate.json")
    spec = load("spec_qw2258_qft_active_single_blocker_discharge_packet.json")

    export_symbol_exists = bool(q2244.get("flags", {}).get("canonical_export_symbol_exists", False))
    active_blockers_empty = q2256.get("active_blockers", []) == []

    obligation_status = {
        "QFT_ACTIVE_CORE_O1": {
            "required": True,
            "satisfied": export_symbol_exists,
            "evidence": {
                "canonical_export_symbol_exists": export_symbol_exists,
                "source_report": "report_qw2244_qft_dax1_non_axiomatic_provider_attempt_gate.json",
            },
        },
        "QFT_ACTIVE_CORE_O2": {
            "required": True,
            "satisfied": export_symbol_exists and active_blockers_empty,
            "evidence": {
                "canonical_export_symbol_exists": export_symbol_exists,
                "active_blockers_empty": active_blockers_empty,
                "source_reports": [
                    "report_qw2244_qft_dax1_non_axiomatic_provider_attempt_gate.json",
                    "report_qw2256_qft_active_path_blocker_reduction_gate.json",
                ],
            },
        },
    }

    n_total = len(obligation_status)
    n_satisfied = sum(1 for v in obligation_status.values() if v["satisfied"])

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2258": "report_qw2258_qft_active_single_blocker_discharge_packet_gate.json",
            "q2256": "report_qw2256_qft_active_path_blocker_reduction_gate.json",
            "q2244": "report_qw2244_qft_dax1_non_axiomatic_provider_attempt_gate.json",
            "spec": "spec_qw2258_qft_active_single_blocker_discharge_packet.json",
        },
        "obligation_status": obligation_status,
        "n_obligations_total": n_total,
        "n_obligations_satisfied": n_satisfied,
        "scope_boundary": {
            "single_core_blocker_eliminated": False,
            "dax1_non_axiomatic_provider_completed": False,
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2260_qft_active_single_blocker_execution_status.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    flags = {
        "q2258_packet_present": q2258.get("verdict")
        == "QFT_ACTIVE_SINGLE_BLOCKER_DISCHARGE_PACKET_GATE_PASS_PACKET_READY_CORE_BLOCKER_PENDING",
        "q2256_active_reduction_present": q2256.get("verdict")
        == "QFT_ACTIVE_PATH_BLOCKER_REDUCTION_GATE_PASS_PARTIAL_SINGLE_CORE_BLOCKER",
        "q2244_attempt_report_present": q2244.get("verdict")
        == "QFT_DAX1_NON_AXIOMATIC_PROVIDER_ATTEMPT_GATE_PASS_PARTIAL_CANONICAL_EXPORT_MISSING",
        "spec_packet_consistent": spec.get("core_blocker") == "PositivityToReconstruction_DerivedOrPending",
        "obligation_status_computed": True,
        "all_obligations_satisfied": n_satisfied == n_total,
        "single_core_blocker_eliminated": False,
        "dax1_non_axiomatic_provider_completed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2258_packet_present"]
        and flags["q2256_active_reduction_present"]
        and flags["q2244_attempt_report_present"]
        and flags["spec_packet_consistent"]
        and flags["obligation_status_computed"]
    )

    verdict = (
        "QFT_ACTIVE_SINGLE_BLOCKER_EXECUTION_STATUS_GATE_PASS_PARTIAL_CORE_BLOCKER_PENDING"
        if core_ok
        else "QFT_ACTIVE_SINGLE_BLOCKER_EXECUTION_STATUS_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "obligation_status": obligation_status,
        "n_obligations_total": n_total,
        "n_obligations_satisfied": n_satisfied,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "SATISFY_QFT_ACTIVE_CORE_O1_O2_BY_CONSTRUCTING_NON_AXIOMATIC_EXPORT_AND_EMPTYING_ACTIVE_BLOCKERS",
    }

    out_json = ROOT / "report_qw2260_qft_active_single_blocker_execution_status_gate.json"
    out_md = ROOT / "RAPORT_QW2260_QFT_ACTIVE_SINGLE_BLOCKER_EXECUTION_STATUS_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2260: QFT ACTIVE SINGLE-BLOCKER EXECUTION STATUS GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- obligations_satisfied: `{n_satisfied}/{n_total}`",
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
                "n_obligations_satisfied": n_satisfied,
                "n_obligations_total": n_total,
            }
        )
    )


if __name__ == "__main__":
    main()
