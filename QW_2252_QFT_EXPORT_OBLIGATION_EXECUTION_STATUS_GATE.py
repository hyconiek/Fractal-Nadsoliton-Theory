#!/usr/bin/env python3
"""QW-2252: QFT export obligation execution status gate."""

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
    q2244 = load("report_qw2244_qft_dax1_non_axiomatic_provider_attempt_gate.json")
    q2246 = load("report_qw2246_qft_dax1_axiom_free_candidate_scan_gate.json")
    q2248 = load("report_qw2248_qft_export_axiomatic_dependency_gate.json")
    q2250_obj = load("proof_object_qw2250_qft_export_axiom_free_obligation_packet.json")

    obligations = q2250_obj.get("obligations", {})

    o1 = bool(q2244.get("flags", {}).get("canonical_export_symbol_exists", False))
    o2 = not bool(q2248.get("flags", {}).get("dependency_chain_hits_derived_or_pending", False))
    o3 = (
        (not bool(q2248.get("flags", {}).get("dependency_chain_hits_axiom_layer", False)))
        and bool(q2248.get("flags", {}).get("non_axiomatic_candidate_exists", False))
    )
    o4 = bool(q2244.get("flags", {}).get("canonical_export_symbol_exists", False))

    obligation_status = {
        "QFT_EXPORT_O1": {"required": True, "satisfied": o1},
        "QFT_EXPORT_O2": {"required": True, "satisfied": o2},
        "QFT_EXPORT_O3": {"required": True, "satisfied": o3},
        "QFT_EXPORT_O4": {"required": True, "satisfied": o4},
    }

    n_total = len(obligation_status)
    n_satisfied = sum(1 for v in obligation_status.values() if v["satisfied"])

    flags = {
        "q2250_packet_present": len(obligations) == 4,
        "q2244_attempt_report_present": q2244.get("verdict")
        == "QFT_DAX1_NON_AXIOMATIC_PROVIDER_ATTEMPT_GATE_PASS_PARTIAL_CANONICAL_EXPORT_MISSING",
        "q2246_scan_report_present": q2246.get("verdict")
        == "QFT_DAX1_AXIOM_FREE_CANDIDATE_SCAN_GATE_PASS_PARTIAL_NO_AXIOM_FREE_CANDIDATE",
        "q2248_dependency_report_present": q2248.get("verdict")
        == "QFT_EXPORT_AXIOMATIC_DEPENDENCY_GATE_PASS_PARTIAL_AXIOM_FREE_EXPORT_ABSENT",
        "obligation_status_computed": True,
        "all_obligations_satisfied": n_satisfied == n_total,
        "dax1_non_axiomatic_provider_completed": False,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2250_packet_present"]
        and flags["q2244_attempt_report_present"]
        and flags["q2246_scan_report_present"]
        and flags["q2248_dependency_report_present"]
        and flags["obligation_status_computed"]
    )

    verdict = (
        "QFT_EXPORT_OBLIGATION_EXECUTION_STATUS_GATE_PASS_PARTIAL_EXPORT_PENDING"
        if core_ok
        else "QFT_EXPORT_OBLIGATION_EXECUTION_STATUS_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2244": "report_qw2244_qft_dax1_non_axiomatic_provider_attempt_gate.json",
            "q2246": "report_qw2246_qft_dax1_axiom_free_candidate_scan_gate.json",
            "q2248": "report_qw2248_qft_export_axiomatic_dependency_gate.json",
            "q2250": "proof_object_qw2250_qft_export_axiom_free_obligation_packet.json",
        },
        "obligation_status": obligation_status,
        "n_obligations_total": n_total,
        "n_obligations_satisfied": n_satisfied,
        "scope_boundary": {
            "dax1_non_axiomatic_provider_completed": False,
            "c1_theorem_discharge_completed": False,
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2252_qft_export_obligation_execution_status.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

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
        "required_next_step": "SATISFY_QFT_EXPORT_O1_O2_O3_O4_BY_CONSTRUCTING_NON_AXIOMATIC_QFT_EXPORT_THEOREM",
    }

    out_json = ROOT / "report_qw2252_qft_export_obligation_execution_status_gate.json"
    out_md = ROOT / "RAPORT_QW2252_QFT_EXPORT_OBLIGATION_EXECUTION_STATUS_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2252: QFT EXPORT OBLIGATION EXECUTION STATUS GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- obligations satisfied: `{n_satisfied}/{n_total}`",
                f"- pass_count: `{pass_count}/{total_flags}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "obligations_satisfied": f"{n_satisfied}/{n_total}", "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
