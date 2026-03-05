#!/usr/bin/env python3
"""
QW-2145: Final external-proof pending gate.

Purpose:
- aggregate QW-2141..QW-2144 and expose minimal remaining blocker count,
- verify that all preparatory theorem-level machinery is closed locally,
- keep final boundary explicit: external proof object attachment pending.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2145_final_external_proof_pending_gate.json"
OUT_MD = ROOT / "RAPORT_QW2145_FINAL_EXTERNAL_PROOF_PENDING_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    r2141 = load("report_qw2141_continuum_weak_distribution_proxy_gate.json")
    r2142 = load("report_qw2142_l13_formal_proof_obligation_export_gate.json")
    r2143 = load("report_qw2143_external_machine_check_packet_gate.json")
    r2144 = load("report_qw2144_local_machine_check_proof_object_gate.json")
    r2146_path = ROOT / "report_qw2146_external_machine_check_execution_gate.json"
    r2146 = json.loads(r2146_path.read_text(encoding="utf-8")) if r2146_path.exists() else None

    flags = {
        "l14_weak_distribution_proxy_closed": str(r2141.get("verdict", "")).startswith(
            "CONTINUUM_WEAK_DISTRIBUTION_PROXY_GATE_PASS"
        ),
        "l13_formal_obligation_export_closed": str(r2142.get("verdict", "")).startswith(
            "L13_FORMAL_PROOF_OBLIGATION_EXPORT_GATE_PASS"
        ),
        "external_machine_check_packet_ready": str(r2143.get("verdict", "")).startswith(
            "EXTERNAL_MACHINE_CHECK_PACKET_GATE_PASS"
        ),
        "local_machine_check_proof_object_closed": str(r2144.get("verdict", "")) in {
            "LOCAL_MACHINE_CHECK_PROOF_OBJECT_GATE_PASS_PARTIAL",
            "LOCAL_MACHINE_CHECK_PROOF_OBJECT_GATE_CLOSED_PASS",
        },
        "external_machine_checked_proof_object_attached": bool(
            r2144["flags"]["full_external_machine_checked_proof_attached"]
        ),
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = len(flags)

    pending = [k for k, v in flags.items() if not bool(v)]

    if all(bool(v) for v in flags.values()):
        verdict = "FINAL_EXTERNAL_PROOF_PENDING_GATE_CLOSED_PASS"
    elif (
        flags["l14_weak_distribution_proxy_closed"]
        and flags["l13_formal_obligation_export_closed"]
        and flags["external_machine_check_packet_ready"]
        and flags["local_machine_check_proof_object_closed"]
        and (not flags["external_machine_checked_proof_object_attached"])
    ):
        verdict = "FINAL_EXTERNAL_PROOF_PENDING_GATE_PASS_PARTIAL_ONE_BLOCKER"
    else:
        verdict = "FINAL_EXTERNAL_PROOF_PENDING_GATE_FAIL_OR_MULTI_BLOCKER"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2141": "report_qw2141_continuum_weak_distribution_proxy_gate.json",
            "q2142": "report_qw2142_l13_formal_proof_obligation_export_gate.json",
            "q2143": "report_qw2143_external_machine_check_packet_gate.json",
            "q2144": "report_qw2144_local_machine_check_proof_object_gate.json",
            "q2146": "report_qw2146_external_machine_check_execution_gate.json" if r2146 is not None else None,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "pending_blockers": pending,
        "n_pending_blockers": len(pending),
        "verdict": verdict,
        "required_next_step": (
            "NO_BLOCKER_REMAINING_IN_QW2145_SCOPE"
            if verdict == "FINAL_EXTERNAL_PROOF_PENDING_GATE_CLOSED_PASS"
            else (
                "ATTACH_EXTERNAL_MACHINE_CHECKED_PROOF_OBJECT_AND_HASH_TO_CLOSE_PENDING_BLOCKER"
                if verdict.startswith("FINAL_EXTERNAL_PROOF_PENDING_GATE_PASS")
                else "REPAIR_PREPARATORY_CHAIN_AND_RERUN_QW2145"
            )
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md = [
        "# RAPORT QW-2145: FINAL EXTERNAL PROOF PENDING GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        f"- pending blockers: `{len(pending)}`",
        "",
        "## Pending list",
    ]
    for p in pending:
        md.append(f"- `{p}`")
    md.extend(
        [
            "",
            "## Boundary",
            "- In `*_PASS_PARTIAL_ONE_BLOCKER`: external proof object is the only final blocker.",
            "- In `*_CLOSED_PASS`: no blocker remains in QW-2145 scope.",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(md), encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pending": pending}, ensure_ascii=False))


if __name__ == "__main__":
    main()
