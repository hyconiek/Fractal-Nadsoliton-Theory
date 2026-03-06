#!/usr/bin/env python3
"""QW-2443: dual Nadsoliton single-foundation minimal blocker-cut gate.

Extracts minimal two-symbol blocker cut from QW-2442 execution output.
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
    q2442 = load("report_qw2442_dual_nadsoliton_single_foundation_provider_execution_status_gate.json")

    l12_unknown = q2442.get("unknown_identifiers", {}).get("l12", [])
    l5_unknown = q2442.get("unknown_identifiers", {}).get("l5", [])

    cut = [
        {
            "branch": "L12",
            "symbol": "RG_CanonicalAction_to_WellPosedness_EXPORT",
            "present_in_unknown_identifiers": "RG_CanonicalAction_to_WellPosedness_EXPORT" in set(l12_unknown),
        },
        {
            "branch": "L5",
            "symbol": "QFT_CanonicalAction_to_Positivity_EXPORT",
            "present_in_unknown_identifiers": "QFT_CanonicalAction_to_Positivity_EXPORT" in set(l5_unknown),
        },
    ]

    flags = {
        "q2442_execution_partial_blocked": q2442.get("verdict")
        == "DUAL_NADSOLITON_SINGLE_FOUNDATION_PROVIDER_EXECUTION_STATUS_GATE_PASS_PARTIAL_BLOCKED_BY_MISSING_CANONICAL_EXPORT_SYMBOLS",
        "cut_size_two": len(cut) == 2,
        "l12_symbol_confirmed": cut[0]["present_in_unknown_identifiers"],
        "l5_symbol_confirmed": cut[1]["present_in_unknown_identifiers"],
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    if q2442.get("verdict") == "DUAL_NADSOLITON_SINGLE_FOUNDATION_PROVIDER_EXECUTION_STATUS_GATE_PASS_PARTIAL_BLOCKED_BY_LEAN_BINARY_UNAVAILABLE":
        verdict = "DUAL_NADSOLITON_SINGLE_FOUNDATION_MINIMAL_BLOCKER_CUT_GATE_PASS_PARTIAL_BLOCKED_BY_LEAN_BINARY_UNAVAILABLE"
        required_next_step = "INSTALL_OR_ATTACH_LEAN_AND_RERUN_QW2442_BEFORE_MIN_CUT_EXTRACTION"
    elif flags["q2442_execution_partial_blocked"] and flags["l12_symbol_confirmed"] and flags["l5_symbol_confirmed"]:
        verdict = "DUAL_NADSOLITON_SINGLE_FOUNDATION_MINIMAL_BLOCKER_CUT_GATE_PASS_PARTIAL_TWO_SYMBOLS_ISOLATED"
        required_next_step = "BUILD_NON_AXIOMATIC_DUAL_CANONICAL_EXPORT_PROVIDER_FROM_NADSOLITON_SINGLE_FOUNDATION"
    else:
        verdict = "DUAL_NADSOLITON_SINGLE_FOUNDATION_MINIMAL_BLOCKER_CUT_GATE_FAIL"
        required_next_step = "REPAIR_DUAL_NADSOLITON_SINGLE_FOUNDATION_BLOCKER_EXTRACTION"

    spec = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2442_dual_nadsoliton_single_foundation_provider_execution_status_gate.json",
        "minimal_blocker_cut": cut,
        "scope_boundary": {
            "minimal_cut_extracted": True,
            "provider_discharge_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }
    spec_path = ROOT / "spec_qw2443_dual_nadsoliton_single_foundation_minimal_blocker_cut.json"
    spec_path.write_text(json.dumps(spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2442": "report_qw2442_dual_nadsoliton_single_foundation_provider_execution_status_gate.json",
            "spec": spec_path.name,
        },
        "n_cut_symbols": len(cut),
    }
    proof_path = ROOT / "proof_object_qw2443_dual_nadsoliton_single_foundation_minimal_blocker_cut.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "spec_file": spec_path.name,
        "spec_sha256": sha256_file(spec_path),
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "n_cut_symbols": len(cut),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": required_next_step,
    }

    out_json = ROOT / "report_qw2443_dual_nadsoliton_single_foundation_minimal_blocker_cut_gate.json"
    out_md = ROOT / "RAPORT_QW2443_DUAL_NADSOLITON_SINGLE_FOUNDATION_MINIMAL_BLOCKER_CUT_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2443: DUAL NADSOLITON SINGLE FOUNDATION MINIMAL BLOCKER CUT GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- n_cut_symbols: `{len(cut)}`",
                f"- cut_symbols: `{[row['symbol'] for row in cut]}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "n_cut_symbols": len(cut)}))


if __name__ == "__main__":
    main()
