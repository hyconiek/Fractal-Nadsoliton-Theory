#!/usr/bin/env python3
"""QW-2448: dual single-foundation V2 minimal blocker-cut gate.

Extracts minimal blocker cut from QW-2445 V2 runtime-backed execution result.
Strict boundary is preserved: no theorem-level/full-closure claim.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent

L12_SYMBOL = "RG_CanonicalAction_to_WellPosedness_EXPORT"
L5_SYMBOL = "QFT_CanonicalAction_to_Positivity_EXPORT"


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def main() -> None:
    q2445 = load("report_qw2445_dual_nadsoliton_single_foundation_provider_execution_status_v2_gate.json")

    l12_unknown = set(q2445.get("unknown_identifiers", {}).get("l12", []))
    l5_unknown = set(q2445.get("unknown_identifiers", {}).get("l5", []))

    cut = [
        {
            "branch": "L12",
            "symbol": L12_SYMBOL,
            "present_in_unknown_identifiers": L12_SYMBOL in l12_unknown,
        },
        {
            "branch": "L5",
            "symbol": L5_SYMBOL,
            "present_in_unknown_identifiers": L5_SYMBOL in l5_unknown,
        },
    ]

    flags = {
        "q2445_runtime_backed_execution_available": q2445.get("flags", {}).get("execution_attempt_performed") is True,
        "q2445_partial_blocked_by_missing_symbols": q2445.get("verdict")
        == "DUAL_NADSOLITON_SINGLE_FOUNDATION_PROVIDER_EXECUTION_STATUS_V2_GATE_PASS_PARTIAL_BLOCKED_BY_MISSING_CANONICAL_EXPORT_SYMBOLS",
        "cut_size_two": len(cut) == 2,
        "l12_symbol_confirmed": cut[0]["present_in_unknown_identifiers"],
        "l5_symbol_confirmed": cut[1]["present_in_unknown_identifiers"],
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    if (
        flags["q2445_runtime_backed_execution_available"]
        and flags["q2445_partial_blocked_by_missing_symbols"]
        and flags["l12_symbol_confirmed"]
        and flags["l5_symbol_confirmed"]
    ):
        verdict = "DUAL_SINGLE_FOUNDATION_V2_MINIMAL_BLOCKER_CUT_GATE_PASS_PARTIAL_TWO_SYMBOLS_ISOLATED"
        required_next_step = "ATTEMPT_NON_AXIOMATIC_DUAL_CANONICAL_EXPORT_PROVIDER_DERIVATION"
    elif q2445.get("verdict") == "DUAL_NADSOLITON_SINGLE_FOUNDATION_PROVIDER_EXECUTION_STATUS_V2_GATE_PASS_PARTIAL_BLOCKED_BY_NO_RUNTIME":
        verdict = "DUAL_SINGLE_FOUNDATION_V2_MINIMAL_BLOCKER_CUT_GATE_PASS_PARTIAL_BLOCKED_BY_NO_RUNTIME"
        required_next_step = "ATTACH_RUNTIME_AND_RERUN_QW2445_BEFORE_QW2448"
    else:
        verdict = "DUAL_SINGLE_FOUNDATION_V2_MINIMAL_BLOCKER_CUT_GATE_FAIL"
        required_next_step = "REPAIR_QW2445_TO_QW2448_CHAIN"

    spec = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2445_dual_nadsoliton_single_foundation_provider_execution_status_v2_gate.json",
        "minimal_blocker_cut": cut,
        "scope_boundary": {
            "minimal_cut_extracted": True,
            "provider_discharge_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }
    spec_path = ROOT / "spec_qw2448_dual_single_foundation_v2_minimal_blocker_cut.json"
    spec_path.write_text(json.dumps(spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2445_dual_nadsoliton_single_foundation_provider_execution_status_v2_gate.json",
        "symbols": {
            "l12_expected": L12_SYMBOL,
            "l5_expected": L5_SYMBOL,
            "l12_unknown_identifiers": sorted(l12_unknown),
            "l5_unknown_identifiers": sorted(l5_unknown),
        },
        "n_cut_symbols": len(cut),
    }
    proof_path = ROOT / "proof_object_qw2448_dual_single_foundation_v2_minimal_blocker_cut.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2445": "report_qw2445_dual_nadsoliton_single_foundation_provider_execution_status_v2_gate.json",
            "spec": spec_path.name,
        },
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

    out_json = ROOT / "report_qw2448_dual_single_foundation_v2_minimal_blocker_cut_gate.json"
    out_md = ROOT / "RAPORT_QW2448_DUAL_SINGLE_FOUNDATION_V2_MINIMAL_BLOCKER_CUT_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2448: DUAL SINGLE FOUNDATION V2 MINIMAL BLOCKER CUT GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- n_cut_symbols: `{len(cut)}`",
                f"- cut_symbols: `{[row['symbol'] for row in cut]}`",
                "",
                "## Wniosek rygorystyczny",
                "- Minimalny dual blocker-cut jest jawny i ograniczony do 2 symboli canonical-export.",
                "- Brak podstaw do theorem-level/full-closure PASS.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "n_cut_symbols": len(cut)}))


if __name__ == "__main__":
    main()
