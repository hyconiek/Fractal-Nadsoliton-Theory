#!/usr/bin/env python3
"""QW-2452: dual deeper-provider minimal blocker-cut gate.

Extracts minimal blocker cut after QW-2451 strict non-axiomatic authoring/discharge
attempt and keeps strict anti-overclaim boundary.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent

L12_SYMBOL = "RG_KernelOperatorClosureToWellPosedness_Theorem"
L5_SYMBOL = "QFT_KernelOperatorClosureToPositivity_Theorem"


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def main() -> None:
    q2451 = load("report_qw2451_strict_non_axiomatic_dual_export_provider_authoring_and_discharge_attempt_gate.json")

    l12_unknown = set(q2451.get("unknown_identifiers", {}).get("l12", []))
    l5_unknown = set(q2451.get("unknown_identifiers", {}).get("l5", []))

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
        "q2451_execution_attempt_performed": q2451.get("flags", {}).get("execution_attempt_performed") is True,
        "q2451_blocked_by_deeper_provider_theorems": q2451.get("verdict")
        == "STRICT_NON_AXIOMATIC_DUAL_EXPORT_PROVIDER_AUTHORING_AND_DISCHARGE_ATTEMPT_GATE_PASS_PARTIAL_BLOCKED_BY_DEEPER_NON_AXIOMATIC_PROVIDER_THEOREMS",
        "cut_size_two": len(cut) == 2,
        "l12_symbol_confirmed": cut[0]["present_in_unknown_identifiers"],
        "l5_symbol_confirmed": cut[1]["present_in_unknown_identifiers"],
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    if (
        flags["q2451_execution_attempt_performed"]
        and flags["q2451_blocked_by_deeper_provider_theorems"]
        and flags["l12_symbol_confirmed"]
        and flags["l5_symbol_confirmed"]
    ):
        verdict = "DUAL_DEEPER_PROVIDER_MINIMAL_BLOCKER_CUT_GATE_PASS_PARTIAL_TWO_SYMBOLS_ISOLATED"
        required_next_step = "ATTEMPT_NON_AXIOMATIC_DUAL_KERNEL_OPERATOR_CLOSURE_PROVIDER_DERIVATION"
    elif q2451.get("verdict") == "STRICT_NON_AXIOMATIC_DUAL_EXPORT_PROVIDER_AUTHORING_AND_DISCHARGE_ATTEMPT_GATE_PASS_PARTIAL_BLOCKED_BY_NO_RUNTIME":
        verdict = "DUAL_DEEPER_PROVIDER_MINIMAL_BLOCKER_CUT_GATE_PASS_PARTIAL_BLOCKED_BY_NO_RUNTIME"
        required_next_step = "ATTACH_RUNTIME_AND_RERUN_QW2451_BEFORE_QW2452"
    else:
        verdict = "DUAL_DEEPER_PROVIDER_MINIMAL_BLOCKER_CUT_GATE_FAIL"
        required_next_step = "REPAIR_QW2451_TO_QW2452_CHAIN"

    spec = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2451_strict_non_axiomatic_dual_export_provider_authoring_and_discharge_attempt_gate.json",
        "minimal_blocker_cut": cut,
        "scope_boundary": {
            "minimal_cut_extracted": True,
            "provider_discharge_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }
    spec_path = ROOT / "spec_qw2452_dual_deeper_provider_minimal_blocker_cut.json"
    spec_path.write_text(json.dumps(spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2451_strict_non_axiomatic_dual_export_provider_authoring_and_discharge_attempt_gate.json",
        "symbols": {
            "l12_expected": L12_SYMBOL,
            "l5_expected": L5_SYMBOL,
            "l12_unknown_identifiers": sorted(l12_unknown),
            "l5_unknown_identifiers": sorted(l5_unknown),
        },
        "n_cut_symbols": len(cut),
    }
    proof_path = ROOT / "proof_object_qw2452_dual_deeper_provider_minimal_blocker_cut.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2451": "report_qw2451_strict_non_axiomatic_dual_export_provider_authoring_and_discharge_attempt_gate.json",
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

    out_json = ROOT / "report_qw2452_dual_deeper_provider_minimal_blocker_cut_gate.json"
    out_md = ROOT / "RAPORT_QW2452_DUAL_DEEPER_PROVIDER_MINIMAL_BLOCKER_CUT_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2452: DUAL DEEPER PROVIDER MINIMAL BLOCKER CUT GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- n_cut_symbols: `{len(cut)}`",
                f"- cut_symbols: `{[row['symbol'] for row in cut]}`",
                "",
                "## Wniosek rygorystyczny",
                "- Minimalny dual deeper-provider blocker-cut jest jawny i ograniczony do 2 symboli kernel-operator frontier.",
                "- Brak podstaw do theorem-level/full-closure PASS.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "n_cut_symbols": len(cut)}))


if __name__ == "__main__":
    main()
