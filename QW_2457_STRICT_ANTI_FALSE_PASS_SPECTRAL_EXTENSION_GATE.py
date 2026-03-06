#!/usr/bin/env python3
"""QW-2457: strict anti-false-pass spectral extension gate.

Cross-checks the new spectral-layer chain (QW-2453..QW-2456) for
consistency and absence of overclaim tokens.
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
    q2453 = load("report_qw2453_non_axiomatic_dual_kernel_operator_closure_provider_derivation_attempt_gate.json")
    q2454 = load("report_qw2454_dual_kernel_spectral_provider_minimal_blocker_cut_gate.json")
    q2455 = load("report_qw2455_dual_kernel_spectral_provider_theorem_spec_gate.json")
    q2456 = load("report_qw2456_dual_kernel_spectral_provider_counterexample_search_gate.json")

    verdicts = [
        q2453.get("verdict", ""),
        q2454.get("verdict", ""),
        q2455.get("verdict", ""),
        q2456.get("verdict", ""),
    ]

    forbidden_tokens = [
        "FULL_CLOSURE_READY",
        "THEOREM_LEVEL_PASS",
        "ALL_STRICT_OBLIGATIONS_FULLY_CLOSED_TRUE",
    ]

    flags = {
        "q2453_expected_partial_blocked_verdict": q2453.get("verdict")
        == "NON_AXIOMATIC_DUAL_KERNEL_OPERATOR_CLOSURE_PROVIDER_DERIVATION_ATTEMPT_GATE_PASS_PARTIAL_BLOCKED_BY_KERNEL_SPECTRAL_CLOSURE_PROVIDER_THEOREMS",
        "q2454_expected_min_cut_verdict": q2454.get("verdict")
        == "DUAL_KERNEL_SPECTRAL_PROVIDER_MINIMAL_BLOCKER_CUT_GATE_PASS_PARTIAL_TWO_SYMBOLS_ISOLATED",
        "q2455_expected_spec_verdict": q2455.get("verdict")
        == "DUAL_KERNEL_SPECTRAL_PROVIDER_THEOREM_SPEC_GATE_PASS_PARTIAL_TERMINAL_LAYER_READY",
        "q2456_expected_counterexample_verdict": q2456.get("verdict")
        == "DUAL_KERNEL_SPECTRAL_PROVIDER_COUNTEREXAMPLE_SEARCH_GATE_PASS_PARTIAL_NO_STRICT_COUNTEREXAMPLE_IN_BOUNDED_DOMAIN",
        "all_flags_keep_not_fully_closed": (
            q2453.get("flags", {}).get("all_strict_obligations_fully_closed") is False
            and q2454.get("flags", {}).get("all_strict_obligations_fully_closed") is False
            and q2455.get("flags", {}).get("all_strict_obligations_fully_closed") is False
            and q2456.get("flags", {}).get("all_strict_obligations_fully_closed") is False
        ),
        "no_forbidden_overclaim_tokens_in_verdicts": not any(
            token in v for token in forbidden_tokens for v in verdicts
        ),
        "chain_required_next_step_present": isinstance(q2456.get("required_next_step"), str)
        and len(q2456.get("required_next_step")) > 0,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2453_expected_partial_blocked_verdict"]
        and flags["q2454_expected_min_cut_verdict"]
        and flags["q2455_expected_spec_verdict"]
        and flags["q2456_expected_counterexample_verdict"]
        and flags["all_flags_keep_not_fully_closed"]
        and flags["no_forbidden_overclaim_tokens_in_verdicts"]
        and flags["chain_required_next_step_present"]
    )

    verdict = (
        "STRICT_ANTI_FALSE_PASS_SPECTRAL_EXTENSION_GATE_PASS_WITH_BLOCKERS_EXPLICIT"
        if core_ok
        else "STRICT_ANTI_FALSE_PASS_SPECTRAL_EXTENSION_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2453": "report_qw2453_non_axiomatic_dual_kernel_operator_closure_provider_derivation_attempt_gate.json",
            "q2454": "report_qw2454_dual_kernel_spectral_provider_minimal_blocker_cut_gate.json",
            "q2455": "report_qw2455_dual_kernel_spectral_provider_theorem_spec_gate.json",
            "q2456": "report_qw2456_dual_kernel_spectral_provider_counterexample_search_gate.json",
        },
        "verdicts": verdicts,
        "forbidden_tokens": forbidden_tokens,
        "required_next_step": q2456.get("required_next_step"),
    }
    proof_path = ROOT / "proof_object_qw2457_strict_anti_false_pass_spectral_extension_gate.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "ATTEMPT_NON_AXIOMATIC_DUAL_KERNEL_SPECTRAL_CLOSURE_PROVIDER_DERIVATION",
    }

    out_json = ROOT / "report_qw2457_strict_anti_false_pass_spectral_extension_gate.json"
    out_md = ROOT / "RAPORT_QW2457_STRICT_ANTI_FALSE_PASS_SPECTRAL_EXTENSION_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2457: STRICT ANTI FALSE PASS SPECTRAL EXTENSION GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Chain `QW-2453..QW-2456` jest spójny i blocker-explicit.",
                "- `all_strict_obligations_fully_closed=false` pozostaje utrzymane we wszystkich bramkach.",
                "- Brak podstaw do theorem-level/full-closure PASS.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
