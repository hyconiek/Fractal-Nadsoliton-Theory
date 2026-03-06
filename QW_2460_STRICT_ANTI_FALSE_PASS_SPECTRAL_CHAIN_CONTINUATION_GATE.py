#!/usr/bin/env python3
"""QW-2460: strict anti-false-pass spectral chain continuation gate."""

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
    q2457 = load("report_qw2457_strict_anti_false_pass_spectral_extension_gate.json")
    q2458 = load("report_qw2458_non_axiomatic_dual_kernel_spectral_closure_provider_derivation_attempt_gate.json")
    q2459 = load("report_qw2459_dual_kernel_spectral_invariance_provider_minimal_blocker_cut_gate.json")

    verdicts = [
        q2457.get("verdict", ""),
        q2458.get("verdict", ""),
        q2459.get("verdict", ""),
    ]
    forbidden_tokens = [
        "FULL_CLOSURE_READY",
        "THEOREM_LEVEL_PASS",
        "ALL_STRICT_OBLIGATIONS_FULLY_CLOSED_TRUE",
    ]

    flags = {
        "q2457_integrity_baseline_ok": q2457.get("verdict")
        == "STRICT_ANTI_FALSE_PASS_SPECTRAL_EXTENSION_GATE_PASS_WITH_BLOCKERS_EXPLICIT",
        "q2458_expected_partial_blocked_verdict": q2458.get("verdict")
        == "NON_AXIOMATIC_DUAL_KERNEL_SPECTRAL_CLOSURE_PROVIDER_DERIVATION_ATTEMPT_GATE_PASS_PARTIAL_BLOCKED_BY_KERNEL_SPECTRAL_INVARIANCE_PROVIDER_THEOREMS",
        "q2459_expected_min_cut_verdict": q2459.get("verdict")
        == "DUAL_KERNEL_SPECTRAL_INVARIANCE_PROVIDER_MINIMAL_BLOCKER_CUT_GATE_PASS_PARTIAL_TWO_SYMBOLS_ISOLATED",
        "all_flags_keep_not_fully_closed": (
            q2458.get("flags", {}).get("all_strict_obligations_fully_closed") is False
            and q2459.get("flags", {}).get("all_strict_obligations_fully_closed") is False
        ),
        "no_forbidden_overclaim_tokens_in_verdicts": not any(
            token in v for token in forbidden_tokens for v in verdicts
        ),
        "next_step_declared_on_q2459": isinstance(q2459.get("required_next_step"), str)
        and len(q2459.get("required_next_step")) > 0,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2457_integrity_baseline_ok"]
        and flags["q2458_expected_partial_blocked_verdict"]
        and flags["q2459_expected_min_cut_verdict"]
        and flags["all_flags_keep_not_fully_closed"]
        and flags["no_forbidden_overclaim_tokens_in_verdicts"]
        and flags["next_step_declared_on_q2459"]
    )

    verdict = (
        "STRICT_ANTI_FALSE_PASS_SPECTRAL_CHAIN_CONTINUATION_GATE_PASS_WITH_BLOCKERS_EXPLICIT"
        if core_ok
        else "STRICT_ANTI_FALSE_PASS_SPECTRAL_CHAIN_CONTINUATION_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2457": "report_qw2457_strict_anti_false_pass_spectral_extension_gate.json",
            "q2458": "report_qw2458_non_axiomatic_dual_kernel_spectral_closure_provider_derivation_attempt_gate.json",
            "q2459": "report_qw2459_dual_kernel_spectral_invariance_provider_minimal_blocker_cut_gate.json",
        },
        "verdicts": verdicts,
        "forbidden_tokens": forbidden_tokens,
        "required_next_step": q2459.get("required_next_step"),
    }
    proof_path = ROOT / "proof_object_qw2460_strict_anti_false_pass_spectral_chain_continuation_gate.json"
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
        "required_next_step": "ATTEMPT_NON_AXIOMATIC_DUAL_KERNEL_SPECTRAL_INVARIANCE_PROVIDER_DERIVATION",
    }

    out_json = ROOT / "report_qw2460_strict_anti_false_pass_spectral_chain_continuation_gate.json"
    out_md = ROOT / "RAPORT_QW2460_STRICT_ANTI_FALSE_PASS_SPECTRAL_CHAIN_CONTINUATION_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2460: STRICT ANTI FALSE PASS SPECTRAL CHAIN CONTINUATION GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Kontynuacja chainu spectral pozostaje blocker-explicit i bez overclaimu.",
                "- Brak podstaw do theorem-level/full-closure PASS.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
