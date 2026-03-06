#!/usr/bin/env python3
"""QW-2464: strict anti-false-pass invariance frontier gate."""

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
    q2461 = load("report_qw2461_dual_kernel_spectral_invariance_provider_theorem_spec_gate.json")
    q2462 = load("report_qw2462_dual_kernel_spectral_invariance_provider_counterexample_search_gate.json")
    q2463 = load("report_qw2463_non_axiomatic_dual_kernel_spectral_invariance_provider_derivation_attempt_gate.json")

    verdicts = [q2461.get("verdict", ""), q2462.get("verdict", ""), q2463.get("verdict", "")]
    forbidden_tokens = [
        "FULL_CLOSURE_READY",
        "THEOREM_LEVEL_PASS",
        "ALL_STRICT_OBLIGATIONS_FULLY_CLOSED_TRUE",
    ]

    flags = {
        "q2461_theorem_spec_partial_ready": q2461.get("verdict")
        == "DUAL_KERNEL_SPECTRAL_INVARIANCE_PROVIDER_THEOREM_SPEC_GATE_PASS_PARTIAL_TERMINAL_LAYER_READY",
        "q2462_counterexample_search_partial_ready": q2462.get("verdict")
        == "DUAL_KERNEL_SPECTRAL_INVARIANCE_PROVIDER_COUNTEREXAMPLE_SEARCH_GATE_PASS_PARTIAL_NO_STRICT_COUNTEREXAMPLE_IN_BOUNDED_DOMAIN",
        "q2463_execution_partial_blocked": q2463.get("verdict")
        == "NON_AXIOMATIC_DUAL_KERNEL_SPECTRAL_INVARIANCE_PROVIDER_DERIVATION_ATTEMPT_GATE_PASS_PARTIAL_BLOCKED_BY_KERNEL_INVARIANCE_IDENTITY_PROVIDER_THEOREMS",
        "all_flags_keep_not_fully_closed": (
            q2461.get("flags", {}).get("all_strict_obligations_fully_closed") is False
            and q2462.get("flags", {}).get("all_strict_obligations_fully_closed") is False
            and q2463.get("flags", {}).get("all_strict_obligations_fully_closed") is False
        ),
        "no_forbidden_overclaim_tokens_in_verdicts": not any(
            token in v for token in forbidden_tokens for v in verdicts
        ),
        "next_step_declared_on_q2463": isinstance(q2463.get("required_next_step"), str)
        and len(q2463.get("required_next_step")) > 0,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2461_theorem_spec_partial_ready"]
        and flags["q2462_counterexample_search_partial_ready"]
        and flags["q2463_execution_partial_blocked"]
        and flags["all_flags_keep_not_fully_closed"]
        and flags["no_forbidden_overclaim_tokens_in_verdicts"]
        and flags["next_step_declared_on_q2463"]
    )

    verdict = (
        "STRICT_ANTI_FALSE_PASS_INVARIANCE_FRONTIER_GATE_PASS_WITH_BLOCKERS_EXPLICIT"
        if core_ok
        else "STRICT_ANTI_FALSE_PASS_INVARIANCE_FRONTIER_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2461": "report_qw2461_dual_kernel_spectral_invariance_provider_theorem_spec_gate.json",
            "q2462": "report_qw2462_dual_kernel_spectral_invariance_provider_counterexample_search_gate.json",
            "q2463": "report_qw2463_non_axiomatic_dual_kernel_spectral_invariance_provider_derivation_attempt_gate.json",
        },
        "verdicts": verdicts,
        "forbidden_tokens": forbidden_tokens,
        "required_next_step": q2463.get("required_next_step"),
    }
    proof_path = ROOT / "proof_object_qw2464_strict_anti_false_pass_invariance_frontier_gate.json"
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
        "required_next_step": "EXTRACT_DUAL_KERNEL_INVARIANCE_IDENTITY_PROVIDER_MINIMAL_BLOCKER_CUT",
    }

    out_json = ROOT / "report_qw2464_strict_anti_false_pass_invariance_frontier_gate.json"
    out_md = ROOT / "RAPORT_QW2464_STRICT_ANTI_FALSE_PASS_INVARIANCE_FRONTIER_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2464: STRICT ANTI FALSE PASS INVARIANCE FRONTIER GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Theorem-spec + counterexample-search + execution-attempt sa spojne i blocker-explicit.",
                "- `all_strict_obligations_fully_closed=false` utrzymane twardo.",
                "- Brak podstaw do theorem-level/full-closure PASS.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
