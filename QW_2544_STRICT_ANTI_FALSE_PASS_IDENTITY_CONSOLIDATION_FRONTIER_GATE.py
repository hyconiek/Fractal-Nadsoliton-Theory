#!/usr/bin/env python3
"""QW-2544: strict anti-false-pass identity-consolidation frontier gate."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def load(path: str) -> dict:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def main() -> None:
    q2541 = load("report_qw2541_dual_kernel_identity_resilience_provider_theorem_spec_gate.json")
    q2542 = load("report_qw2542_dual_kernel_identity_resilience_provider_counterexample_search_gate.json")
    q2543 = load("report_qw2543_non_axiomatic_dual_kernel_identity_resilience_provider_derivation_attempt_gate.json")

    verdicts = [q2541.get("verdict", ""), q2542.get("verdict", ""), q2543.get("verdict", "")]
    forbidden_tokens = [
        "FULL_CLOSURE_READY",
        "THEOREM_LEVEL_PASS",
        "ALL_STRICT_OBLIGATIONS_FULLY_CLOSED_TRUE",
    ]

    flags = {
        "q2541_theorem_spec_partial_ready": q2541.get("verdict")
        == "DUAL_KERNEL_IDENTITY_RESILIENCE_PROVIDER_THEOREM_SPEC_GATE_PASS_PARTIAL_TERMINAL_LAYER_READY",
        "q2542_counterexample_search_partial_ready": q2542.get("verdict")
        == "DUAL_KERNEL_IDENTITY_RESILIENCE_PROVIDER_COUNTEREXAMPLE_SEARCH_GATE_PASS_PARTIAL_NO_STRICT_COUNTEREXAMPLE_IN_BOUNDED_DOMAIN",
        "q2543_execution_partial_blocked": q2543.get("verdict")
        == "NON_AXIOMATIC_DUAL_KERNEL_IDENTITY_RESILIENCE_PROVIDER_DERIVATION_ATTEMPT_GATE_PASS_PARTIAL_BLOCKED_BY_KERNEL_IDENTITY_CONSOLIDATION_PROVIDER_THEOREMS",
        "all_flags_keep_not_fully_closed": (
            q2541.get("flags", {}).get("all_strict_obligations_fully_closed") is False
            and q2542.get("flags", {}).get("all_strict_obligations_fully_closed") is False
            and q2543.get("flags", {}).get("all_strict_obligations_fully_closed") is False
        ),
        "no_forbidden_overclaim_tokens_in_verdicts": not any(
            token in verdict for token in forbidden_tokens for verdict in verdicts
        ),
        "next_step_declared_on_q2543": isinstance(q2543.get("required_next_step"), str)
        and len(q2543.get("required_next_step")) > 0,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for value in flags.values() if value)
    total_flags = len(flags)

    core_ok = (
        flags["q2541_theorem_spec_partial_ready"]
        and flags["q2542_counterexample_search_partial_ready"]
        and flags["q2543_execution_partial_blocked"]
        and flags["all_flags_keep_not_fully_closed"]
        and flags["no_forbidden_overclaim_tokens_in_verdicts"]
        and flags["next_step_declared_on_q2543"]
    )

    verdict = (
        "STRICT_ANTI_FALSE_PASS_IDENTITY_CONSOLIDATION_FRONTIER_GATE_PASS_WITH_BLOCKERS_EXPLICIT"
        if core_ok
        else "STRICT_ANTI_FALSE_PASS_IDENTITY_CONSOLIDATION_FRONTIER_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2541": "report_qw2541_dual_kernel_identity_resilience_provider_theorem_spec_gate.json",
            "q2542": "report_qw2542_dual_kernel_identity_resilience_provider_counterexample_search_gate.json",
            "q2543": "report_qw2543_non_axiomatic_dual_kernel_identity_resilience_provider_derivation_attempt_gate.json",
        },
        "verdicts": verdicts,
        "forbidden_tokens": forbidden_tokens,
        "required_next_step": q2543.get("required_next_step"),
    }
    proof_path = ROOT / "proof_object_qw2544_strict_anti_false_pass_identity_consolidation_frontier_gate.json"
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
        "required_next_step": "EXTRACT_DUAL_KERNEL_IDENTITY_CONSOLIDATION_PROVIDER_MINIMAL_BLOCKER_CUT",
    }

    out_json = ROOT / "report_qw2544_strict_anti_false_pass_identity_consolidation_frontier_gate.json"
    out_md = ROOT / "RAPORT_QW2544_STRICT_ANTI_FALSE_PASS_IDENTITY_CONSOLIDATION_FRONTIER_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2544: STRICT ANTI FALSE PASS IDENTITY CONSOLIDATION FRONTIER GATE",
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
