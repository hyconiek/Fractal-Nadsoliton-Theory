#!/usr/bin/env python3
"""QW-2569: strict anti-false-pass identity-finality frontier gate."""

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
    q2566 = load("report_qw2566_dual_kernel_identity_totality_provider_theorem_spec_gate.json")
    q2567 = load("report_qw2567_dual_kernel_identity_totality_provider_counterexample_search_gate.json")
    q2568 = load("report_qw2568_non_axiomatic_dual_kernel_identity_totality_provider_derivation_attempt_gate.json")

    verdicts = [q2566.get("verdict", ""), q2567.get("verdict", ""), q2568.get("verdict", "")]
    forbidden_tokens = [
        "FULL_CLOSURE_READY",
        "THEOREM_LEVEL_PASS",
        "ALL_STRICT_OBLIGATIONS_FULLY_CLOSED_TRUE",
    ]

    flags = {
        "q2566_theorem_spec_partial_ready": q2566.get("verdict")
        == "DUAL_KERNEL_IDENTITY_TOTALITY_PROVIDER_THEOREM_SPEC_GATE_PASS_PARTIAL_TERMINAL_LAYER_READY",
        "q2567_counterexample_search_partial_ready": q2567.get("verdict")
        == "DUAL_KERNEL_IDENTITY_TOTALITY_PROVIDER_COUNTEREXAMPLE_SEARCH_GATE_PASS_PARTIAL_NO_STRICT_COUNTEREXAMPLE_IN_BOUNDED_DOMAIN",
        "q2568_execution_partial_blocked": q2568.get("verdict")
        == "NON_AXIOMATIC_DUAL_KERNEL_IDENTITY_TOTALITY_PROVIDER_DERIVATION_ATTEMPT_GATE_PASS_PARTIAL_BLOCKED_BY_KERNEL_IDENTITY_FINALITY_PROVIDER_THEOREMS",
        "all_flags_keep_not_fully_closed": (
            q2566.get("flags", {}).get("all_strict_obligations_fully_closed") is False
            and q2567.get("flags", {}).get("all_strict_obligations_fully_closed") is False
            and q2568.get("flags", {}).get("all_strict_obligations_fully_closed") is False
        ),
        "no_forbidden_overclaim_tokens_in_verdicts": not any(
            token in verdict for token in forbidden_tokens for verdict in verdicts
        ),
        "next_step_declared_on_q2568": isinstance(q2568.get("required_next_step"), str)
        and len(q2568.get("required_next_step")) > 0,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for value in flags.values() if value)
    total_flags = len(flags)

    core_ok = (
        flags["q2566_theorem_spec_partial_ready"]
        and flags["q2567_counterexample_search_partial_ready"]
        and flags["q2568_execution_partial_blocked"]
        and flags["all_flags_keep_not_fully_closed"]
        and flags["no_forbidden_overclaim_tokens_in_verdicts"]
        and flags["next_step_declared_on_q2568"]
    )

    verdict = (
        "STRICT_ANTI_FALSE_PASS_IDENTITY_FINALITY_FRONTIER_GATE_PASS_WITH_BLOCKERS_EXPLICIT"
        if core_ok
        else "STRICT_ANTI_FALSE_PASS_IDENTITY_FINALITY_FRONTIER_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2566": "report_qw2566_dual_kernel_identity_totality_provider_theorem_spec_gate.json",
            "q2567": "report_qw2567_dual_kernel_identity_totality_provider_counterexample_search_gate.json",
            "q2568": "report_qw2568_non_axiomatic_dual_kernel_identity_totality_provider_derivation_attempt_gate.json",
        },
        "verdicts": verdicts,
        "forbidden_tokens": forbidden_tokens,
        "required_next_step": q2568.get("required_next_step"),
    }
    proof_path = ROOT / "proof_object_qw2569_strict_anti_false_pass_identity_finality_frontier_gate.json"
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
        "required_next_step": "EXTRACT_DUAL_KERNEL_IDENTITY_FINALITY_PROVIDER_MINIMAL_BLOCKER_CUT",
    }

    out_json = ROOT / "report_qw2569_strict_anti_false_pass_identity_finality_frontier_gate.json"
    out_md = ROOT / "RAPORT_QW2569_STRICT_ANTI_FALSE_PASS_IDENTITY_FINALITY_FRONTIER_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2569: STRICT ANTI FALSE PASS IDENTITY FINALITY FRONTIER GATE",
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
