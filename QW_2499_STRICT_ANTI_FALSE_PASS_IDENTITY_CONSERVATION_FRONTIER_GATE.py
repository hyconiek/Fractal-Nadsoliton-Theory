#!/usr/bin/env python3
"""QW-2499: strict anti-false-pass identity-conservation frontier gate."""

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
    q2496 = load("report_qw2496_dual_kernel_identity_regularity_provider_theorem_spec_gate.json")
    q2497 = load("report_qw2497_dual_kernel_identity_regularity_provider_counterexample_search_gate.json")
    q2498 = load("report_qw2498_non_axiomatic_dual_kernel_identity_regularity_provider_derivation_attempt_gate.json")

    verdicts = [q2496.get("verdict", ""), q2497.get("verdict", ""), q2498.get("verdict", "")]
    forbidden_tokens = [
        "FULL_CLOSURE_READY",
        "THEOREM_LEVEL_PASS",
        "ALL_STRICT_OBLIGATIONS_FULLY_CLOSED_TRUE",
    ]

    flags = {
        "q2496_theorem_spec_partial_ready": q2496.get("verdict")
        == "DUAL_KERNEL_IDENTITY_REGULARITY_PROVIDER_THEOREM_SPEC_GATE_PASS_PARTIAL_TERMINAL_LAYER_READY",
        "q2497_counterexample_search_partial_ready": q2497.get("verdict")
        == "DUAL_KERNEL_IDENTITY_REGULARITY_PROVIDER_COUNTEREXAMPLE_SEARCH_GATE_PASS_PARTIAL_NO_STRICT_COUNTEREXAMPLE_IN_BOUNDED_DOMAIN",
        "q2498_execution_partial_blocked": q2498.get("verdict")
        == "NON_AXIOMATIC_DUAL_KERNEL_IDENTITY_REGULARITY_PROVIDER_DERIVATION_ATTEMPT_GATE_PASS_PARTIAL_BLOCKED_BY_KERNEL_IDENTITY_CONSERVATION_PROVIDER_THEOREMS",
        "all_flags_keep_not_fully_closed": (
            q2496.get("flags", {}).get("all_strict_obligations_fully_closed") is False
            and q2497.get("flags", {}).get("all_strict_obligations_fully_closed") is False
            and q2498.get("flags", {}).get("all_strict_obligations_fully_closed") is False
        ),
        "no_forbidden_overclaim_tokens_in_verdicts": not any(
            token in verdict for token in forbidden_tokens for verdict in verdicts
        ),
        "next_step_declared_on_q2498": isinstance(q2498.get("required_next_step"), str)
        and len(q2498.get("required_next_step")) > 0,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for value in flags.values() if value)
    total_flags = len(flags)

    core_ok = (
        flags["q2496_theorem_spec_partial_ready"]
        and flags["q2497_counterexample_search_partial_ready"]
        and flags["q2498_execution_partial_blocked"]
        and flags["all_flags_keep_not_fully_closed"]
        and flags["no_forbidden_overclaim_tokens_in_verdicts"]
        and flags["next_step_declared_on_q2498"]
    )

    verdict = (
        "STRICT_ANTI_FALSE_PASS_IDENTITY_CONSERVATION_FRONTIER_GATE_PASS_WITH_BLOCKERS_EXPLICIT"
        if core_ok
        else "STRICT_ANTI_FALSE_PASS_IDENTITY_CONSERVATION_FRONTIER_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2496": "report_qw2496_dual_kernel_identity_regularity_provider_theorem_spec_gate.json",
            "q2497": "report_qw2497_dual_kernel_identity_regularity_provider_counterexample_search_gate.json",
            "q2498": "report_qw2498_non_axiomatic_dual_kernel_identity_regularity_provider_derivation_attempt_gate.json",
        },
        "verdicts": verdicts,
        "forbidden_tokens": forbidden_tokens,
        "required_next_step": q2498.get("required_next_step"),
    }
    proof_path = ROOT / "proof_object_qw2499_strict_anti_false_pass_identity_conservation_frontier_gate.json"
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
        "required_next_step": "EXTRACT_DUAL_KERNEL_IDENTITY_CONSERVATION_PROVIDER_MINIMAL_BLOCKER_CUT",
    }

    out_json = ROOT / "report_qw2499_strict_anti_false_pass_identity_conservation_frontier_gate.json"
    out_md = ROOT / "RAPORT_QW2499_STRICT_ANTI_FALSE_PASS_IDENTITY_CONSERVATION_FRONTIER_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2499: STRICT ANTI FALSE PASS IDENTITY CONSERVATION FRONTIER GATE",
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
