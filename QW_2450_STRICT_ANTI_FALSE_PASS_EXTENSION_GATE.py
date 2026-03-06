#!/usr/bin/env python3
"""QW-2450: strict anti-false-pass extension gate.

Extends anti-overclaim control to the chain:
QW-2447 -> QW-2448 -> QW-2449.
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
    q2447 = load("report_qw2447_strict_anti_false_pass_integrity_gate.json")
    q2448 = load("report_qw2448_dual_single_foundation_v2_minimal_blocker_cut_gate.json")
    q2449 = load("report_qw2449_non_axiomatic_dual_canonical_export_provider_derivation_attempt_gate.json")

    verdicts = {
        "q2447": q2447.get("verdict", ""),
        "q2448": q2448.get("verdict", ""),
        "q2449": q2449.get("verdict", ""),
    }

    overclaim_tokens = [
        "FULL_CLOSURE_READY",
        "THEORY_FULLY_CLOSED",
        "THEOREM_LEVEL_FULL_PASS",
        "PASS_TERMINAL_THEORY_CLOSED",
    ]

    flags = {
        "q2447_integrity_pass_with_blockers_explicit": verdicts["q2447"]
        == "STRICT_ANTI_FALSE_PASS_INTEGRITY_GATE_PASS_WITH_BLOCKERS_EXPLICIT",
        "q2448_dual_minimal_cut_isolated": verdicts["q2448"]
        == "DUAL_SINGLE_FOUNDATION_V2_MINIMAL_BLOCKER_CUT_GATE_PASS_PARTIAL_TWO_SYMBOLS_ISOLATED",
        "q2449_blocked_by_no_non_axiomatic_provider_definition": verdicts["q2449"]
        == "NON_AXIOMATIC_DUAL_CANONICAL_EXPORT_PROVIDER_DERIVATION_ATTEMPT_GATE_PASS_PARTIAL_BLOCKED_BY_NO_NON_AXIOMATIC_PROVIDER_DEFINITION",
        "q2449_zero_non_axiomatic_definition_counts": q2449.get("counts", {}).get("n_rg_non_axiomatic_definitions") == 0
        and q2449.get("counts", {}).get("n_qft_non_axiomatic_definitions") == 0,
        "q2448_all_strict_obligations_fully_closed_false": q2448.get("flags", {}).get("all_strict_obligations_fully_closed") is False,
        "q2449_all_strict_obligations_fully_closed_false": q2449.get("flags", {}).get("all_strict_obligations_fully_closed") is False,
        "no_overclaim_token_in_extended_chain_verdicts": all(
            all(tok not in verdict for tok in overclaim_tokens) for verdict in verdicts.values()
        ),
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    hard_requirements = [
        "q2447_integrity_pass_with_blockers_explicit",
        "q2448_dual_minimal_cut_isolated",
        "q2449_blocked_by_no_non_axiomatic_provider_definition",
        "q2449_zero_non_axiomatic_definition_counts",
        "q2448_all_strict_obligations_fully_closed_false",
        "q2449_all_strict_obligations_fully_closed_false",
        "no_overclaim_token_in_extended_chain_verdicts",
    ]

    if all(flags[k] for k in hard_requirements):
        verdict = "STRICT_ANTI_FALSE_PASS_EXTENSION_GATE_PASS_WITH_BLOCKERS_EXPLICIT"
        required_next_step = "AUTHOR_AND_DISCHARGE_STRICT_NON_AXIOMATIC_DUAL_EXPORT_PROVIDERS"
    else:
        verdict = "STRICT_ANTI_FALSE_PASS_EXTENSION_GATE_FAIL_INCONSISTENT_CHAIN_OR_OVERCLAIM_RISK"
        required_next_step = "REPAIR_QW2447_QW2449_CHAIN_AND_RERUN"

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2447": "report_qw2447_strict_anti_false_pass_integrity_gate.json",
            "q2448": "report_qw2448_dual_single_foundation_v2_minimal_blocker_cut_gate.json",
            "q2449": "report_qw2449_non_axiomatic_dual_canonical_export_provider_derivation_attempt_gate.json",
        },
        "verdicts": verdicts,
        "hard_requirements": hard_requirements,
        "flags": flags,
        "scope_boundary": {
            "integrity_extension_only": True,
            "theorem_discharge_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }
    proof_path = ROOT / "proof_object_qw2450_strict_anti_false_pass_extension_gate.json"
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
        "required_next_step": required_next_step,
    }

    out_json = ROOT / "report_qw2450_strict_anti_false_pass_extension_gate.json"
    out_md = ROOT / "RAPORT_QW2450_STRICT_ANTI_FALSE_PASS_EXTENSION_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2450: STRICT ANTI FALSE PASS EXTENSION GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- extended_chain_verdicts: `{verdicts}`",
                "",
                "## Wniosek rygorystyczny",
                "- Rozszerzony łańcuch zachowuje status blocker-explicit.",
                "- Brak podstaw do full-closure PASS; granica theorem-level pozostaje jawna.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
