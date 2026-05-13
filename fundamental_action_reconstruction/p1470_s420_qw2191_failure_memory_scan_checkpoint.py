#!/usr/bin/env python3
"""P1470 S4.20: scan known failed QW-2191 lanes and export anti-pattern memory."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
SUMMARY = GEN / "p1470_s420_qw2191_failure_memory_summary.json"

SOURCES = [
    "V1_INFORMATIONAL_VISCOSITY_HYPOTHESIS_AUDIT.md",
    "V2_VISCOSITY_PROXY_TO_PAIR1_REDUCTION_AUDIT.md",
    "V3_MINIMAL_PAIR_LEVEL_VISCOSITY_OPERATOR_AUDIT.md",
    "V4_PSI0_COUPLED_VISCOSITY_EFFECT_AUDIT.md",
    "V5_PSI0_VISCOSITY_ANTI_OVERCLAIM_BOUNDARY_AUDIT.md",
    "V6_PSI0_VISCOSITY_DISTINCT_EFFECT_AUDIT.md",
    "V7_INFORMATIONAL_VISCOSITY_CURRENT_BEST_SUPPORTED_CONCLUSION.md",
    "T6_ROUTE_FAMILY_CLOSURE_CERTIFICATE_THEOREM_SPEC.md",
    "T7_ROUTE_FAMILY_CLOSURE_CERTIFICATE_DISCHARGE_ATTEMPT.md",
    "T9_ROUTE_ADMISSIBILITY_GRAMMAR_DISCHARGE_ATTEMPT.md",
]


def main() -> None:
    anti_patterns = [
        {
            "id": "AP1_viscosity_as_primary_selector_source",
            "why_failed": "Viscosity lanes stayed secondary/anchor-imported and did not generate selector source.",
            "avoid_repetition_rule": "Do not claim viscosity-only closure for QW-2191.",
        },
        {
            "id": "AP2_proxy_without_pair1_selector_block",
            "why_failed": "Proxy reductions did not export a concrete pair1 selector-sector witness.",
            "avoid_repetition_rule": "Require explicit pair1 selector block before closure talk.",
        },
        {
            "id": "AP3_grammar_or_certificate_without_discharge",
            "why_failed": "Route grammar/certificate specs existed but discharge was not completed.",
            "avoid_repetition_rule": "Do not equate theorem spec with executed discharge.",
        },
        {
            "id": "AP4_local_pass_overclaim",
            "why_failed": "Local improvements can fail holdout/edge checks if generalized too early.",
            "avoid_repetition_rule": "Keep claims local until stress-edge and guardrails remain stable.",
        },
    ]

    source_hits = []
    for name in SOURCES:
        path = ROOT / name
        txt = path.read_text(encoding="utf-8")
        hit = {
            "file": name,
            "mentions_qw2191": "QW-2191" in txt,
            "mentions_no_claim": ("no claim" in txt.lower()) or ("not discharge" in txt.lower()) or ("not discharged" in txt.lower()),
        }
        source_hits.append(hit)

    summary = {
        "packet": "P1470",
        "status": "PASS_QW2191_FAILURE_MEMORY_LOCAL_ONLY",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "anti_patterns": anti_patterns,
        "source_hits": source_hits,
        "next_rule": "Every new QW-2191 proposal must reference P1470 anti-pattern compliance.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1470] status={summary['status']} anti_patterns={len(anti_patterns)}")


if __name__ == "__main__":
    main()
