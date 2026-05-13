#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

FAR = Path(__file__).resolve().parent
GEN = FAR / "generated"


def main() -> None:
    matrix = [
        {
            "review_claim": "No full SM derivation yet",
            "far_status": "consistent_with_repo_state",
            "evidence": [
                "A6_GAUGE_RECONSTRUCTION_SPEC.md",
                "P1361_FAR_CONSTANT_CLAIMS_SCOREBOARD_PACKET_PL.md",
            ],
            "note": "FAR has partial strict candidates/scaffold, not full strict-verified closure.",
        },
        {
            "review_claim": "No full GR derivation yet",
            "far_status": "consistent_with_repo_state",
            "evidence": [
                "P1360_FAR_PHYSICAL_CONSTANT_DERIVATIONS_SENSIBILITY_AND_POST_R81_PRECISION_PACKET_PL.md",
                "P1368_GAIN_SUPPRESS_SWEEP_AND_HOLDOUT_READINESS_PACKET_PL.md",
            ],
            "note": "GR lane still open at strict-verified level.",
        },
        {
            "review_claim": "Framework has formal discipline",
            "far_status": "confirmed",
            "evidence": [
                "P1362_STRICT_CANDIDATE_RESIDUAL_BENCHMARK_AND_UNCERTAINTY_BUDGET_PACKET_PL.md",
                "P1363_STRICT_UPGRADE_GATE_AND_BLOCKER_REGISTRY_PACKET_PL.md",
                "P1365_PULLREQUEST_ARTIFACT_COMPLETENESS_PACKET_PL.md",
            ],
            "note": "No-false-pass, upgrade gate, artifact completeness checks are operationalized.",
        },
        {
            "review_claim": "Need falsifiable predictions",
            "far_status": "partially_addressed",
            "evidence": [
                "P1358_STRICT_KERNEL_ONLY_CANDIDATE_VALUE_GENERATION_AND_FIRST_RESIDUAL_VERDICT_PL.md",
                "P1367_HEBBIAN_MODE_GAIN_TEST_AND_GOVERNANCE_LOOP_PACKET_PL.md",
                "P1368_GAIN_SUPPRESS_SWEEP_AND_HOLDOUT_READINESS_PACKET_PL.md",
            ],
            "note": "FAR now runs falsifiable residual tests, but current class still fails thresholds.",
        },
    ]

    out = {
        "packet": "P1369",
        "as_of": "2026-05-12",
        "matrix": matrix,
        "bridge_verdict": "FAR_HAS_BRIDGE_ATTEMPTS_AND_TESTS_BUT_NOT_YET_STRICT_VERIFIED_SM_GR_CLOSURE",
        "next_priority": "P1370_BRIDGE_EVIDENCE_COMPILER_WITH_SM_GR_MINIMUM_SUCCESS_CRITERIA",
    }

    path = GEN / "p1369_external_review_vs_far_bridge_status_summary.json"
    path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1369] wrote {path}")


if __name__ == "__main__":
    main()
