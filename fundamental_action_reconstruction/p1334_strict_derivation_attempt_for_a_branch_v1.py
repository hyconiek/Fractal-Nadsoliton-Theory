"""P1334: strict-derivation attempt for A_branch_v1."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1334_strict_derivation_attempt_for_a_branch_v1_report_v1.json"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    p1326 = load("p1326_v4_replay_adversarial_reintegration_report_v1.json")
    p1328 = load("p1328_formal_export_registry_report_v1.json")
    p1333 = load("p1333_branch_admissibility_axiom_test_report_v1.json")
    p1336 = load("p1336_internal_source_v2_tiebreak_family_search_report_v1.json")

    checks = {
        "v4_replay_adversarial_ready": p1326.get("o3_reintegration_ready") is True,
        "l1_exported": p1328.get("exports", {}).get("theorem_V4_selector_consistency_formally_exported", False),
        "axiom_non_strict_closure_exists": p1333.get("qw2191_status_under_axiom") == "CLOSED_NON_STRICT",
        "strict_internal_source_for_A_branch_v1_exported": p1336.get("strict_internal_source_exportable", False),
    }

    derivation_succeeded = all(
        checks[k]
        for k in [
            "v4_replay_adversarial_ready",
            "l1_exported",
            "strict_internal_source_for_A_branch_v1_exported",
        ]
    )

    payload = {
        "packet_id": "P1334_STRICT_DERIVATION_ATTEMPT_FOR_A_BRANCH_V1_REPORT_V1",
        "date_utc": "2026-05-12",
        "checks": checks,
        "strict_derivation_succeeded": derivation_succeeded,
        "status": "STRICT_DERIVATION_ESTABLISHED" if derivation_succeeded else "STRICT_DERIVATION_NOT_ESTABLISHED",
        "qw2191_strict_status": "NOT_CLOSED",
        "qw2191_non_strict_status": "CLOSED_NON_STRICT" if checks["axiom_non_strict_closure_exists"] else "OPEN",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
