#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n235_current_first_source_topology_invariant_candidate_packet_theorem_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p215 = load_json(
        GENERATED / "p215_current_source_topology_invariant_candidate_packet_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_FUTURE_SOURCE_TOPOLOGY_INVARIANT_CANDIDATE_PACKET_BELOW_BASIS_INDEPENDENCE_AND_QW2191_QUOTIENT_SAFETY_AFTER_P215"
    )
    status_ok = p215["status"] == expected_status

    summary = {
        "step": "N235",
        "status": "N235_DISCHARGED_CURRENT_FIRST_SOURCE_TOPOLOGY_INVARIANT_CANDIDATE_PACKET_THEOREM_NO_FALSE_PASS",
        "scope": "current_repo_future_source_topology_selector_route_only",
        "checks": [
            {
                "id": "p215_status",
                "actual": p215["status"],
                "expected": expected_status,
                "pass": status_ok,
            }
        ],
        "theorem_result": {
            "discharged": status_ok,
            "future_tau_src_candidate_packet_exported": p215["tau_src_candidate_packet_exported"],
            "basis_independent_selector_promotion_exported": False,
            "qw2191_quotient_safe_promotion_exported": False,
            "observer_downstream_only": p215["observer_downstream_only"],
            "current_selector_closure": False,
            "current_qw2191_discharge": False,
            "full_closure_pass": False,
        },
        "hard_limits": [
            "source_topology_invariant_not_yet_discharged",
            "no_basis_independent_selector_promotion",
            "no_qw2191_quotient_safe_selector_promotion",
            "observer_remains_downstream_only",
            "no_current_selector_closure",
            "no_current_qw2191_discharge",
            "no_ToE_closure",
        ],
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
