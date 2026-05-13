#!/usr/bin/env python3
"""P1441: bifurcation checkpoint after strict v3 failure (v4 strict vs non-strict premise branch)."""

from __future__ import annotations

import json
from pathlib import Path


def main() -> None:
    base = Path(__file__).resolve().parent
    gen = base / "generated"
    gen.mkdir(parents=True, exist_ok=True)

    summary = {
        "packet": "P1441",
        "status": "BIFURCATION_DECLARED",
        "route": "F-Nadsoliton=>L_SM+L_GR",
        "input_state": "P1440_FAIL_V3_ASYMMETRY_OPERATOR_INSUFFICIENT",
        "branch_A_strict": {
            "name": "v4_strict_internal_selector_source_attempt",
            "allowed": True,
            "requires_new_strict_operator": True,
            "global_closure_claim_allowed": False,
        },
        "branch_B_non_strict": {
            "name": "explicit_non_strict_selector_premise_branch",
            "allowed": True,
            "must_be_explicitly_labeled_non_strict": True,
            "global_closure_claim_allowed": False,
        },
        "legacy_import_used": False,
        "qw2191_boundary_state": "OPEN_UNIQUENESS_OBSTRUCTION",
        "scope_of_pass": "local_contract",
        "strict_core_qw2191_closed": False,
        "next_honest_step": "execute branch_A_v4_strict first; keep branch_B as clearly separated fallback",
    }

    out = gen / "p1441_selector_source_bifurcation_summary.json"
    out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print('[P1441] status=BIFURCATION_DECLARED')


if __name__ == '__main__':
    main()
