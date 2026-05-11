#!/usr/bin/env python3
"""P1299: R15 axiom-tagged closure policy checkpoint.

This packet allows recording a user/program preference that legacy->strict bridging
is not required for *axiom-augmented* closure narratives, while preserving strict-core
closure guardrails.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1298", type=Path, default=GEN / "p1298_qw2191_r14_b1_nb1_interface_theorem_draft_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1299_qw2191_r15_axiom_tagged_closure_policy_summary.json")
    args = parser.parse_args()

    p1298 = _read(args.p1298)
    if p1298.get("next_priority") != "R15_B1_NB1_OBLIGATION_MATRIX_AND_PROOF_PLAN":
        raise SystemExit("P1299 requires next_priority=R15_B1_NB1_OBLIGATION_MATRIX_AND_PROOF_PLAN from P1298.")

    if p1298.get("r14_interface_theorem", {}).get("status") != "THEOREM_INTERFACE_DRAFTED":
        raise SystemExit("P1299 requires THEOREM_INTERFACE_DRAFTED status from P1298.")

    out = {
        "packet": "P1299",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_PLUS_AXIOM_TAGGED_POLICY",
        "input": {"p1298": str(args.p1298)},
        "policy": {
            "legacy_kernel_status": "HISTORICAL_REFERENCE_ONLY",
            "legacy_kernel_class": "NON_STRICT_HISTORICAL",
            "bridge_to_legacy_allowed": False,
            "preferred_resolution_path": "NB1_NONBRIDGE",
            "bridge_required_for_strict_core_closure": False,
            "bridge_required_for_axiom_augmented_closure": False,
            "axiom_augmented_closure_label": "NON_STRICT",
            "note": "Legacy is treated as non-strict historical; bridge-to-legacy claims are disallowed. Closure statements must remain explicitly non-strict.",
        },
        "next_priority": "R16_NONBRIDGE_SEPARATION_SCOPE_AND_LIMIT_THEOREM",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1299] wrote {args.out}; label={out['policy']['axiom_augmented_closure_label']}")


if __name__ == "__main__":
    main()
