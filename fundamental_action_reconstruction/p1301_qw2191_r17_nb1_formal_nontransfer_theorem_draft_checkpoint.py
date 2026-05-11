#!/usr/bin/env python3
"""P1301: R17 NB1 formal non-transfer theorem draft checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1300", type=Path, default=GEN / "p1300_qw2191_r16_nonbridge_separation_scope_and_limit_theorem_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1301_qw2191_r17_nb1_formal_nontransfer_theorem_draft_summary.json")
    args = parser.parse_args()

    p1300 = _read(args.p1300)
    if p1300.get("next_priority") != "R17_NB1_FORMAL_NONTRANSFER_THEOREM_DRAFT":
        raise SystemExit("P1301 requires next_priority=R17_NB1_FORMAL_NONTRANSFER_THEOREM_DRAFT from P1300.")

    r16 = p1300.get("r16_nonbridge_scope_limit_theorem", {})
    if r16.get("status") != "SCOPE_LIMIT_DRAFTED":
        raise SystemExit("P1301 requires SCOPE_LIMIT_DRAFTED status from P1300.")

    theorem = {
        "id": "NB1_NONTRANSFER_R17",
        "hypotheses": [
            "legacy_kernel_class_is_nonstrict_historical",
            "strict_lane_claims_use_only_strict_operational_primitives",
            "no_bridge_axiom_between_legacy_and_strict_kernels",
        ],
        "claim": "No closure-relevant role transfer from legacy parameter layer to strict-core theorem layer is admissible.",
        "proof_obligations": [
            "show_parameter_role_domain_mismatch",
            "show_absence_of_admissible_transport_functor",
            "show_qw2191_obstruction_persists_without_new_selector_source",
        ],
        "falsification_conditions": [
            "construct_strict_admissible_transport_map_with_preserved_roles",
            "export_internal_selector_source_eliminating_qw2191_obstruction",
        ],
        "status": "DRAFT_WITH_OBLIGATIONS",
    }

    out = {
        "packet": "P1301",
        "as_of": "2026-05-11",
        "lane": "NB1_NONBRIDGE_TRACK",
        "input": {"p1300": str(args.p1300)},
        "r17_nb1_nontransfer_theorem": theorem,
        "next_priority": "R18_NB1_NONTRANSFER_OBLIGATION_MATRIX_AND_PROOF_SKETCH",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1301] wrote {args.out}; status={theorem['status']}")


if __name__ == "__main__":
    main()
