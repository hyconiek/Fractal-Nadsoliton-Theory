#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "f652_first_exported_s_sel_int_strict_core_source_object_selector_acceptance_independence_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f647 = load_json(
        "fundamental_action_reconstruction/generated/f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt_summary.json"
    )
    n125 = load_json(
        "fundamental_action_reconstruction/generated/n125_current_selector_requirement_theory_level_acceptance_theorem_summary.json"
    )

    props = f647["strict_core_export_properties"]
    theorem = n125["theorem_result"]

    summary = {
        "stage": "F652",
        "lane": "first_exported_s_sel_int_strict_core_source_object_selector_acceptance_independence_only",
        "status": "F652_EXECUTED_FIRST_EXPORTED_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_SELECTOR_ACCEPTANCE_INDEPENDENCE_PACKET_NO_FALSE_PASS",
        "exported_object": "S_sel_int_strict_core_source_object_v1",
        "uses_axiom_lane_artifact": bool(props["uses_axiom_lane_artifact"]),
        "strict_core_only": bool(props["strict_core_only"]),
        "selector_requirement_accepted_at_theory_level": bool(theorem["selector_requirement_accepted_at_theory_level"]),
        "accepted_scope": str(theorem["accepted_scope"]),
        "strict_core_changed": bool(theorem["strict_core_changed"]),
        "full_admissibility_pass": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

