#!/usr/bin/env python3
"""P1289: R8 new-data or formal-strengthening checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1288", type=Path, default=GEN / "p1288_qw2191_r7_provider_validation_and_selector_split_test_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1289_qw2191_r8_new_data_or_formal_strengthening_summary.json")
    args = parser.parse_args()

    p1288 = _read(args.p1288)
    if p1288.get("next_priority") != "R8_NEW_DATA_OR_FORMAL_STRENGTHENING":
        raise SystemExit("P1289 requires next_priority=R8_NEW_DATA_OR_FORMAL_STRENGTHENING from P1288.")

    split_result = p1288.get("r7", {}).get("selector_split_test", {}).get("result")
    if split_result != "INCONCLUSIVE":
        raise SystemExit("P1289 requires unresolved selector split result INCONCLUSIVE from P1288.")

    plan = {
        "path_A_new_data": {
            "observable_extensions": ["delta_phase_curvature_v2", "cross_sector_noise_reduction"],
            "expected_margin_gain": 0.008,
            "status": "OPEN",
        },
        "path_B_formal_strengthening": {
            "lemma_target": "margin_amplification_under_strict_transport_control",
            "required_assumptions": ["transport_stability_uniformity", "bounded_curvature_fluctuation"],
            "status": "OPEN",
        },
    }

    out = {
        "packet": "P1289",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1288": str(args.p1288)},
        "r8_program": {
            "branching_plan": plan,
            "recommended_first_move": "path_A_new_data",
            "status": "PARTIAL_DISCHARGE",
        },
        "closure_policy": {
            "strict_core_selector_closure_allowed": False,
            "global_qw2191_closure_allowed": False,
        },
        "next_priority": "R8A_DATA_ACQUISITION_PROTOCOL_V2",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1289] wrote {args.out}; first_move={out['r8_program']['recommended_first_move']}")


if __name__ == "__main__":
    main()
