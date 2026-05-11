#!/usr/bin/env python3
"""P1290: R8A data-acquisition protocol v2 checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1289", type=Path, default=GEN / "p1289_qw2191_r8_new_data_or_formal_strengthening_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1290_qw2191_r8a_data_acquisition_protocol_v2_summary.json")
    args = parser.parse_args()

    p1289 = _read(args.p1289)
    if p1289.get("next_priority") != "R8A_DATA_ACQUISITION_PROTOCOL_V2":
        raise SystemExit("P1290 requires next_priority=R8A_DATA_ACQUISITION_PROTOCOL_V2 from P1289.")

    rec = p1289.get("r8_program", {}).get("recommended_first_move")
    if rec != "path_A_new_data":
        raise SystemExit("P1290 requires recommended_first_move=path_A_new_data from P1289.")

    protocol = {
        "id": "R8A_DAQ_V2",
        "observable_targets": ["delta_phase_curvature_v2", "cross_sector_noise_reduction"],
        "sample_size_min": 48,
        "acceptance_criteria": {
            "noise_floor_max": 0.006,
            "margin_gain_min": 0.008,
            "cross_sector_consistency_min": 0.95,
        },
        "status": "DECLARED",
    }

    out = {
        "packet": "P1290",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1289": str(args.p1289)},
        "r8a_protocol": protocol,
        "closure_policy": {
            "strict_core_selector_closure_allowed": False,
            "global_qw2191_closure_allowed": False,
        },
        "next_priority": "R8A_EXECUTION_AND_MARGIN_REEVALUATION",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1290] wrote {args.out}; protocol={protocol['id']} n_min={protocol['sample_size_min']}")


if __name__ == "__main__":
    main()
