#!/usr/bin/env python3
"""P1293: R9 formal selector-source theorem draft checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1292", type=Path, default=GEN / "p1292_qw2191_r8b_selector_split_rerun_with_updated_margin_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1293_qw2191_r9_formal_selector_source_theorem_draft_summary.json")
    args = parser.parse_args()

    p1292 = _read(args.p1292)
    if p1292.get("next_priority") != "R9_FORMAL_SELECTOR_SOURCE_THEOREM_DRAFT":
        raise SystemExit("P1293 requires next_priority=R9_FORMAL_SELECTOR_SOURCE_THEOREM_DRAFT from P1292.")

    result = p1292.get("r8b_split_rerun", {}).get("report", {}).get("result")
    if result != "DECISIVE_FOR_SSEL_SRC_A":
        raise SystemExit("P1293 requires decisive R8B result DECISIVE_FOR_SSEL_SRC_A from P1292.")

    theorem = {
        "id": "THM_R9_STRICT_SELECTOR_SOURCE_A",
        "statement": "Under R8A/R8B acceptance conditions, strict selector source reduces to SSEL_SRC_A over the declared domain.",
        "status": "DRAFT",
        "assumptions": [
            "transport_stability_uniformity",
            "bounded_curvature_fluctuation",
            "observable_provider_calibration_pass",
        ],
    }

    out = {
        "packet": "P1293",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1292": str(args.p1292)},
        "r9_theorem_draft": theorem,
        "closure_policy": {
            "strict_core_selector_closure_allowed": False,
            "global_qw2191_closure_allowed": False,
        },
        "next_priority": "R10_FORMAL_PROOF_CHAIN_AND_COUNTERMODEL_SWEEP",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1293] wrote {args.out}; theorem={theorem['id']} status={theorem['status']}")


if __name__ == "__main__":
    main()
