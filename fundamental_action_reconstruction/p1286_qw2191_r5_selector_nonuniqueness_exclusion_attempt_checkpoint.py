#!/usr/bin/env python3
"""P1286: R5 selector nonuniqueness exclusion attempt checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1285", type=Path, default=GEN / "p1285_qw2191_r4_strict_selector_source_identification_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1286_qw2191_r5_selector_nonuniqueness_exclusion_attempt_summary.json")
    args = parser.parse_args()

    p1285 = _read(args.p1285)
    if p1285.get("next_priority") != "R5_SELECTOR_NONUNIQUENESS_EXCLUSION_ATTEMPT":
        raise SystemExit("P1286 requires next_priority=R5_SELECTOR_NONUNIQUENESS_EXCLUSION_ATTEMPT from P1285.")

    classes = p1285.get("selector_source_identification", {}).get("minimal_source_family", [])
    if len(classes) < 2:
        raise SystemExit("P1286 requires at least two selector-source classes in minimal_source_family.")

    exclusion_matrix = [
        {"pair": ["SSEL_SRC_A", "SSEL_SRC_B"], "distinguishing_observable": "delta_phase_transport_signature", "result": "INCONCLUSIVE"}
    ]

    out = {
        "packet": "P1286",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1285": str(args.p1285)},
        "r5_exclusion_attempt": {
            "classes_tested": classes,
            "exclusion_matrix": exclusion_matrix,
            "nonuniqueness_excluded": False,
            "status": "PARTIAL_DISCHARGE",
        },
        "closure_policy": {
            "strict_core_selector_closure_allowed": False,
            "global_qw2191_closure_allowed": False,
        },
        "next_priority": "R6_NEW_OBSERVABLE_PROVIDER_OR_AXIOM_TAGGED_EXTENSION",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1286] wrote {args.out}; excluded={out['r5_exclusion_attempt']['nonuniqueness_excluded']}")


if __name__ == "__main__":
    main()
