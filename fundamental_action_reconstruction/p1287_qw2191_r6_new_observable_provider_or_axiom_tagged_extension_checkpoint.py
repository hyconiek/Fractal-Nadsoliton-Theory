#!/usr/bin/env python3
"""P1287: R6 new observable provider / axiom-tagged extension checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1286", type=Path, default=GEN / "p1286_qw2191_r5_selector_nonuniqueness_exclusion_attempt_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1287_qw2191_r6_new_observable_provider_or_axiom_tagged_extension_summary.json")
    args = parser.parse_args()

    p1286 = _read(args.p1286)
    if p1286.get("next_priority") != "R6_NEW_OBSERVABLE_PROVIDER_OR_AXIOM_TAGGED_EXTENSION":
        raise SystemExit("P1287 requires next_priority=R6_NEW_OBSERVABLE_PROVIDER_OR_AXIOM_TAGGED_EXTENSION from P1286.")

    obs_provider = {
        "id": "OBS_DELTA_PHASE_CURVATURE_V1",
        "class": "strict_lane_new_observable_provider",
        "target_pairs": [["SSEL_SRC_A", "SSEL_SRC_B"]],
        "discrimination_power": "CANDIDATE_POSITIVE",
        "status": "DECLARED",
    }

    axiom_extension = {
        "id": "AX_TAGGED_EXTENSION_R6_V1",
        "tag": "NON_STRICT_AUXILIARY",
        "enabled": False,
        "note": "Keep disabled while strict observable provider is under evaluation.",
    }

    out = {
        "packet": "P1287",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1286": str(args.p1286)},
        "r6_program": {
            "observable_provider": obs_provider,
            "axiom_tagged_extension": axiom_extension,
            "status": "PARTIAL_DISCHARGE",
        },
        "closure_policy": {
            "strict_core_selector_closure_allowed": False,
            "global_qw2191_closure_allowed": False,
        },
        "next_priority": "R7_PROVIDER_VALIDATION_AND_SELECTOR_SPLIT_TEST",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1287] wrote {args.out}; provider={obs_provider['id']} status={out['r6_program']['status']}")


if __name__ == "__main__":
    main()
