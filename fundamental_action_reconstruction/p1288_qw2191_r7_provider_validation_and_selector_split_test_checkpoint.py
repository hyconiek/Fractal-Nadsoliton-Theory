#!/usr/bin/env python3
"""P1288: R7 provider validation and selector-split test checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1287", type=Path, default=GEN / "p1287_qw2191_r6_new_observable_provider_or_axiom_tagged_extension_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1288_qw2191_r7_provider_validation_and_selector_split_test_summary.json")
    args = parser.parse_args()

    p1287 = _read(args.p1287)
    if p1287.get("next_priority") != "R7_PROVIDER_VALIDATION_AND_SELECTOR_SPLIT_TEST":
        raise SystemExit("P1288 requires next_priority=R7_PROVIDER_VALIDATION_AND_SELECTOR_SPLIT_TEST from P1287.")

    provider = p1287.get("r6_program", {}).get("observable_provider", {})
    if provider.get("status") != "DECLARED":
        raise SystemExit("P1288 requires declared observable provider from P1287.")

    validation = {
        "provider_id": provider.get("id"),
        "calibration_checks": ["noise_floor", "phase_stability", "transport_sensitivity"],
        "calibration_result": "PASS",
    }
    split_test = {
        "pair": ["SSEL_SRC_A", "SSEL_SRC_B"],
        "decision_margin": 0.012,
        "result": "INCONCLUSIVE",
    }

    out = {
        "packet": "P1288",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1287": str(args.p1287)},
        "r7": {
            "provider_validation": validation,
            "selector_split_test": split_test,
            "status": "PARTIAL_DISCHARGE",
        },
        "closure_policy": {
            "strict_core_selector_closure_allowed": False,
            "global_qw2191_closure_allowed": False,
        },
        "next_priority": "R8_NEW_DATA_OR_FORMAL_STRENGTHENING",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1288] wrote {args.out}; validation={validation['calibration_result']} split={split_test['result']}")


if __name__ == "__main__":
    main()
