#!/usr/bin/env python3
"""P1448 checkpoint: chain consistency certificate for P1446 + P1447."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
P1446 = ROOT / "generated" / "p1446_s1_selector_margin_monotonicity_summary.json"
P1447 = ROOT / "generated" / "p1447_s2_transport_robustness_summary.json"
OUT = ROOT / "generated" / "p1448_s3_chain_consistency_certificate_summary.json"


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    a = load(P1446)
    b = load(P1447)

    status_ok = (a.get("status") == "PASS_LOCAL_ONLY") and (b.get("status") == "PASS_LOCAL_TRANSPORT_ONLY")

    scope_ok = (a.get("scope") == "LOCAL_ONLY_NON_GLOBAL_CLAIM") and (
        b.get("scope_of_pass") == "LOCAL_ONLY_NON_GLOBAL_CLAIM"
    )

    legacy_ok = (a.get("legacy_bridge_used") is False) and (b.get("legacy_bridge_used") is False)

    qw_flag = b.get("strict_core_qw2191_closed")
    qw_ok = qw_flag is False

    if not status_ok:
        verdict = "FAIL_CHAIN_STATUS"
    elif not scope_ok:
        verdict = "FAIL_CHAIN_SCOPE"
    elif not legacy_ok:
        verdict = "FAIL_CHAIN_LEGACY_FLAG"
    elif not qw_ok:
        verdict = "FAIL_CHAIN_QW2191_FLAG"
    else:
        verdict = "PASS_CHAIN_LOCAL_ONLY"

    out = {
        "packet": "P1448",
        "status": verdict,
        "inputs": {
            "p1446_status": a.get("status"),
            "p1447_status": b.get("status"),
        },
        "checks": {
            "status_alignment": status_ok,
            "scope_alignment_local_only": scope_ok,
            "legacy_bridge_absent": legacy_ok,
            "qw2191_not_closed": qw_ok,
        },
        "scope_of_pass": "LOCAL_ONLY_NON_GLOBAL_CLAIM" if verdict == "PASS_CHAIN_LOCAL_ONLY" else "NO_PASS",
        "legacy_bridge_used": False,
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1448] status={verdict} -> {OUT}")


if __name__ == "__main__":
    main()
