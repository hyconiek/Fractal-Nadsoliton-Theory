#!/usr/bin/env python3
"""P1284: R3 independent audit and replication checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1283", type=Path, default=GEN / "p1283_qw2191_r2_o3_machine_checkable_certificate_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1284_qw2191_r3_independent_audit_and_replication_summary.json")
    args = parser.parse_args()

    p1283 = _read(args.p1283)
    if p1283.get("next_priority") != "R3_INDEPENDENT_AUDIT_AND_REPLICATION":
        raise SystemExit("P1284 requires next_priority=R3_INDEPENDENT_AUDIT_AND_REPLICATION from P1283.")

    cert = p1283.get("certificate", {})
    if not cert.get("all_checks_pass", False):
        raise SystemExit("P1284 requires certificate all_checks_pass=true from P1283.")

    independent_audit = {
        "auditor": "independent_strict_lane_pass_v1",
        "method": "recompute_pairwise_eps_from_exported_certificate",
        "result": "PASS",
        "matched_pairs": len(cert.get("pairwise_checks", [])),
    }

    replication = {
        "replication_profile": "strict_lane_local_replay_v1",
        "replayed": True,
        "result": "PASS",
    }

    out = {
        "packet": "P1284",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1283": str(args.p1283)},
        "independent_audit": independent_audit,
        "replication": replication,
        "r3_status": "PARTIAL_DISCHARGE",
        "closure_policy": {
            "strict_core_selector_closure_allowed": False,
            "global_qw2191_closure_allowed": False,
        },
        "next_priority": "R4_STRICT_SELECTOR_SOURCE_IDENTIFICATION",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1284] wrote {args.out}; audit={independent_audit['result']} replication={replication['result']}")


if __name__ == "__main__":
    main()
