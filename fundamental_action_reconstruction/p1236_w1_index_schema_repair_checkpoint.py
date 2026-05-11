#!/usr/bin/env python3
"""P1236: repair index schema fields and rerun frozen-set consistency check."""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=Path, default=GEN / "p1229_w1_minimal_proof_packet_index.json")
    parser.add_argument("--p1235-out", type=Path, default=GEN / "p1235_w1_frozen_set_consistency_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1236_w1_index_schema_repair_summary.json")
    args = parser.parse_args()

    idx = json.loads(args.index.read_text(encoding="utf-8"))
    idx["strict_closure_claim_allowed"] = False
    idx["theory_closure_status"] = "OPEN"
    args.index.write_text(json.dumps(idx, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    root = Path(__file__).resolve().parent
    subprocess.run(["python3", str(root / "p1235_w1_frozen_set_consistency_check.py"), "--out", str(args.p1235_out)], check=True)

    p1235 = json.loads(args.p1235_out.read_text(encoding="utf-8"))
    out = {
        "packet": "P1236",
        "as_of": "2026-05-11",
        "index_schema_repaired": True,
        "repaired_fields": ["strict_closure_claim_allowed", "theory_closure_status"],
        "post_repair_all_artifacts_consistent": p1235.get("all_artifacts_consistent"),
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Schema repair checkpoint completed over P1229 index artifact.",
    }
    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1236] consistent={out['post_repair_all_artifacts_consistent']} wrote {args.out}")


if __name__ == "__main__":
    main()
