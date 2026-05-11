#!/usr/bin/env python3
"""P1262: build strict-only external replication pack manifest with hashes."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"
ROOT = Path(__file__).resolve().parent

FILES = [
    GEN / "p1259_strict_only_prediction_risk_ledger_summary.json",
    GEN / "p1260_strict_only_benchmark_execution_summary.json",
    GEN / "p1261_strict_only_uncertainty_decomposition_summary.json",
]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=GEN / "p1262_strict_only_external_replication_pack_summary.json")
    args = parser.parse_args()

    manifest = []
    for p in FILES:
        if not p.exists():
            raise SystemExit(f"Missing required artifact: {p}")
        manifest.append({
            "path": str(p.relative_to(ROOT)),
            "sha256": _sha256(p),
        })

    out = {
        "packet": "P1262",
        "as_of": "2026-05-11",
        "lane": "STRICT_ONLY",
        "replication_manifest": manifest,
        "execution_recipe": [
            "python3 fundamental_action_reconstruction/p1259_strict_only_prediction_risk_ledger_checkpoint.py",
            "python3 fundamental_action_reconstruction/p1260_strict_only_benchmark_execution_checkpoint.py",
            "python3 fundamental_action_reconstruction/p1261_strict_only_uncertainty_decomposition_checkpoint.py",
            "python3 fundamental_action_reconstruction/p1262_strict_only_external_replication_pack_checkpoint.py",
        ],
        "closure_note": "Replication readiness does not imply strict-core/global theory closure.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1262] wrote {args.out} with {len(manifest)} hashed artifacts")


if __name__ == "__main__":
    main()
