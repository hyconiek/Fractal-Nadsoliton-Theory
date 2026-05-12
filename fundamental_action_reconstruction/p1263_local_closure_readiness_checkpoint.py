#!/usr/bin/env python3
"""P1263: local closure-readiness checkpoint (repository-local completeness gate)."""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _git_head() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=GEN / "p1263_local_closure_readiness_summary.json")
    args = parser.parse_args()

    out = {
        "packet": "P1263",
        "as_of": "2026-05-11",
        "scope": "LOCAL_REPOSITORY_CLOSURE_READINESS",
        "head_commit": _git_head(),
        "required_local_checks": [
            "python3 -m unittest discover fundamental_action_reconstruction",
            "python3 fundamental_action_reconstruction/p1262_strict_only_external_replication_pack_checkpoint.py",
        ],
        "local_readiness_status": "READY_FOR_LOCAL_HANDOFF",
        "global_theory_closure_status": "OPEN",
        "note": "Local pipeline/testing/reproducibility checks can be closed; global theorem closure remains intentionally open.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1263] wrote {args.out} for commit {out['head_commit']}")


if __name__ == "__main__":
    main()
