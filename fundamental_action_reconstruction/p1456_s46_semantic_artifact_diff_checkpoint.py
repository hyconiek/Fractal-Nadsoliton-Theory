#!/usr/bin/env python3
"""P1456 S4.6: local-only semantic guardrail diff for selected strict artifacts."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
SUMMARY = GEN / "p1456_s46_semantic_diff_summary.json"
OBSTRUCTION = GEN / "p1456_s46_semantic_diff_obstruction.json"

TARGETS = {
    "p1452": GEN / "p1452_s43_holdout_summary.json",
    "p1453": GEN / "p1453_s44_h2_rerun_summary.json",
    "p1455": GEN / "p1455_s45_generated_artifact_diff_summary.json",
}


def readj(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def extract_scope(x: dict) -> str:
    return str(x.get("scope_of_pass", x.get("scope", "")))


def main() -> None:
    violations = []
    rows = []

    for tag, path in TARGETS.items():
        data = readj(path)
        scope = extract_scope(data)
        strict_closed = bool(data.get("strict_core_qw2191_closed", False))
        legacy_bridge = bool(data.get("legacy_bridge_used", False))
        status = str(data.get("status", "UNKNOWN"))

        row = {
            "tag": tag,
            "file": path.name,
            "status": status,
            "scope": scope,
            "strict_core_qw2191_closed": strict_closed,
            "legacy_bridge_used": legacy_bridge,
        }
        rows.append(row)

        if "LOCAL_ONLY_NON_GLOBAL_CLAIM" not in scope:
            violations.append({"type": "scope_violation", **row})
        if strict_closed:
            violations.append({"type": "qw2191_premature_closure", **row})
        if legacy_bridge:
            violations.append({"type": "legacy_bridge_violation", **row})

    status = "PASS_SEMANTIC_GUARDRAIL_LOCAL_ONLY" if not violations else "FAIL_SEMANTIC_GUARDRAIL"
    summary = {
        "packet": "P1456",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "checked_artifacts": rows,
        "violation_count": len(violations),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if violations:
        obstruction = {
            "packet": "P1456",
            "status": status,
            "violations": violations,
            "rule": "immediate obstruction export on first semantic guardrail fail",
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"[P1456] status={status} violations={len(violations)} -> {OBSTRUCTION}")
    else:
        if OBSTRUCTION.exists():
            OBSTRUCTION.unlink()
        print(f"[P1456] status={status} checked={len(rows)}")


if __name__ == "__main__":
    main()
