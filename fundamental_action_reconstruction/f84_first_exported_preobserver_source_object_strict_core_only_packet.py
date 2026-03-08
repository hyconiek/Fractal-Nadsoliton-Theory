#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f84_first_exported_preobserver_source_object_strict_core_only_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f81 = load_json(
        "fundamental_action_reconstruction/generated/f81_first_additive_preobserver_strict_core_source_object_export_packet_summary.json"
    )

    props = f81["strict_core_export_properties"]
    carrier = f81["carrier"]

    summary = {
        "stage": "F84",
        "lane": "first_exported_preobserver_source_object_strict_core_only_only",
        "status": "F84_EXECUTED_FIRST_EXPORTED_PREOBSERVER_SOURCE_OBJECT_STRICT_CORE_ONLY_PACKET_NO_FALSE_PASS",
        "exported_object": "S_preLM_strict_core_source_object_v1",
        "strict_core_only": props["strict_core_only"],
        "uses_axiom_lane_artifact": props["uses_axiom_lane_artifact"],
        "observer_excluded": carrier["observer_excluded"],
        "full_admissibility_pass": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
