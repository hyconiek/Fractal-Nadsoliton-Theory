#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f83_first_exported_preobserver_source_object_source_seed_only_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f81 = load_json(
        "fundamental_action_reconstruction/generated/f81_first_additive_preobserver_strict_core_source_object_export_packet_summary.json"
    )
    f82 = load_json(
        "fundamental_action_reconstruction/generated/f82_first_exported_preobserver_source_object_second_clause_typing_packet_summary.json"
    )

    summary = {
        "stage": "F83",
        "lane": "third_admissibility_clause_source_seed_only_only",
        "goal": "freeze_the_source_seed_only_status_of_the_exported_preobserver_source_object",
        "status": "F83_EXECUTED_FIRST_EXPORTED_PREOBSERVER_SOURCE_OBJECT_SOURCE_SEED_ONLY_PACKET_NO_FALSE_PASS",
        "exported_object": "S_preLM_strict_core_source_object_v1",
        "source_seed_only": f81["strict_core_export_properties"]["source_seed_only"],
        "E_orient_exported": f82["E_orient_exported"],
        "B_sel_exported": False,
        "R_sel_exported": False,
        "O_sel_exported": False,
        "full_admissibility_pass": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
