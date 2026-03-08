#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f85_first_exported_preobserver_source_object_kernel_split_safety_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f77 = load_json(
        "fundamental_action_reconstruction/generated/f77_first_additive_preobserver_source_object_admissibility_upgrade_target_packet_summary.json"
    )
    f81 = load_json(
        "fundamental_action_reconstruction/generated/f81_first_additive_preobserver_strict_core_source_object_export_packet_summary.json"
    )

    props = f81["strict_core_export_properties"]
    guardrails = f77["guardrails"]

    summary = {
        "stage": "F85",
        "lane": "first_exported_preobserver_source_object_kernel_split_safety_only",
        "status": "F85_EXECUTED_FIRST_EXPORTED_PREOBSERVER_SOURCE_OBJECT_KERNEL_SPLIT_SAFETY_PACKET_NO_FALSE_PASS",
        "exported_object": "S_preLM_strict_core_source_object_v1",
        "kernel_split_safe": props["kernel_split_safe"],
        "no_external_selector_import": props["no_external_selector_import"],
        "guardrail_kernel_split_safe": guardrails["kernel_split_safe"],
        "full_admissibility_pass": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
