#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f79_first_additive_preobserver_source_object_nonreduction_comparison_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f76 = load_json(
        "fundamental_action_reconstruction/generated/f76_first_additive_preobserver_source_object_construction_attempt_packet_summary.json"
    )
    candidate = f76["attempt_profile"]["closed_form"]
    cos_phi = f76["generator"]["matrix"][1][0]
    i_mat = f76["strict_kernel_role"]["I_mat(d_star)"]

    packaged = [
        cos_phi,
        cos_phi,
        cos_phi * i_mat,
    ]
    delta = [
        candidate[0] - packaged[0],
        candidate[1] - packaged[1],
        candidate[2] - packaged[2],
    ]
    delta_norm = math.sqrt(sum(x * x for x in delta))

    summary = {
        "stage": "F79",
        "lane": "first_additive_preobserver_source_object_nonreduction_comparison_only",
        "goal": "compare_the_additive_attempt_against_the_same_basis_packaged_realization_of_the_F75_target",
        "status": "F79_EXECUTED_FIRST_ADDITIVE_PREOBSERVER_SOURCE_OBJECT_NONREDUCTION_COMPARISON_PACKET_NO_FALSE_PASS",
        "construction_attempt": "S_preLM_additive_candidate_v1",
        "packaged_realization": "S_preLM_target_packaged_realization_v1(d_*=1)",
        "same_basis": ["u_T", "u_L", "u_M"],
        "candidate_closed_form": candidate,
        "packaged_closed_form": packaged,
        "delta": delta,
        "delta_norm": delta_norm,
        "delta_is_nonzero": delta_norm > 0.0,
        "future_only": True,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
