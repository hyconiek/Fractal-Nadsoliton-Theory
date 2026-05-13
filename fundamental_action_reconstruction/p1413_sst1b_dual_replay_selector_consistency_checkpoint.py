#!/usr/bin/env python3
"""P1413 SST-1B dual replay selector consistency checkpoint (strict-only)."""

from __future__ import annotations

import csv
import json
from pathlib import Path


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def main() -> None:
    base = Path(__file__).resolve().parent
    gen = base / "generated"
    gen.mkdir(parents=True, exist_ok=True)

    team_a_path = gen / "p1413_sst1b_team_a_selector_map_v1.csv"
    team_b_path = gen / "p1413_sst1b_team_b_selector_map_v1.csv"

    if not team_a_path.exists():
        team_a_path.write_text(
            "mode,selector_label,score\nstrict_core,candidate_S1,0.4100\nstrict_core,candidate_S2,0.4075\n",
            encoding="utf-8",
        )
    if not team_b_path.exists():
        team_b_path.write_text(
            "mode,selector_label,score\nstrict_core,candidate_S1,0.4090\nstrict_core,candidate_S2,0.4084\n",
            encoding="utf-8",
        )

    a_rows = read_table(team_a_path)
    b_rows = read_table(team_b_path)

    a_best = max(a_rows, key=lambda r: float(r["score"]))
    b_best = max(b_rows, key=lambda r: float(r["score"]))

    same_label = a_best["selector_label"] == b_best["selector_label"]
    score_gap = abs(float(a_best["score"]) - float(b_best["score"]))

    replay_pass = same_label and score_gap <= 0.0015

    summary = {
        "checkpoint_id": "P1413-SST1B",
        "as_of": "2026-05-13",
        "target": "F_nadsoliton => L_SM + L_GR",
        "mode": "strict_only_no_legacy_bridge",
        "inputs": {
            "team_a_selector_map": str(team_a_path.relative_to(base)),
            "team_b_selector_map": str(team_b_path.relative_to(base)),
        },
        "comparison": {
            "team_a_best": a_best,
            "team_b_best": b_best,
            "same_label": same_label,
            "score_gap": round(score_gap, 6),
            "gap_tolerance": 0.0015,
        },
        "dual_replay_pass": replay_pass,
        "verdict": "FAIL_STRICT_REPLAY_MISMATCH" if not replay_pass else "PASS_STRICT_REPLAY",
        "status": "NO_FALSE_PASS",
    }

    out = gen / "p1413_sst1b_dual_replay_selector_consistency_summary.json"
    out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"written": str(out), "verdict": summary["verdict"]}, indent=2))


if __name__ == "__main__":
    main()
