#!/usr/bin/env python3
"""
QW-2083: Missing-14 epistemic status gate.

Purpose:
- consume QW-2081 strict-rigor frontier and produce deterministic update rows
  for QW-2069 package integration,
- explicitly separate "addressed with non-closing epistemic label" from
  "still strictly missing".
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
IN_2081 = ROOT / "report_qw2081_missing14_strict_rigor_frontier.json"
OUT_JSON = ROOT / "report_qw2083_missing14_epistemic_status_gate.json"
OUT_MD = ROOT / "RAPORT_QW2083_MISSING14_EPISTEMIC_STATUS_GATE.md"


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    r2081 = load_json(IN_2081)
    frontier: List[Dict] = r2081["frontier"]

    updates: List[Dict] = []
    n_still_missing = 0
    n_addressed_nonclosing = 0
    n_strict_target_miss = 0

    for row in frontier:
        pid = row["id"]
        cls = row["classification"]
        pred_val = row.get("proposed_value")
        method = row.get("method", "qw2081_frontier_map")
        notes = row.get("notes", "")

        if cls == "strict_candidate_target_miss":
            n_strict_target_miss += 1
            updates.append(
                {
                    "id": pid,
                    "predicted_value": pred_val,
                    "method": method,
                    "status": "derived_strict_target_miss",
                    "strict_level": "strict_internal_gate",
                    "notes": notes,
                }
            )
            continue

        if cls == "anchor_dependent_baseline_only":
            n_addressed_nonclosing += 1
            updates.append(
                {
                    "id": pid,
                    "predicted_value": pred_val,
                    "method": method,
                    "status": "derived_nofit_anchor_dependent",
                    "strict_level": "physical_relation_anchor_dependent",
                    "notes": notes,
                }
            )
            continue

        if cls == "model_assumption_required":
            n_addressed_nonclosing += 1
            updates.append(
                {
                    "id": pid,
                    "predicted_value": pred_val,
                    "method": method,
                    "status": "derived_model_assumption_nonclosing",
                    "strict_level": "model_assumption_anchor",
                    "notes": notes,
                }
            )
            continue

        if cls == "coupled_underdetermined_pair":
            n_addressed_nonclosing += 1
            updates.append(
                {
                    "id": pid,
                    "predicted_value": pred_val,
                    "method": method,
                    "status": "derived_coupled_anchor_dependent",
                    "strict_level": "coupled_anchor_dependent",
                    "notes": notes,
                }
            )
            continue

        if cls == "underdetermined_without_new_observable":
            n_still_missing += 1
            updates.append(
                {
                    "id": pid,
                    "predicted_value": None,
                    "method": method,
                    "status": "missing",
                    "strict_level": "not_derived",
                    "notes": notes,
                }
            )
            continue

        raise RuntimeError(f"Unhandled QW-2081 class: {cls} (id={pid})")

    if len(updates) == 0:
        verdict = "MISSING14_EPISTEMIC_STATUS_GATE_PASS_NOTHING_TO_UPDATE"
    elif n_still_missing == 0 and n_strict_target_miss == 0:
        verdict = "MISSING14_EPISTEMIC_STATUS_GATE_PASS_ALL_MAPPED_NONCLOSING_ALLOWED"
    elif n_still_missing == 0:
        verdict = "MISSING14_EPISTEMIC_STATUS_GATE_PASS_WITH_TARGET_MISS"
    else:
        verdict = "MISSING14_EPISTEMIC_STATUS_GATE_PARTIAL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_frontier": IN_2081.name,
        "input_verdict": r2081.get("verdict"),
        "updates": updates,
        "summary": {
            "n_updates": len(updates),
            "n_strict_target_miss": n_strict_target_miss,
            "n_addressed_nonclosing": n_addressed_nonclosing,
            "n_still_missing": n_still_missing,
        },
        "verdict": verdict,
        "required_next_step": "RERUN_QW2069_QW2070_QW2071_WITH_QW2083_UPDATES",
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2083: MISSING-14 EPISTEMIC STATUS GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Input frontier: `{IN_2081.name}`",
        f"- Verdict: **{out['verdict']}**",
        "",
        "## Summary",
        f"- n_updates: {out['summary']['n_updates']}",
        f"- n_strict_target_miss: {n_strict_target_miss}",
        f"- n_addressed_nonclosing: {n_addressed_nonclosing}",
        f"- n_still_missing: {n_still_missing}",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2083] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2083] Saved MD:   {OUT_MD.name}")
    print(
        f"[QW-2083] verdict={out['verdict']} "
        f"updates={len(updates)} still_missing={n_still_missing}"
    )


if __name__ == "__main__":
    main()
