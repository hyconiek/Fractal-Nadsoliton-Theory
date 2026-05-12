#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

FAR = Path(__file__).resolve().parent
GEN = FAR / "generated"


def _load_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _candidate_entry(claim_id: str, residual_file: str | None, uncertainty_budget: dict) -> dict:
    entry = {
        "claim_id": claim_id,
        "residual_file": residual_file,
        "uncertainty_budget": uncertainty_budget,
        "residual_status": "MISSING",
        "max_abs_z": None,
        "strict_upgrade_ready": False,
        "notes": "",
    }
    if residual_file is None:
        entry["notes"] = "No residual artifact declared yet."
        return entry

    rpath = GEN / Path(residual_file).name
    if not rpath.exists():
        entry["notes"] = "Residual artifact path missing on current repo state."
        return entry

    data = _load_json(rpath)
    rows = data.get("rows", [])
    max_abs_z = max((float(r.get("abs_z", 0.0)) for r in rows), default=0.0)
    status = data.get("pass_fail_status_v1", "UNKNOWN")

    entry["residual_status"] = status
    entry["max_abs_z"] = max_abs_z
    entry["strict_upgrade_ready"] = status == "PASS" and uncertainty_budget.get("declared") is True
    entry["notes"] = "Residual evaluated from artifact." if rows else "No rows in residual artifact."
    return entry


def main() -> None:
    score = _load_json(GEN / "p1361_far_constant_claims_scoreboard_summary.json")
    claims = score.get("claims", [])

    strict_candidates = [c for c in claims if c.get("status") == "strict_candidate"]

    # Minimal explicit uncertainty declarations (placeholder discipline, not tuned values).
    budgets = {
        "C1_gauge_couplings_g_gp_g3": {"declared": True, "type": "single_scale_sigma_budget_v1", "components": ["numerical_discretization", "mapping_model_error", "reference_dataset_sigma"]},
        "C3_fine_structure_successor": {"declared": False, "type": "missing", "components": []},
        "C5_kernel_only_first_prediction_run": {"declared": True, "type": "single_scale_sigma_budget_v1", "components": ["numerical_discretization", "kernel_parameter_lock", "reference_dataset_sigma"]},
    }

    benchmark = []
    for c in strict_candidates:
        cid = c["id"]
        residual_file = c.get("residual_artifact")
        benchmark.append(_candidate_entry(cid, residual_file, budgets.get(cid, {"declared": False, "type": "missing", "components": []})))

    summary = {
        "packet": "P1362",
        "as_of": "2026-05-12",
        "strict_candidates_benchmark": benchmark,
        "ready_count": sum(1 for b in benchmark if b["strict_upgrade_ready"]),
        "candidate_count": len(benchmark),
        "next_priority": "P1363_UPGRADE_ONLY_READY_CANDIDATES_AND_OPEN_BLOCKERS_FOR_THE_REST",
    }

    out = GEN / "p1362_strict_candidate_residual_benchmark_summary.json"
    out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1362] wrote {out}; ready={summary['ready_count']}/{summary['candidate_count']}")


if __name__ == "__main__":
    main()
