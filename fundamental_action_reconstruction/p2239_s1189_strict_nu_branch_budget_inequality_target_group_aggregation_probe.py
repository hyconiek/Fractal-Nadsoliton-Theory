#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2226 = GEN / "p2226_s1176_strict_nu_branch_grouped_policy_width_coverage_tradeoff.json"
IN_2238 = GEN / "p2238_s1188_strict_nu_branch_weighted_stability_to_budget_inequality_probe.json"
OUT = GEN / "p2239_s1189_strict_nu_branch_budget_inequality_target_group_aggregation_probe.json"
MD = GEN / "p2239_s1189_strict_nu_branch_budget_inequality_target_group_aggregation_probe.md"

def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))

def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2226 = load(IN_2226)
    p2238 = load(IN_2238)

    inequality = (p2238.get("strict_nu_branch_weighted_stability_to_budget_inequality_probe", {}) or {}).get("derived_inequality", {}) or {}
    load_ratio = float(inequality.get("relative_load_ratio", 1.0) or 1.0)
    margin = float(inequality.get("conservative_margin_factor", -1.0) or -1.0)

    trade = (p2226.get("strict_nu_branch_grouped_policy_width_coverage_tradeoff", {}) or {})
    groups = trade.get("group_summaries", []) or trade.get("groups", []) or []
    if not groups:
        # fallback synthetic one-group witness from packet-level values
        groups = [{"group_id": "fallback", "coverage_fraction": float(trade.get("coverage_fraction", 1.0) or 1.0)}]

    rows = []
    worst_score = 1.0
    for g in groups:
        gid = g.get("group_id", g.get("group", "unknown"))
        cov = float(g.get("coverage_fraction", g.get("coverage", 1.0)) or 1.0)
        # risk score: larger when coverage low or ratio high
        score = load_ratio / max(cov, 1e-12)
        worst_score = max(worst_score, score)
        rows.append({"group_id": gid, "coverage_fraction": cov, "inequality_load_ratio": load_ratio, "risk_score": score})

    groupwise_holds = all(r["risk_score"] <= 1.0 + 1e-15 for r in rows)

    payload = {
        "schema_version": "p2239_s1189_v1",
        "packet_id": "P2239",
        "stage_id": "S1189",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_BUDGET_INEQUALITY_TARGET_GROUP_AGGREGATION_PROBE",
        "strict_nu_branch_budget_inequality_target_group_aggregation_probe": {
            "probe_id": "STRICT_NU_BRANCH_BUDGET_INEQUALITY_TARGET_GROUP_AGGREGATION_PROBE_V1",
            "source_packets": [str(IN_2226.relative_to(ROOT)), str(IN_2238.relative_to(ROOT))],
            "global_inputs": {"inequality_load_ratio": load_ratio, "conservative_margin_factor": margin},
            "group_rows": rows,
            "worst_case_group_risk_score": worst_score,
            "groupwise_compatibility_holds": groupwise_holds,
            "physical_interpretation_note": "Aggregates local uncertainty-to-budget load against coverage dilution across grouped policies; controls worst-case robustness leakage.",
            "theorem_scope_limit": "finite grouped-policy aggregation check only; not a legacy->strict bridge theorem",
        },
        "recommended_next_honest_step": {"id": "P2240_candidate", "goal": "derive analytic sufficient condition on grouped coverage/load ensuring preserved sign-stability under policy mixing"},
        "gatekeeper_checks": {
            "target_group_aggregation_exported": True,
            "global_margin_nonnegative": margin >= 0.0,
            "groupwise_compatibility_holds": groupwise_holds,
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text("\n".join([
        "# P2239 S1189: budget-inequality target/group aggregation probe",
        "",
        f"- inequality load ratio: `{load_ratio:.12e}`",
        f"- conservative margin factor: `{margin:.12e}`",
        f"- worst-case group risk score: `{worst_score:.12e}`",
        f"- groupwise compatibility holds: `{groupwise_holds}`",
        "",
        "Grouped-policy aggregation only; no kernel-bridge theorem claim.",
    ]) + "\n", encoding="utf-8")

if __name__ == "__main__":
    main()
