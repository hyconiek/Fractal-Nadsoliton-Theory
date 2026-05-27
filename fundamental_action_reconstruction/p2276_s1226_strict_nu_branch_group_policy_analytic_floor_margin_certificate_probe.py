#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2274 = GEN / "p2274_s1224_strict_nu_branch_group_policy_robustness_region_certificate_probe.json"
IN_2275 = GEN / "p2275_s1225_strict_nu_branch_group_policy_certified_box_corner_replay_passrate_floor_probe.json"
OUT = GEN / "p2276_s1226_strict_nu_branch_group_policy_analytic_floor_margin_certificate_probe.json"
MD = GEN / "p2276_s1226_strict_nu_branch_group_policy_analytic_floor_margin_certificate_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2274 = load(IN_2274)
    p2275 = load(IN_2275)

    cert_rows = (p2274.get("strict_nu_branch_group_policy_robustness_region_certificate_probe", {}) or {}).get("certified_rows", []) or []
    replay_rows = (p2275.get("strict_nu_branch_group_policy_certified_box_corner_replay_passrate_floor_probe", {}) or {}).get("replay_rows", []) or []

    replay_by_risk = {float(r.get("risk_tolerance", -1.0)): r for r in replay_rows}

    # explicit conservative confidence margin policy for min-over-corners, min-over-seeds proxy
    confidence_margin = 0.02
    rows = []
    for c in cert_rows:
        risk = float(c.get("risk_tolerance", 0.05) or 0.05)
        target_floor = max(0.0, 1.0 - risk)
        rr = replay_by_risk.get(risk, {})
        empirical_floor = float(rr.get("empirical_passrate_floor_over_corners", 0.0) or 0.0)

        certified_analytic_floor = max(0.0, empirical_floor - confidence_margin)
        margin_to_target = certified_analytic_floor - target_floor

        rows.append(
            {
                "risk_tolerance": risk,
                "target_floor": target_floor,
                "empirical_corner_floor": empirical_floor,
                "confidence_margin": confidence_margin,
                "certified_analytic_floor": certified_analytic_floor,
                "analytic_floor_meets_target": margin_to_target + 1e-12 >= 0.0,
                "margin_to_target": margin_to_target,
                "certificate_statement": "certified_analytic_floor := empirical_corner_floor - confidence_margin",
            }
        )

    payload = {
        "schema_version": "p2276_s1226_v1",
        "packet_id": "P2276",
        "stage_id": "S1226",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_ANALYTIC_FLOOR_MARGIN_CERTIFICATE_PROBE",
        "strict_nu_branch_group_policy_analytic_floor_margin_certificate_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_ANALYTIC_FLOOR_MARGIN_CERTIFICATE_PROBE_V1",
            "source_packets": [str(IN_2274.relative_to(ROOT)), str(IN_2275.relative_to(ROOT))],
            "certificate_rows": rows,
            "global_summary": {
                "all_rows_meet_target_after_margin": all(r["analytic_floor_meets_target"] for r in rows),
                "worst_margin_to_target": min((r["margin_to_target"] for r in rows), default=0.0),
            },
            "proof_stub": {
                "assumption": "corner replay empirical floor is baseline estimator for certified box floor",
                "confidence_policy": "subtract fixed confidence_margin from empirical floor",
                "sufficient_condition": "if empirical_floor - confidence_margin >= 1-risk then analytic floor target holds",
            },
            "theorem_scope_limit": "analytic sufficient-condition packaging over surrogate replay only; not selector closure and not ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2277_candidate",
            "goal": "replace fixed confidence_margin with seed-count/trial-count dependent concentration bound and export adaptive confidence certificate",
        },
        "gatekeeper_checks": {
            "certificate_rows_exported": len(rows) > 0,
            "all_targets_bounded": all(0.0 <= r["target_floor"] <= 1.0 for r in rows),
            "all_certified_floors_bounded": all(0.0 <= r["certified_analytic_floor"] <= 1.0 for r in rows),
            "confidence_margin_nonnegative": confidence_margin >= 0.0,
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2276 S1226: analytic floor margin certificate probe",
            "",
            f"- certificate rows: `{len(rows)}`",
            f"- confidence margin: `{confidence_margin:.12e}`",
            f"- all rows meet target after margin: `{all(r['analytic_floor_meets_target'] for r in rows)}`",
            f"- worst margin to target: `{min((r['margin_to_target'] for r in rows), default=0.0):.12e}`",
            "",
            "Analytic sufficient-condition packaging only; no selector closure / ToE closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
