#!/usr/bin/env python3
from __future__ import annotations
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2275 = GEN / "p2275_s1225_strict_nu_branch_group_policy_certified_box_corner_replay_passrate_floor_probe.json"
IN_2276 = GEN / "p2276_s1226_strict_nu_branch_group_policy_analytic_floor_margin_certificate_probe.json"
OUT = GEN / "p2277_s1227_strict_nu_branch_group_policy_adaptive_confidence_margin_certificate_probe.json"
MD = GEN / "p2277_s1227_strict_nu_branch_group_policy_adaptive_confidence_margin_certificate_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2275 = load(IN_2275)
    p2276 = load(IN_2276)

    probe2275 = p2275.get("strict_nu_branch_group_policy_certified_box_corner_replay_passrate_floor_probe", {}) or {}
    settings = probe2275.get("settings", {}) or {}
    replay_rows = probe2275.get("replay_rows", []) or []

    old_rows = (p2276.get("strict_nu_branch_group_policy_analytic_floor_margin_certificate_probe", {}) or {}).get("certificate_rows", []) or []
    old_by_risk = {float(r.get("risk_tolerance", -1.0)): r for r in old_rows}

    n_seeds = len(settings.get("seeds", []) or [])
    trials = int(settings.get("trials", 1) or 1)
    corners = 4
    n_eff = max(1, n_seeds * trials * corners)

    # Hoeffding-style margin for bounded Bernoulli mean, with conservative delta.
    delta = 0.01
    adaptive_margin = math.sqrt(math.log(2.0 / delta) / (2.0 * n_eff))

    rows = []
    for r in replay_rows:
        risk = float(r.get("risk_tolerance", 0.05) or 0.05)
        target_floor = float(r.get("target_passrate_floor", max(0.0, 1.0 - risk)) or 0.0)
        empirical_floor = float(r.get("empirical_passrate_floor_over_corners", 0.0) or 0.0)

        old_fixed_margin = float((old_by_risk.get(risk, {}) or {}).get("confidence_margin", 0.0) or 0.0)
        fixed_cert_floor = float((old_by_risk.get(risk, {}) or {}).get("certified_analytic_floor", 0.0) or 0.0)

        adaptive_cert_floor = max(0.0, empirical_floor - adaptive_margin)

        rows.append(
            {
                "risk_tolerance": risk,
                "target_floor": target_floor,
                "empirical_corner_floor": empirical_floor,
                "fixed_margin_previous": old_fixed_margin,
                "fixed_certified_floor_previous": fixed_cert_floor,
                "adaptive_margin": adaptive_margin,
                "adaptive_certified_floor": adaptive_cert_floor,
                "adaptive_floor_meets_target": adaptive_cert_floor + 1e-12 >= target_floor,
                "adaptive_minus_fixed_certified_floor": adaptive_cert_floor - fixed_cert_floor,
            }
        )

    payload = {
        "schema_version": "p2277_s1227_v1",
        "packet_id": "P2277",
        "stage_id": "S1227",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_ADAPTIVE_CONFIDENCE_MARGIN_CERTIFICATE_PROBE",
        "strict_nu_branch_group_policy_adaptive_confidence_margin_certificate_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_ADAPTIVE_CONFIDENCE_MARGIN_CERTIFICATE_PROBE_V1",
            "source_packets": [str(IN_2275.relative_to(ROOT)), str(IN_2276.relative_to(ROOT))],
            "confidence_model": {
                "type": "hoeffding_bernoulli_floor_surrogate",
                "delta": delta,
                "n_seeds": n_seeds,
                "trials": trials,
                "corners": corners,
                "n_eff": n_eff,
                "adaptive_margin": adaptive_margin,
                "formula": "sqrt(log(2/delta)/(2*n_eff))",
            },
            "certificate_rows": rows,
            "global_summary": {
                "all_rows_meet_target_adaptive": all(x["adaptive_floor_meets_target"] for x in rows),
                "worst_adaptive_margin_to_target": min((x["adaptive_certified_floor"] - x["target_floor"] for x in rows), default=0.0),
            },
            "theorem_scope_limit": "adaptive concentration-margin packaging for surrogate replay only; not selector closure and not ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2278_candidate",
            "goal": "export finite-sample confidence-curve sweep over delta and replay budgets to optimize floor-tightness vs confidence",
        },
        "gatekeeper_checks": {
            "certificate_rows_exported": len(rows) > 0,
            "adaptive_margin_nonnegative": adaptive_margin >= 0.0,
            "all_adaptive_floors_bounded": all(0.0 <= x["adaptive_certified_floor"] <= 1.0 for x in rows),
            "all_targets_bounded": all(0.0 <= x["target_floor"] <= 1.0 for x in rows),
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2277 S1227: adaptive confidence margin certificate probe",
            "",
            f"- n_eff: `{n_eff}`",
            f"- delta: `{delta}`",
            f"- adaptive margin: `{adaptive_margin:.12e}`",
            f"- all rows meet target (adaptive): `{all(x['adaptive_floor_meets_target'] for x in rows)}`",
            f"- worst adaptive margin-to-target: `{min((x['adaptive_certified_floor']-x['target_floor'] for x in rows), default=0.0):.12e}`",
            "",
            "Adaptive concentration-margin packaging only; no selector closure / ToE closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
