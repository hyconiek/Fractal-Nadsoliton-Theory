#!/usr/bin/env python3
"""
QW-2077: Empirical prediction validation gate.

Reads preregistered predictions (QW-2076) and external observation input,
then produces explicit support/falsification/pending statuses.
"""

from __future__ import annotations

import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple


ROOT = Path(__file__).resolve().parent
PREREG_JSON = ROOT / "report_qw2076_empirical_prediction_preregistration.json"
DEFAULT_OBS_INPUT = ROOT / "empirical_observations_input_qw2077.json"
OUT_JSON = ROOT / "report_qw2077_empirical_prediction_validation_gate.json"
OUT_MD = ROOT / "RAPORT_QW2077_EMPIRICAL_PREDICTION_VALIDATION_GATE.md"


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def get_prediction(preds: List[Dict], pred_id: str) -> Dict:
    for p in preds:
        if p.get("id") == pred_id:
            return p
    raise KeyError(f"Missing prediction id: {pred_id}")


def validate_pmns(pred: Dict, obs: Dict) -> Dict:
    needed = ["sin_delta_central", "sin_delta_ci95_low", "sin_delta_ci95_high"]
    if any(obs.get(k) is None for k in needed):
        return {"status": "pending_data", "reason": "missing_pmns_input"}

    c = float(obs["sin_delta_central"])
    lo = float(obs["sin_delta_ci95_low"])
    hi = float(obs["sin_delta_ci95_high"])
    b_lo = float(pred["pre_registered_band"]["sin_delta_min"])
    b_hi = float(pred["pre_registered_band"]["sin_delta_max"])

    ci_outside = lo > b_hi or hi < b_lo
    central_inside = b_lo <= c <= b_hi

    if ci_outside:
        status = "falsified"
    elif central_inside:
        status = "supported"
    else:
        status = "inconclusive"

    return {
        "status": status,
        "central": c,
        "ci95": [lo, hi],
        "band": [b_lo, b_hi],
        "source": obs.get("source", ""),
    }


def validate_cosmo(pred: Dict, obs_nodes: List[Dict]) -> Dict:
    if not isinstance(obs_nodes, list) or len(obs_nodes) == 0:
        return {"status": "pending_data", "reason": "missing_cosmology_nodes"}

    pred_map = {float(x["z"]): float(x["w_eff_pred"]) for x in pred["predicted_nodes"]}
    abs_corridor = float(pred["absolute_corridor"])

    node_checks = []
    missing_nodes = []
    for z, w_pred in pred_map.items():
        cand = [x for x in obs_nodes if x.get("z") is not None and abs(float(x["z"]) - z) < 1e-9]
        if not cand:
            missing_nodes.append(z)
            continue
        node = cand[0]
        w_obs = node.get("w_eff_central")
        sigma = node.get("sigma_total")
        if w_obs is None or sigma is None:
            node_checks.append({"z": z, "status": "pending_data"})
            continue
        w_obs = float(w_obs)
        sigma = abs(float(sigma))
        tol = max(abs_corridor, 3.0 * sigma)
        resid = abs(w_obs - w_pred)
        ok = resid <= tol
        node_checks.append(
            {
                "z": z,
                "w_pred": w_pred,
                "w_obs": w_obs,
                "sigma_total": sigma,
                "tolerance": tol,
                "residual_abs": resid,
                "pass": bool(ok),
            }
        )

    if missing_nodes:
        return {"status": "pending_data", "reason": "missing_z_nodes", "missing_nodes": missing_nodes, "nodes": node_checks}
    if any(n.get("status") == "pending_data" for n in node_checks):
        return {"status": "pending_data", "reason": "incomplete_node_values", "nodes": node_checks}

    all_pass = all(bool(n.get("pass")) for n in node_checks)
    return {
        "status": "supported" if all_pass else "falsified",
        "nodes": node_checks,
    }


def validate_gw(pred: Dict, obs: Dict) -> Dict:
    needed = [
        "auc_h1l1_vs_ctrl",
        "adv_shared_minus_ctrl_q90",
        "sep_median_h1l1_minus_ctrl",
        "control_median_gap",
    ]
    if any(obs.get(k) is None for k in needed):
        return {"status": "pending_data", "reason": "missing_gw_metrics"}

    th = pred["hard_thresholds"]
    auc = float(obs["auc_h1l1_vs_ctrl"])
    adv = float(obs["adv_shared_minus_ctrl_q90"])
    sep = float(obs["sep_median_h1l1_minus_ctrl"])
    gap = float(obs["control_median_gap"])

    checks = {
        "auc_ge_min": auc >= float(th["auc_h1l1_vs_ctrl_min"]),
        "adv_ge_min": adv >= float(th["adv_shared_minus_ctrl_q90_min"]),
        "sep_ge_min": sep >= float(th["sep_median_h1l1_minus_ctrl_min"]),
        "gap_le_max": gap <= float(th["control_median_gap_max"]),
    }
    all_pass = all(checks.values())
    return {
        "status": "supported" if all_pass else "falsified",
        "metrics": {
            "auc_h1l1_vs_ctrl": auc,
            "adv_shared_minus_ctrl_q90": adv,
            "sep_median_h1l1_minus_ctrl": sep,
            "control_median_gap": gap,
            "n_windows": obs.get("n_windows"),
            "source": obs.get("source", ""),
        },
        "checks": checks,
    }


def main() -> None:
    obs_path = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else DEFAULT_OBS_INPUT.resolve()

    prereg = load_json(PREREG_JSON)
    preds = prereg["predictions"]

    if not obs_path.exists():
        out = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "verdict": "EMPIRICAL_VALIDATION_PENDING_NO_OBSERVATION_FILE",
            "preregistration": PREREG_JSON.name,
            "observation_input": str(obs_path),
            "required_next_step": "PROVIDE_OBSERVATION_INPUT_AND_RERUN_QW2077",
        }
        OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
        OUT_MD.write_text(
            "\n".join(
                [
                    "# RAPORT QW-2077: EMPIRICAL PREDICTION VALIDATION GATE",
                    "",
                    f"- Data UTC: {out['generated_utc']}",
                    f"- Verdict: **{out['verdict']}**",
                    f"- Missing observation file: `{obs_path}`",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        print(f"[QW-2077] verdict={out['verdict']}")
        return

    obs = load_json(obs_path)
    pmns_pred = get_prediction(preds, "PRED-2076-PMNS-CP")
    cosmo_pred = get_prediction(preds, "PRED-2076-COSMO-WEFF")
    gw_pred = get_prediction(preds, "PRED-2076-GW-HOLDOUT")

    pmns_res = validate_pmns(pmns_pred, obs.get("pmns_cp", {}))
    cosmo_res = validate_cosmo(cosmo_pred, obs.get("cosmology_weff_nodes", []))
    gw_res = validate_gw(gw_pred, obs.get("gw_future_holdout", {}))

    statuses = [pmns_res["status"], cosmo_res["status"], gw_res["status"]]
    if any(s == "falsified" for s in statuses):
        verdict = "EMPIRICAL_VALIDATION_PARTIAL_FALSIFICATION"
    elif all(s == "supported" for s in statuses):
        verdict = "EMPIRICAL_VALIDATION_SUPPORT_ALL_REGISTERED_PREDICTIONS"
    elif all(s == "pending_data" for s in statuses):
        verdict = "EMPIRICAL_VALIDATION_PENDING_DATA"
    else:
        verdict = "EMPIRICAL_VALIDATION_MIXED_OR_INCONCLUSIVE"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "preregistration": PREREG_JSON.name,
        "observation_input": str(obs_path),
        "results": {
            "PRED-2076-PMNS-CP": pmns_res,
            "PRED-2076-COSMO-WEFF": cosmo_res,
            "PRED-2076-GW-HOLDOUT": gw_res,
        },
        "status_counts": {
            "supported": sum(1 for s in statuses if s == "supported"),
            "falsified": sum(1 for s in statuses if s == "falsified"),
            "pending_data": sum(1 for s in statuses if s == "pending_data"),
            "inconclusive": sum(1 for s in statuses if s == "inconclusive"),
        },
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2077: EMPIRICAL PREDICTION VALIDATION GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- observation input: `{obs_path}`",
        "",
        "## Status Counts",
        f"- supported: {out['status_counts']['supported']}",
        f"- falsified: {out['status_counts']['falsified']}",
        f"- pending_data: {out['status_counts']['pending_data']}",
        f"- inconclusive: {out['status_counts']['inconclusive']}",
        "",
        "## Per Prediction",
        f"- PRED-2076-PMNS-CP: {pmns_res['status']}",
        f"- PRED-2076-COSMO-WEFF: {cosmo_res['status']}",
        f"- PRED-2076-GW-HOLDOUT: {gw_res['status']}",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2077] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2077] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2077] verdict={verdict}")


if __name__ == "__main__":
    main()
