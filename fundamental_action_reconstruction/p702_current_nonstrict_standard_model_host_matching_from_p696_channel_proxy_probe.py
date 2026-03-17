#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linear_sum_assignment


AS_OF = "2026-03-17"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P696_SUMMARY = (
    GENERATED
    / "p696_current_strict_physical_computability_selector_aligned_channel_spectrum_proxy_from_projective_selector_closure_probe_summary.json"
)

IN_SM_TARGETS = ROOT / "external_data" / "sm_mass_targets_v1.json"
TEMPLATE_SM_TARGETS = ROOT / "external_data" / "SM_MASS_TARGETS_TEMPLATE.json"

OUT_JSON = (
    GENERATED
    / "p702_current_nonstrict_standard_model_host_matching_from_p696_channel_proxy_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p702_current_nonstrict_standard_model_host_matching_from_p696_channel_proxy_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_finite_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P696_SUMMARY, IN_SM_TARGETS]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P702",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "hints": {
                "sm_targets_template": str(TEMPLATE_SM_TARGETS.relative_to(REPO)),
                "sm_targets_required_path": str(IN_SM_TARGETS.relative_to(REPO)),
            },
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p696 = load_json(IN_P696_SUMMARY)
    channels_sorted = p696.get("channel_m2_proxy_sorted_ascending")
    if not (
        isinstance(channels_sorted, list)
        and len(channels_sorted) >= 1
        and all(isinstance(x, dict) and isinstance(x.get("channel"), str) and is_finite_number(x.get("m2")) for x in channels_sorted)
    ):
        artifact = {
            "stage": "P702",
            "status": "INVALID_P696_SUMMARY_SHAPE",
            "as_of": AS_OF,
            "error": "P696 summary must export channel_m2_proxy_sorted_ascending as list[{channel:str,m2:number}]",
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    channel_to_m2 = {str(x["channel"]): float(x["m2"]) for x in channels_sorted}

    sm = load_json(IN_SM_TARGETS)
    if not isinstance(sm, dict):
        artifact = {
            "stage": "P702",
            "status": "INVALID_SM_TARGETS_SHAPE",
            "as_of": AS_OF,
            "error": "SM targets file must be a JSON object",
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    units = sm.get("units")
    targets = sm.get("targets")
    predicted_channels = sm.get("predicted_channels")
    policy = sm.get("matching_policy") or {}

    if not isinstance(units, str) or not units:
        artifact = {
            "stage": "P702",
            "status": "INVALID_SM_TARGETS_UNITS",
            "as_of": AS_OF,
            "error": "sm_mass_targets_v1.json must define nonempty string units (e.g. 'GeV')",
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    if predicted_channels is None:
        predicted_channels = [str(x["channel"]) for x in channels_sorted]
    if not (isinstance(predicted_channels, list) and all(isinstance(c, str) and c for c in predicted_channels)):
        artifact = {
            "stage": "P702",
            "status": "INVALID_SM_TARGETS_PREDICTED_CHANNELS",
            "as_of": AS_OF,
            "error": "sm_mass_targets_v1.json predicted_channels must be a list[str] (or omit it to use P696 sorted channels)",
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    missing_channels = [c for c in predicted_channels if c not in channel_to_m2]
    if missing_channels:
        artifact = {
            "stage": "P702",
            "status": "NOT_COMPUTABLE_MISSING_CHANNELS_IN_P696",
            "as_of": AS_OF,
            "missing_channels": missing_channels,
            "available_channels": sorted(channel_to_m2.keys()),
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    if not (isinstance(targets, list) and all(isinstance(t, dict) for t in targets)):
        artifact = {
            "stage": "P702",
            "status": "INVALID_SM_TARGETS_TARGETS_LIST",
            "as_of": AS_OF,
            "error": "sm_mass_targets_v1.json must define targets as list[object]",
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    target_ids: list[str] = []
    target_mass: list[float] = []
    invalid_targets: list[str] = []
    for i, t in enumerate(targets):
        tid = t.get("id")
        m = t.get("mass")
        if not isinstance(tid, str) or not tid:
            invalid_targets.append(f"targets[{i}].id")
            continue
        if not (is_finite_number(m) and float(m) > 0.0):
            invalid_targets.append(f"targets[{i}].mass:{tid}")
            continue
        target_ids.append(tid)
        target_mass.append(float(m))

    if invalid_targets:
        artifact = {
            "stage": "P702",
            "status": "NOT_COMPUTABLE_INVALID_SM_TARGETS",
            "as_of": AS_OF,
            "invalid_fields": invalid_targets,
            "hint": f"Fill {IN_SM_TARGETS.relative_to(REPO)} using {TEMPLATE_SM_TARGETS.relative_to(REPO)} as a template.",
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    n_pred = len(predicted_channels)
    n_tgt = len(target_ids)
    if n_pred != n_tgt:
        artifact = {
            "stage": "P702",
            "status": "NOT_COMPUTABLE_TARGET_COUNT_MISMATCH",
            "as_of": AS_OF,
            "predicted_channel_count": n_pred,
            "target_count": n_tgt,
            "error": "Probe requires a bijection: predicted_channels and targets must have equal length.",
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    pred_m2 = np.array([float(channel_to_m2[c]) for c in predicted_channels], dtype=float)
    if not np.all(np.isfinite(pred_m2)) or not np.all(pred_m2 > 0.0):
        artifact = {
            "stage": "P702",
            "status": "NOT_COMPUTABLE_INVALID_P696_M2_VALUES",
            "as_of": AS_OF,
            "error": "All selected P696 m2 proxy values must be finite and >0.",
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    tgt_mass = np.array(target_mass, dtype=float)
    tgt_m2 = tgt_mass * tgt_mass

    log_pred = np.log(pred_m2)
    log_tgt = np.log(tgt_m2)
    c_pred = log_pred - float(np.mean(log_pred))
    c_tgt = log_tgt - float(np.mean(log_tgt))

    cost = (c_pred.reshape(-1, 1) - c_tgt.reshape(1, -1)) ** 2
    row_ind, col_ind = linear_sum_assignment(cost)

    # Global best-fit scale (in units^2) after assignment: log(scale) = mean(log(tgt_m2) - log(pred_m2)).
    log_scale = float(np.mean(log_tgt[col_ind] - log_pred[row_ind]))
    scale_m2 = float(math.exp(log_scale))
    scale_m = float(math.sqrt(scale_m2))

    pred_m2_scaled = scale_m2 * pred_m2
    pred_mass_scaled = np.sqrt(pred_m2_scaled)

    log_err_m2 = (np.log(pred_m2_scaled[row_ind]) - np.log(tgt_m2[col_ind])).astype(float)
    rmse_log_m2 = float(np.sqrt(np.mean(log_err_m2 * log_err_m2)))
    max_abs_log_m2 = float(np.max(np.abs(log_err_m2)))

    rel_err_mass = (pred_mass_scaled[row_ind] / tgt_mass[col_ind] - 1.0).astype(float)
    rmse_rel_mass = float(np.sqrt(np.mean(rel_err_mass * rel_err_mass)))
    max_abs_rel_mass = float(np.max(np.abs(rel_err_mass)))

    assignment = []
    for r, c in zip(row_ind.tolist(), col_ind.tolist(), strict=True):
        assignment.append(
            {
                "predicted_channel": predicted_channels[int(r)],
                "target_id": target_ids[int(c)],
                "proxy_m2": float(pred_m2[int(r)]),
                "target_mass": float(tgt_mass[int(c)]),
                "predicted_mass": float(pred_mass_scaled[int(r)]),
                "relative_error_mass": float(pred_mass_scaled[int(r)] / tgt_mass[int(c)] - 1.0),
                "log_error_m2": float(np.log(pred_m2_scaled[int(r)]) - np.log(tgt_m2[int(c)])),
            }
        )

    assignment_sorted_by_target = sorted(assignment, key=lambda a: str(a["target_id"]))

    pass_criteria = (policy.get("pass_criteria") or None) if isinstance(policy, dict) else None
    pass_claim_available = isinstance(pass_criteria, dict)
    pass_under_user_criteria = None
    failed_criteria: list[str] = []
    if pass_claim_available:
        pass_under_user_criteria = True
        if "max_abs_rel_mass" in pass_criteria:
            thr = pass_criteria.get("max_abs_rel_mass")
            if not (is_finite_number(thr) and float(thr) >= 0.0 and max_abs_rel_mass <= float(thr)):
                pass_under_user_criteria = False
                failed_criteria.append("max_abs_rel_mass")
        if "rmse_rel_mass" in pass_criteria:
            thr = pass_criteria.get("rmse_rel_mass")
            if not (is_finite_number(thr) and float(thr) >= 0.0 and rmse_rel_mass <= float(thr)):
                pass_under_user_criteria = False
                failed_criteria.append("rmse_rel_mass")
        if "rmse_log_m2" in pass_criteria:
            thr = pass_criteria.get("rmse_log_m2")
            if not (is_finite_number(thr) and float(thr) >= 0.0 and rmse_log_m2 <= float(thr)):
                pass_under_user_criteria = False
                failed_criteria.append("rmse_log_m2")

    if pass_claim_available and pass_under_user_criteria is True:
        status = "PASS_NONSTRICT_SM_HOST_MATCHING_UNDER_USER_THRESHOLDS"
    elif pass_claim_available and pass_under_user_criteria is False:
        status = "FAIL_NONSTRICT_SM_HOST_MATCHING_UNDER_USER_THRESHOLDS"
    else:
        status = "COMPUTED_NONSTRICT_SM_HOST_MATCHING_METRICS_NO_PASS_CLAIM"

    artifact = {
        "stage": "P702",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "nonstrict_external_host_matching_only",
        "goal": "compute_non_strict_host_matching_metrics_between_P696_channel_m2_proxies_and_an_explicit_external_standard_model_mass_target_dataset__no_false_pass",
        "inputs": {
            "P696_summary": str(IN_P696_SUMMARY.relative_to(REPO)),
            "sm_mass_targets": str(IN_SM_TARGETS.relative_to(REPO)),
        },
        "p696_mixing_audit": {
            "H_in_selector_aligned_basis_offdiag_to_diag_ratio": p696.get("H_in_selector_aligned_basis_offdiag_to_diag_ratio"),
            "offblock_max_fro": p696.get("offblock_max_fro"),
        },
        "method": {
            "assignment": "linear_sum_assignment on squared differences of centered log(mass^2)",
            "scale": "global scale set by mean log(mass^2/ proxy_m2) after assignment",
        },
        "units": units,
        "scale": {
            "scale_mass2_per_proxy_unit": scale_m2,
            "scale_mass_per_sqrt_proxy_unit": scale_m,
        },
        "metrics": {
            "rmse_log_m2": rmse_log_m2,
            "max_abs_log_m2": max_abs_log_m2,
            "rmse_rel_mass": rmse_rel_mass,
            "max_abs_rel_mass": max_abs_rel_mass,
        },
        "assignment": assignment_sorted_by_target,
        "pass_policy": {
            "pass_criteria": pass_criteria,
            "pass_claim_available": pass_claim_available,
            "pass_under_user_criteria": pass_under_user_criteria,
            "failed_criteria": failed_criteria,
        },
        "hard_limits": [
            "Non-strict: depends on an explicit external target dataset and identification policy.",
            "Does not claim any Standard Model identification or host matching discharge in strict scope.",
            "Does not upgrade any proxy spectrum into physical masses without an explicit model-level map.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P702",
        "status": status,
        "units": units,
        "scale_mass2_per_proxy_unit": scale_m2,
        "rmse_log_m2": rmse_log_m2,
        "max_abs_rel_mass": max_abs_rel_mass,
        "pass_under_user_criteria": pass_under_user_criteria,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

