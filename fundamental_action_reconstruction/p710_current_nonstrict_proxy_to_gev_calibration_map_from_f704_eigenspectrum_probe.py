#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

AS_OF = "2026-03-17"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F704_OBJECT = GENERATED / "mass_observable_diagonal_local_strict_derived_v1.json"
IN_TARGETS = ROOT / "external_data" / "sm_mass_targets_v1.json"
IN_POLICY = ROOT / "external_data" / "proxy_to_gev_calibration_policy_v1.json"

OUT = GENERATED / "p710_current_nonstrict_proxy_to_gev_calibration_map_from_f704_eigenspectrum_probe.json"
OUT_SUMMARY = GENERATED / "p710_current_nonstrict_proxy_to_gev_calibration_map_from_f704_eigenspectrum_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_finite_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F704_OBJECT, IN_TARGETS, IN_POLICY]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P710",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f704 = load_json(IN_F704_OBJECT)
    out = (f704.get("outputs") or {}) if isinstance(f704, dict) else {}
    proxy_m2 = out.get("eigenvalues_m2_proxy_ascending")
    if not (
        isinstance(proxy_m2, list)
        and len(proxy_m2) == 12
        and all(is_finite_number(v) and float(v) > 0.0 for v in proxy_m2)
    ):
        artifact = {
            "stage": "P710",
            "status": "INVALID_F704_OBJECT_SHAPE",
            "as_of": AS_OF,
            "error": "F704 outputs must export eigenvalues_m2_proxy_ascending as length-12 list[positive finite number].",
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    targets = load_json(IN_TARGETS)
    units = targets.get("units")
    targets_list = targets.get("targets")
    if not (isinstance(units, str) and units):
        artifact = {
            "stage": "P710",
            "status": "INVALID_TARGETS_UNITS",
            "as_of": AS_OF,
            "error": "sm_mass_targets_v1.json must define nonempty string units (e.g. 'GeV').",
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    if not (isinstance(targets_list, list) and all(isinstance(t, dict) for t in targets_list) and len(targets_list) == 12):
        artifact = {
            "stage": "P710",
            "status": "INVALID_TARGETS_LIST",
            "as_of": AS_OF,
            "error": "sm_mass_targets_v1.json must define targets as list[object] of length 12 for this v1 map.",
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    target_ids: list[str] = []
    target_m: list[float] = []
    invalid: list[str] = []
    for i, t in enumerate(targets_list):
        tid = t.get("id")
        m = t.get("mass")
        if not isinstance(tid, str) or not tid:
            invalid.append(f"targets[{i}].id")
            continue
        if not (is_finite_number(m) and float(m) > 0.0):
            invalid.append(f"targets[{i}].mass:{tid}")
            continue
        target_ids.append(tid)
        target_m.append(float(m))

    if invalid:
        artifact = {
            "stage": "P710",
            "status": "INVALID_TARGET_FIELDS",
            "as_of": AS_OF,
            "invalid_fields": invalid,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    policy = load_json(IN_POLICY)
    if not isinstance(policy, dict):
        artifact = {
            "stage": "P710",
            "status": "INVALID_POLICY_SHAPE",
            "as_of": AS_OF,
            "error": "proxy_to_gev_calibration_policy_v1.json must be a JSON object.",
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    expected_dataset_id = policy.get("expected_dataset_id")
    if isinstance(expected_dataset_id, str) and expected_dataset_id:
        if targets.get("dataset_id") != expected_dataset_id:
            artifact = {
                "stage": "P710",
                "status": "NOT_COMPUTABLE_POLICY_DATASET_ID_MISMATCH",
                "as_of": AS_OF,
                "dataset_id": targets.get("dataset_id"),
                "expected_dataset_id": expected_dataset_id,
                "no_false_pass": True,
            }
            OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
            OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
            print(OUT_SUMMARY)
            return

    dof = policy.get("degrees_of_freedom") or {}
    if not (
        isinstance(dof, dict)
        and dof.get("global_scale_only") is True
        and dof.get("sector_scales_allowed") is False
        and dof.get("kinetic_metric_refit_allowed") is False
        and dof.get("nonlinear_mass_map_allowed") is False
    ):
        artifact = {
            "stage": "P710",
            "status": "NOT_COMPUTABLE_UNSUPPORTED_POLICY_DEGREES_OF_FREEDOM",
            "as_of": AS_OF,
            "degrees_of_freedom": dof,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    # Compute the single global scale using only geometric means of m^2.
    proxy_m2_arr = np.array([float(v) for v in proxy_m2], dtype=float)
    target_m2_arr = np.array([float(m) * float(m) for m in target_m], dtype=float)
    log_scale = float(np.mean(np.log(target_m2_arr)) - np.mean(np.log(proxy_m2_arr)))
    scale_m2_per_proxy_unit = float(math.exp(log_scale))
    scale_m_per_sqrt_proxy_unit = float(math.sqrt(scale_m2_per_proxy_unit))

    calibrated_m2 = (scale_m2_per_proxy_unit * proxy_m2_arr).tolist()
    calibrated_m = (np.sqrt(scale_m2_per_proxy_unit * proxy_m2_arr)).tolist()

    # Sanity: geometric mean m^2 should match by construction (within float).
    gm_target_m2 = float(math.exp(float(np.mean(np.log(target_m2_arr)))))
    gm_proxy_m2 = float(math.exp(float(np.mean(np.log(proxy_m2_arr)))))
    gm_calibrated_m2 = float(math.exp(float(np.mean(np.log(np.array(calibrated_m2, dtype=float))))))

    status = "PASS_EXPORTED_NONSTRICT_PROXY_TO_GEV_CALIBRATION_MAP_FROM_F704_EIGENSPECTRUM"
    artifact = {
        "stage": "P710",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "nonstrict_units_calibration_map_only",
        "inputs": {
            "proxy_object": str(IN_F704_OBJECT.relative_to(REPO)),
            "proxy_field": "outputs.eigenvalues_m2_proxy_ascending",
            "targets_dataset": str(IN_TARGETS.relative_to(REPO)),
            "policy": str(IN_POLICY.relative_to(REPO)),
        },
        "dataset": {
            "dataset_id": targets.get("dataset_id"),
            "units": units,
            "target_ids": target_ids,
        },
        "policy": {
            "policy_id": policy.get("policy_id"),
            "degrees_of_freedom": dof,
            "pass_criteria": policy.get("pass_criteria"),
        },
        "calibration_map": {
            "units": units,
            "scale_m2_per_proxy_unit": scale_m2_per_proxy_unit,
            "scale_m_per_sqrt_proxy_unit": scale_m_per_sqrt_proxy_unit,
            "formula": "m2_physical := scale_m2_per_proxy_unit * m2_proxy; m_physical := sqrt(m2_physical)",
            "scale_construction": "geometric_mean_match_on_m2 (assignment-free)",
            "notes": [
                "This is a non-strict calibration interface only.",
                "This map does not identify which proxy modes correspond to which particles.",
                "Use P704/P702 for diagnostic-only matching metrics under a separate non-strict host-matching policy.",
            ],
        },
        "sanity": {
            "geometric_mean_target_m2": gm_target_m2,
            "geometric_mean_proxy_m2": gm_proxy_m2,
            "geometric_mean_calibrated_proxy_m2": gm_calibrated_m2,
            "abs_diff_geometric_mean_m2": abs(gm_target_m2 - gm_calibrated_m2),
        },
        "results": {
            "proxy_m2_proxy_ascending": [float(x) for x in proxy_m2_arr.tolist()],
            "calibrated_m2_physical_ascending": [float(x) for x in calibrated_m2],
            "calibrated_m_physical_ascending": [float(x) for x in calibrated_m],
        },
        "hard_limits": [
            "Non-strict policy layer only; does not claim a strict physical-unit map.",
            "Does not claim Standard Model identification or 'match'.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim a directed/sign-sensitive physical orientation datum in strict core.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P710",
        "status": status,
        "dataset_id": targets.get("dataset_id"),
        "policy_id": policy.get("policy_id"),
        "units": units,
        "scale_m2_per_proxy_unit": scale_m2_per_proxy_unit,
        "scale_m_per_sqrt_proxy_unit": scale_m_per_sqrt_proxy_unit,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

