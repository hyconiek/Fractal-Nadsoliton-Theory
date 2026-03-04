#!/usr/bin/env python3
"""
QW-2076: Empirical prediction preregistration package.

Purpose:
- freeze falsifiable prospective predictions before future data arrives,
- define explicit pass/fail rules to avoid post-hoc reinterpretation.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2076_empirical_prediction_preregistration.json"
OUT_MD = ROOT / "RAPORT_QW2076_EMPIRICAL_PREDICTION_PREREGISTRATION.md"
OBS_TEMPLATE = ROOT / "empirical_observations_input_qw2077.template.json"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    r2049 = load_json("report_qw2049_spectral_micro_stagec_intersection_gate.json")
    r2063 = load_json("report_qw2063_derivational_reconstruction_shared_flavor_basis.json")
    r2073 = load_json("report_qw2073_radiative_channels_closure_upgrade.json")
    r2075 = load_json("report_qw2075_strict_cp_phase_derivation_gate.json")

    pmns = r2075["pmns_phase"]
    gw_metrics = r2063["metrics"]["gw"]
    gw_weights = r2063["derived_gw_weights"]
    gw_thresholds = {
        "auc_h1l1_vs_ctrl_min": float(r2063["thresholds"]["gw_auc_min"]),
        "adv_shared_minus_ctrl_q90_min": float(r2063["thresholds"]["gw_adv_min"]),
        "sep_median_h1l1_minus_ctrl_min": float(r2063["thresholds"]["gw_sep_min"]),
        "control_median_gap_max": float(r2063["thresholds"]["gw_control_gap_max"]),
    }

    z_grid = [0.0, 0.5, 1.0, 2.0, 3.0]
    w_eff_grid = [float(x) for x in r2073["diagnostics"]["w_eff_grid"]]

    prereg_files: List[str] = [
        "report_qw2049_spectral_micro_stagec_intersection_gate.json",
        "report_qw2063_derivational_reconstruction_shared_flavor_basis.json",
        "report_qw2073_radiative_channels_closure_upgrade.json",
        "report_qw2075_strict_cp_phase_derivation_gate.json",
    ]
    frozen_dependencies = {
        name: sha256_file(ROOT / name) for name in prereg_files
    }

    predictions = [
        {
            "id": "PRED-2076-PMNS-CP",
            "domain": "neutrino_oscillation",
            "observable": "sin(delta_cp_pmns)",
            "novelty_status": "not_yet_precisely_resolved_experimentally",
            "predicted_point": float(pmns["sin_delta"]),
            "predicted_phase_branches_rad": [
                float(pmns["delta_primary_rad"]),
                float(pmns["delta_secondary_rad"]),
            ],
            "pre_registered_band": {
                "sin_delta_min": -0.10,
                "sin_delta_max": 0.10,
            },
            "falsification_rule": (
                "FALSIFIED if future independent 95% CI for sin(delta_cp_pmns) is entirely outside [-0.10, 0.10]."
            ),
            "support_rule": (
                "SUPPORTED if central estimate is inside [-0.10, 0.10] and CI95 overlaps the same band."
            ),
            "required_future_input": [
                "sin_delta_central",
                "sin_delta_ci95_low",
                "sin_delta_ci95_high",
                "source",
            ],
        },
        {
            "id": "PRED-2076-COSMO-WEFF",
            "domain": "cosmology",
            "observable": "w_eff(z) at fixed nodes",
            "novelty_status": "prospective_shape_test",
            "predicted_nodes": [
                {"z": float(z), "w_eff_pred": float(w)}
                for z, w in zip(z_grid, w_eff_grid)
            ],
            "absolute_corridor": 0.02,
            "falsification_rule": (
                "FALSIFIED if at any node |w_obs - w_pred| > max(0.02, 3*sigma_total) "
                "in an independent joint-fit dataset."
            ),
            "support_rule": (
                "SUPPORTED if all nodes satisfy |w_obs - w_pred| <= max(0.02, 3*sigma_total)."
            ),
            "required_future_input": [
                "z",
                "w_eff_central",
                "sigma_total",
                "source",
            ],
        },
        {
            "id": "PRED-2076-GW-HOLDOUT",
            "domain": "gravitational_waves",
            "observable": "future independent holdout metrics with locked feature pipeline",
            "novelty_status": "future_dataset_prediction",
            "locked_kernel": {k: float(v) for k, v in r2049["stagec_pool"]["selected_kernel"].items()},
            "locked_weights": {k: float(v) for k, v in gw_weights.items()},
            "reference_current_metrics": {k: float(v) for k, v in gw_metrics.items()},
            "hard_thresholds": gw_thresholds,
            "falsification_rule": (
                "FALSIFIED if future independent holdout does not satisfy all hard thresholds "
                "(AUC, ADV, SEP, CONTROL_GAP) under the locked pipeline."
            ),
            "support_rule": "SUPPORTED if all hard thresholds are satisfied on future independent holdout.",
            "required_future_input": [
                "auc_h1l1_vs_ctrl",
                "adv_shared_minus_ctrl_q90",
                "sep_median_h1l1_minus_ctrl",
                "control_median_gap",
                "n_windows",
                "source",
            ],
        },
    ]

    obs_template = {
        "pmns_cp": {
            "sin_delta_central": None,
            "sin_delta_ci95_low": None,
            "sin_delta_ci95_high": None,
            "source": "",
        },
        "cosmology_weff_nodes": [
            {"z": float(z), "w_eff_central": None, "sigma_total": None, "source": ""}
            for z in z_grid
        ],
        "gw_future_holdout": {
            "auc_h1l1_vs_ctrl": None,
            "adv_shared_minus_ctrl_q90": None,
            "sep_median_h1l1_minus_ctrl": None,
            "control_median_gap": None,
            "n_windows": None,
            "source": "",
        },
        "notes": "Fill with truly new independent observations, then run QW_2077.",
    }
    OBS_TEMPLATE.write_text(json.dumps(obs_template, ensure_ascii=False, indent=2), encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "verdict": "EMPIRICAL_PREDICTION_PREREGISTRATION_READY",
        "protocol_version": "1.0",
        "frozen_dependencies_sha256": frozen_dependencies,
        "predictions": predictions,
        "observation_input_template": OBS_TEMPLATE.name,
        "required_next_step": "COLLECT_TRULY_NEW_EXTERNAL_DATA_AND_RUN_QW2077_VALIDATOR",
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2076: EMPIRICAL PREDICTION PREREGISTRATION",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- Predictions preregistered: {len(predictions)}",
        "",
        "## Prediction IDs",
        "- PRED-2076-PMNS-CP",
        "- PRED-2076-COSMO-WEFF",
        "- PRED-2076-GW-HOLDOUT",
        "",
        "## Key PMNS Forecast",
        f"- sin(delta_cp_pmns) point: {pmns['sin_delta']:.9f}",
        "- pre-registered test band: [-0.10, 0.10]",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
        f"- observation template: `{OBS_TEMPLATE.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2076] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2076] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2076] Saved template: {OBS_TEMPLATE.name}")
    print(f"[QW-2076] verdict={out['verdict']} predictions={len(predictions)}")


if __name__ == "__main__":
    main()

