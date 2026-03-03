#!/usr/bin/env python3
"""
QW-1899: External detector protocol design for amplitude-lite branch.

Creates a preregistered, hash-locked protocol file with fixed PASS/FAIL
criteria derived from QW-1897/1898 performance envelope.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
PROTO_DIR = ROOT / "external_confirmatory_v2"
PROTO_JSON = PROTO_DIR / "protocol_qw1899_external_detector.json"
OUT_JSON = ROOT / "report_qw1899_external_detector_protocol_design.json"
OUT_MD = ROOT / "RAPORT_QW1899_EXTERNAL_DETECTOR_PROTOCOL_DESIGN.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def sha256_bytes(b: bytes) -> str:
    h = hashlib.sha256()
    h.update(b)
    return h.hexdigest()


def main() -> None:
    d1897 = read_json("report_qw1897_admissible_amplitude_lite_multisplit.json")
    d1898 = read_json("report_qw1898_empirical_bridge_precondition_gate.json")

    s = d1897.get("summary", {})

    # Conservative prereg thresholds from achieved multisplit medians.
    thr = {
        "success_rate_min": max(0.60, round(0.80 * float(s.get("success_rate", 0.72)), 3)),
        "rmse_gain_median_min": max(0.010, round(0.50 * float(s.get("rmse_gain_median", 0.0265)), 4)),
        "canon_gain_median_min": max(0.025, round(0.50 * float(s.get("canon_gain_median", 0.0677)), 4)),
        "nonboundary_gain_median_min": max(0.250, round(0.50 * float(s.get("nonboundary_gain_median", 0.5000)), 4)),
    }

    protocol = {
        "protocol_id": "QW-1899-EXTERNAL-DETECTOR-PROTOCOL",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "model_lock": {
            "model_name": "admissible_amplitude_lite",
            "source_reports": [
                "report_qw1896_admissible_amplitude_lite_gate.json",
                "report_qw1897_admissible_amplitude_lite_multisplit.json",
            ],
            "complexity": {
                "n_obs_per_profile": 12,
                "k_model": 10,
                "residual_dof": 2,
            },
            "fixed_hyperparams": {
                "lambda_c": 0.2,
                "lambda_b": float(d1897.get("protocol", {}).get("lambda_b", 0.35)),
                "lambda_bias": float(d1897.get("protocol", {}).get("lambda_bias", 0.05)),
                "constraints": {
                    "omega": "pi/4 +/- 0.12",
                    "phi": "{+/- pi/6} +/- 0.22",
                    "beta": "[0.005, 0.05]",
                },
            },
        },
        "dataset_requirements": {
            "externality": "must be independent from all design/tuning datasets",
            "providers": ["LIGO", "Virgo", "KAGRA", "PTA collaboration"],
            "freeze_required": True,
            "required_files": ["manifest.json", "pta_v2_pairs.csv", "gw_windows.csv"],
            "manifest_fields": [
                "dataset.dataset_id",
                "dataset.provider",
                "dataset.externality_statement",
                "files[].sha256",
            ],
        },
        "primary_metrics": {
            "success_rate": "fraction of external splits passing all three median-gain constraints",
            "rmse_gain_median": "median(ctrl_rmse - model_rmse)",
            "canon_gain_median": "median(model_canon - ctrl_canon)",
            "nonboundary_gain_median": "median(model_nonboundary - ctrl_nonboundary)",
        },
        "pass_thresholds": thr,
        "hard_fail_conditions": [
            "residual_dof < 2",
            "any post-hoc threshold change",
            "dataset hash mismatch vs manifest",
            "reuse of design/training dataset",
        ],
        "decision_rule": "PASS only if all thresholds met; else FAIL",
        "anti_p_hacking": {
            "single_protocol_hash": True,
            "no_hyperparam_retuning_on_external": True,
            "no_metric_substitution": True,
        },
        "bridge_context": {
            "precondition_report": "report_qw1898_empirical_bridge_precondition_gate.json",
            "precondition_readiness": d1898.get("readiness", "UNKNOWN"),
        },
    }

    PROTO_DIR.mkdir(parents=True, exist_ok=True)
    proto_bytes = json.dumps(protocol, ensure_ascii=False, indent=2).encode("utf-8")
    proto_hash = sha256_bytes(proto_bytes)
    PROTO_JSON.write_bytes(proto_bytes)

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_path": str(PROTO_JSON.relative_to(ROOT)),
        "protocol_sha256": proto_hash,
        "pass_thresholds": thr,
        "source_summary": s,
        "precondition_readiness": d1898.get("readiness", "UNKNOWN"),
        "verdict": "EXTERNAL_DETECTOR_PROTOCOL_FROZEN",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1899: EXTERNAL DETECTOR PROTOCOL DESIGN",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- Protocol: `{out['protocol_path']}`",
        f"- Protocol SHA256: `{proto_hash}`",
        "",
        "## Locked Thresholds",
        f"- success_rate_min: {thr['success_rate_min']}",
        f"- rmse_gain_median_min: {thr['rmse_gain_median_min']}",
        f"- canon_gain_median_min: {thr['canon_gain_median_min']}",
        f"- nonboundary_gain_median_min: {thr['nonboundary_gain_median_min']}",
        "",
        "## Next",
        "- Prepare external dataset package (`manifest.json`, `pta_v2_pairs.csv`, `gw_windows.csv`) and run QW-1900 frozen external execution.",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1899] Saved protocol: {PROTO_JSON.name}")
    print(f"[QW-1899] Saved JSON:     {OUT_JSON.name}")
    print(f"[QW-1899] Saved MD:       {OUT_MD.name}")


if __name__ == "__main__":
    main()
