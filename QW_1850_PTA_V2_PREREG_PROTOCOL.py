#!/usr/bin/env python3
"""
QW-1850: PTA V2 prereg protocol (versioned criterion reparam).

Goal:
- define a scientifically coherent replacement for PTA V1 probability criterion,
  based on pair-level unit of analysis,
- freeze rules for future EXTERNAL confirmatory execution.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1850_pta_v2_prereg_protocol.json"
OUT_MD = ROOT / "RAPORT_QW1850_PTA_V2_PREREG_PROTOCOL.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1848 = load("report_qw1848_pta_unit_of_analysis_audit.json")
    d1849 = load("report_qw1849_pta_strict_path_selection_gate.json")
    d1839 = load("report_qw1839_joint_confirmatory_prereg_protocol.json")

    if d1849.get("decision", {}).get("best_path") != "PATH_B_VERSIONED_CRITERION_REPARAM_WITH_EXTERNAL_CONFIRMATORY":
        raise RuntimeError("QW-1849 does not select PATH_B; PTA V2 prereg freeze blocked.")

    pta_ref = d1848.get("summary", {})
    pair_ref = pta_ref.get("pair_level", {})
    inf_ref = pta_ref.get("inference", {})

    gw_protocol = d1839.get("protocol", {}).get("gw_protocol", {})

    protocol = {
        "stage": "JOINT_CONFIRMATORY_V2_PTA_REPARAM",
        "frozen_utc": datetime.now(timezone.utc).isoformat(),
        "eligibility": {
            "requires_qw1849_best_path": "PATH_B_VERSIONED_CRITERION_REPARAM_WITH_EXTERNAL_CONFIRMATORY",
            "requires_external_confirmatory_dataset": True,
            "forbid_reusing_qw1848_design_dataset_for_decision": True,
        },
        "pta_v2_protocol": {
            "source_script": "QW_1848_PTA_UNIT_OF_ANALYSIS_AUDIT.py",
            "unit_of_analysis": "pair",
            "split_rule": "repeated_stratified_oos_splits",
            "n_replications_min": 60,
            "metrics": {
                "primary": "mean_pair_mean_gain",
                "primary_significance": "bootstrap_lower95_mean_pair_mean_gain",
                "secondary": "prob_pair_mean_gain_positive",
                "secondary_significance": "one_sided_lower95_prob_pair_mean_gain_positive",
            },
            "thresholds": {
                "mean_pair_mean_gain_min": 0.040,
                "bootstrap_lower95_mean_pair_mean_gain_min": 0.000,
                "prob_pair_mean_gain_positive_min": 0.667,
                "one_sided_lower95_prob_pair_mean_gain_positive_min": 0.600,
            },
            "reference_values_qw1848": {
                "mean_pair_mean_gain": pair_ref.get("mean_pair_mean_gain"),
                "prob_pair_mean_gain_positive": pair_ref.get("prob_pair_mean_gain_positive"),
                "one_sided_lower95_prob_pair_mean_gain_positive": inf_ref.get("lower95_one_sided_for_prob_pair_mean_positive"),
            },
        },
        "gw_protocol": gw_protocol,
        "joint_gate": {
            "logic": "AND",
            "require_all_pta_v2_thresholds": True,
            "require_all_gw_thresholds": True,
            "external_dataset_only": True,
        },
        "anti_leakage": {
            "forbidden": [
                "threshold_tuning_after_external_run",
                "redefining_units_after_results",
                "using_design_dataset_as_confirmatory_dataset",
                "foldwise_information_leakage",
            ],
            "allowed": [
                "single external confirmatory run under frozen protocol",
                "separate exploratory branch with explicit version bump",
            ],
        },
        "status": {
            "is_preregistered": True,
            "is_confirmed": False,
            "confirmatory_requires_new_external_data": True,
        },
    }

    canonical = json.dumps(protocol, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    protocol_hash = hashlib.sha256(canonical.encode("utf-8")).hexdigest()

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": protocol,
        "protocol_sha256": protocol_hash,
        "verdict": "PTA_V2_PREREG_FROZEN_EXTERNAL_CONFIRMATORY_PENDING",
    }

    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    t = protocol["pta_v2_protocol"]["thresholds"]
    ref = protocol["pta_v2_protocol"]["reference_values_qw1848"]

    lines = [
        "# RAPORT QW-1850: PTA V2 PREREG PROTOCOL",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Verdict: **{output['verdict']}**",
        f"- Protocol SHA256: `{protocol_hash}`",
        "",
        "## PTA V2 Thresholds (Pair-Level)",
        f"- mean_pair_mean_gain >= {t['mean_pair_mean_gain_min']}",
        f"- bootstrap_lower95_mean_pair_mean_gain >= {t['bootstrap_lower95_mean_pair_mean_gain_min']}",
        f"- prob_pair_mean_gain_positive >= {t['prob_pair_mean_gain_positive_min']}",
        f"- one_sided_lower95_prob_pair_mean_gain_positive >= {t['one_sided_lower95_prob_pair_mean_gain_positive_min']}",
        "",
        "## Reference (Design Dataset, Non-Confirmatory)",
        f"- mean_pair_mean_gain: {ref.get('mean_pair_mean_gain')}",
        f"- prob_pair_mean_gain_positive: {ref.get('prob_pair_mean_gain_positive')}",
        f"- one_sided_lower95_prob_pair_mean_gain_positive: {ref.get('one_sided_lower95_prob_pair_mean_gain_positive')}",
        "",
        "## Constraints",
        "- External confirmatory dataset required.",
        "- Reusing design dataset for final claim is forbidden.",
        "- GW protocol inherited from QW-1839 (unchanged).",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1850] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1850] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
