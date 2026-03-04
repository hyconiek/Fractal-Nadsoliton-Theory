#!/usr/bin/env python3
"""
QW-2111: strict external acquisition packet for T3/T4 closure.

Purpose:
- convert gap diagnostics into an operational acquisition checklist,
- keep the path falsifiable and deterministic (no hidden retune).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2111_t3t4_strict_external_acquisition_packet.json"
OUT_MD = ROOT / "RAPORT_QW2111_T3T4_STRICT_EXTERNAL_ACQUISITION_PACKET.md"


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _top_pairs_from_qw2107(rep: Dict, k: int = 5) -> List[List[float]]:
    top = rep.get("top_solutions", [])
    out: List[List[float]] = []
    if not isinstance(top, list):
        return out
    for row in top[:k]:
        if not isinstance(row, dict):
            continue
        nodes = row.get("added_nodes", [])
        if isinstance(nodes, list) and len(nodes) == 2:
            out.append([float(nodes[0]), float(nodes[1])])
    return out


def main() -> None:
    r2105 = load_json(ROOT / "report_qw2105_t3t4_strict_input_gap_report.json")
    r2106 = load_json(ROOT / "report_qw2106_strict_external_input_intake_gate.json")
    r2107 = load_json(ROOT / "report_qw2107_hz_strict_design_search.json")
    r2108 = load_json(ROOT / "report_qw2108_gnewton_dimensionless_acquisition_spec.json")
    r2109 = load_json(ROOT / "report_qw2109_strict_external_evidence_manifest_gate.json")

    hz_metrics = (r2105.get("hz_path") or {}).get("current_metrics", {})
    hz_thresholds = (r2105.get("hz_path") or {}).get("thresholds", {})
    hz_gaps = (r2105.get("hz_path") or {}).get("gaps", [])
    hz_pairs = _top_pairs_from_qw2107(r2107, k=10)

    g_gaps = (r2105.get("gnewton_path") or {}).get("gaps", [])
    g_spec = r2108.get("bridge_spec", {})
    g_contract = r2108.get("strict_contract", {})

    packet = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "qw2105": "report_qw2105_t3t4_strict_input_gap_report.json",
            "qw2106": "report_qw2106_strict_external_input_intake_gate.json",
            "qw2107": "report_qw2107_hz_strict_design_search.json",
            "qw2108": "report_qw2108_gnewton_dimensionless_acquisition_spec.json",
            "qw2109": "report_qw2109_strict_external_evidence_manifest_gate.json",
        },
        "status_snapshot": {
            "qw2105_verdict": r2105.get("verdict"),
            "qw2106_verdict": r2106.get("verdict"),
            "qw2106_pass_count": r2106.get("pass_count"),
            "qw2106_total_flags": r2106.get("total_flags"),
            "qw2109_verdict": r2109.get("verdict"),
            "qw2109_pass_count": r2109.get("pass_count"),
            "qw2109_total_flags": r2109.get("total_flags"),
        },
        "hz_acquisition_plan": {
            "objective": "Reach strict-ready H(z) identifiability for QW-2090/QW-2102.",
            "current_metrics": hz_metrics,
            "thresholds": hz_thresholds,
            "gaps": hz_gaps,
            "recommended_added_z_pairs_top10": hz_pairs,
            "required_fields_per_new_node": [
                "z",
                "h_km_s_mpc",
                "sigma_total",
                "citation",
                "reference_url",
                "source_version",
            ],
            "hard_rule": "New nodes must be external measurements; no synthetic/interpolated/manual-fit points.",
        },
        "gnewton_acquisition_plan": {
            "objective": "Replace SI-backsolved bridge with direct external dimensionless observable for QW-2092/QW-2103.",
            "current_gaps": g_gaps,
            "target_mu_ref_gev": g_spec.get("mu_ref_gev"),
            "g_dimensionless_target": g_spec.get("g_dimensionless_target"),
            "g_dimensionless_acceptance_range": g_spec.get("g_dimensionless_acceptance_range"),
            "strict_contract": g_contract,
            "required_payload_keys": [
                "mu_ref_gev",
                "g_dimensionless_mu_ref",
                "bridge_observable_origin",
                "notes",
            ],
            "hard_rule": "bridge_observable_origin must be external_dimensionless_observable (never backsolved_from_g_si).",
        },
        "evidence_sidecar_minimum": {
            "status": "PASS" if r2109.get("verdict") == "STRICT_EXTERNAL_EVIDENCE_MANIFEST_GATE_PASS" else "PENDING",
            "required_fields": [
                "source_dataset_name",
                "citation",
                "reference_url",
                "source_version",
                "acquired_utc",
                "artifact_sha256",
                "acquisition_protocol",
                "acquisition_command",
                "analyst",
            ],
            "consistency_requirements": [
                "artifact_sha256 must equal SHA256(payload)",
                "H(z): metadata record_count and columns_schema must match CSV",
                "G bridge: metadata json_keys must match payload keys",
            ],
        },
        "rerun_protocol_after_data_refresh": [
            "python3 QW_2112_HZ_STRICT_NODE_PACK_GATE.py",
            "python3 QW_2113_GNEWTON_DIRECT_DIMENSIONLESS_PACK_GATE.py",
            "python3 QW_2109_STRICT_EXTERNAL_EVIDENCE_MANIFEST_GATE.py",
            "python3 QW_2106_STRICT_EXTERNAL_INPUT_INTAKE_GATE.py",
            "python3 QW_2099_HZ_EXTERNAL_DECOUPLING_AUTOCOLLECTOR.py --nodes-csv external_hz_nodes_qw2099.csv --citation \"...\" --reference-url \"...\" --source-version \"...\" --require-strict-ready",
            "python3 QW_2101_GNEWTON_BRIDGE_EXTERNAL_AUTOCOLLECTOR.py --source-file external_gnewton_bridge_qw2101.json --citation \"...\" --reference-url \"...\" --source-version \"...\" --strict-dimensionless-only --omit-g-si-optional --require-strict-ready",
            "python3 QW_2102_HZ_DECOUPLING_IDENTIFIABILITY_GATE.py",
            "python3 QW_2103_GNEWTON_DIMENSIONLESS_PROVENANCE_GATE.py",
            "python3 QW_2090_H0_LAMBDA_DECOUPLING_GATE.py",
            "python3 QW_2092_GNEWTON_SI_BRIDGE_GATE.py",
            "python3 QW_2104_T3T4_STRICT_PREFLIGHT_GATE.py",
            "python3 QW_2105_T3T4_STRICT_INPUT_GAP_REPORT.py",
            "python3 QW_2094_STRICT_RIGOR_DEFECT_SWEEP.py",
        ],
        "verdict": "T3T4_STRICT_EXTERNAL_ACQUISITION_PACKET_READY",
        "required_next_step": "COLLECT_REAL_EXTERNAL_MEASUREMENTS_ACCORDING_TO_PACKET_AND_RERUN_PROTOCOL",
    }

    OUT_JSON.write_text(json.dumps(packet, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2111: T3/T4 STRICT EXTERNAL ACQUISITION PACKET",
        "",
        f"- Date UTC: {packet['generated_utc']}",
        f"- Verdict: **{packet['verdict']}**",
        "",
        "## Status Snapshot",
        f"- QW-2105: `{packet['status_snapshot']['qw2105_verdict']}`",
        f"- QW-2106: `{packet['status_snapshot']['qw2106_verdict']}` ({packet['status_snapshot']['qw2106_pass_count']}/{packet['status_snapshot']['qw2106_total_flags']})",
        f"- QW-2109: `{packet['status_snapshot']['qw2109_verdict']}` ({packet['status_snapshot']['qw2109_pass_count']}/{packet['status_snapshot']['qw2109_total_flags']})",
        "",
        "## H(z) Acquisition",
        f"- Current n_nodes: `{hz_metrics.get('n_nodes')}` (threshold `{hz_thresholds.get('min_nodes')}`)",
        f"- Current z_span: `{hz_metrics.get('z_span')}` (threshold `{hz_thresholds.get('min_z_span')}`)",
        f"- Current E_span: `{hz_metrics.get('e_span')}` (threshold `{hz_thresholds.get('min_e_span')}`)",
        f"- Current cond([E,1]): `{hz_metrics.get('design_condition_number')}` (threshold `< {hz_thresholds.get('max_design_condition_number')}`)",
        f"- Suggested added z pairs (top 10): `{hz_pairs}`",
        "",
        "## G Bridge Acquisition",
        f"- Target mu_ref_gev: `{g_spec.get('mu_ref_gev')}`",
        f"- Target g_dimensionless: `{g_spec.get('g_dimensionless_target')}`",
        f"- Accepted range: `{g_spec.get('g_dimensionless_acceptance_range')}`",
        f"- Hard origin requirement: `{g_contract.get('required_bridge_observable_origin')}`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2111] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2111] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2111] verdict={packet['verdict']}")


if __name__ == "__main__":
    main()
