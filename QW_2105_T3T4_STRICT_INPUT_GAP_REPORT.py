#!/usr/bin/env python3
"""
QW-2105: T3/T4 strict input gap report.

Purpose:
- provide a deterministic, machine-readable summary of what exact external
  input gaps still block strict closure for:
  1) h0/lambda (QW-2099 -> QW-2102 -> QW-2090)
  2) G_newton (QW-2101 -> QW-2103 -> QW-2092)
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
R2099 = ROOT / "report_qw2099_hz_external_decoupling_autocollector.json"
R2102 = ROOT / "report_qw2102_hz_decoupling_identifiability_gate.json"
R2101 = ROOT / "report_qw2101_gnewton_bridge_external_autocollector.json"
R2103 = ROOT / "report_qw2103_gnewton_dimensionless_provenance_gate.json"
R2104 = ROOT / "report_qw2104_t3t4_strict_preflight_gate.json"

OUT_JSON = ROOT / "report_qw2105_t3t4_strict_input_gap_report.json"
OUT_MD = ROOT / "RAPORT_QW2105_T3T4_STRICT_INPUT_GAP_REPORT.md"


def load(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def append_if_missing(lst: List[str], item: str) -> None:
    if item not in lst:
        lst.append(item)


def main() -> None:
    d2099 = load(R2099)
    d2102 = load(R2102)
    d2101 = load(R2101)
    d2103 = load(R2103)
    d2104 = load(R2104)

    # --- H(z) gap analysis ---
    m2099 = d2099.get("strict_identifiability_metrics", {})
    t2099 = d2099.get("strict_thresholds", {})
    f2099 = d2099.get("strict_identifiability_flags", {})
    f2102 = d2102.get("flags", {})

    hz_gaps: List[str] = []
    if not bool(f2099.get("n_nodes_ge_5", False)):
        cur = float(m2099.get("n_nodes", 0.0))
        req = int(t2099.get("min_nodes", 5))
        append_if_missing(hz_gaps, f"add_nodes: at least {max(0, req - int(cur))} more valid H(z) nodes")
    if not bool(f2099.get("z_span_ge_0p8", False)):
        z_min = m2099.get("z_min")
        if isinstance(z_min, (int, float)):
            target = float(z_min) + float(t2099.get("min_z_span", 0.8))
            append_if_missing(
                hz_gaps,
                f"extend_redshift_span: include at least one validated node at z >= {target:.3f}",
            )
        else:
            append_if_missing(hz_gaps, "extend_redshift_span: z_span must be >= 0.8")
    if not bool(f2099.get("e_span_ge_1p0", False)):
        append_if_missing(
            hz_gaps,
            "increase_E_span: add validated nodes at sufficiently separated z so E(z) span >= 1.0",
        )
    if not bool(f2099.get("design_condition_lt_8", False)):
        cond = m2099.get("design_condition_number")
        if isinstance(cond, (int, float)) and math.isfinite(float(cond)):
            append_if_missing(
                hz_gaps,
                f"improve_design_condition: reduce cond([E,1]) below 8 (current {float(cond):.3f}) via wider z coverage",
            )
        else:
            append_if_missing(hz_gaps, "improve_design_condition: cond([E,1]) must be finite and < 8")
    if not bool(f2102.get("source_metadata_complete", False)):
        append_if_missing(hz_gaps, "complete_source_metadata: source/citation/reference_url + integrity metadata")
    if not bool(f2102.get("provenance_anchor_free", False)):
        append_if_missing(hz_gaps, "set_provenance_anchor_free_true")

    hz_ready = len(hz_gaps) == 0

    # --- G_newton gap analysis ---
    s2101 = d2101.get("bridge", {})
    f2101 = d2101.get("flags", {})
    f2103 = d2103.get("flags", {})
    g_gaps: List[str] = []

    if not bool(f2101.get("dimensionless_directly_provided", False)):
        append_if_missing(
            g_gaps,
            "provide_direct_dimensionless_bridge: g_dimensionless_mu_ref must come from external observable, not backsolve",
        )
    if not bool(f2103.get("bridge_origin_external_dimensionless", False)):
        append_if_missing(g_gaps, "set_bridge_observable_origin_external_dimensionless")
    if not bool(f2103.get("provenance_anchor_free", False)):
        append_if_missing(g_gaps, "set_provenance_anchor_free_true")
    if not bool(f2103.get("not_seeded_from_registry", False)):
        append_if_missing(g_gaps, "set_seeded_from_registry_false")
    if not bool(f2103.get("g_si_not_primary_input", False)):
        append_if_missing(g_gaps, "set_g_si_input_optional_null")
    if not bool(f2103.get("source_metadata_complete", False)):
        append_if_missing(g_gaps, "complete_source_metadata: source/citation/reference_url + integrity metadata")

    g_ready = len(g_gaps) == 0

    verdict = (
        "T3T4_STRICT_INPUT_GAPS_CLOSED_READY_FOR_STRICT_RERUN"
        if hz_ready and g_ready
        else "T3T4_STRICT_INPUT_GAPS_PRESENT"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "qw2099": R2099.name,
            "qw2102": R2102.name,
            "qw2101": R2101.name,
            "qw2103": R2103.name,
            "qw2104": R2104.name,
        },
        "hz_path": {
            "strict_ready": hz_ready,
            "gaps": hz_gaps,
            "current_metrics": {
                "n_nodes": m2099.get("n_nodes"),
                "z_span": m2099.get("z_span"),
                "e_span": m2099.get("e_span"),
                "design_condition_number": m2099.get("design_condition_number"),
            },
            "thresholds": {
                "min_nodes": t2099.get("min_nodes"),
                "min_z_span": t2099.get("min_z_span"),
                "min_e_span": t2099.get("min_e_span"),
                "max_design_condition_number": t2099.get("max_design_condition_number"),
            },
        },
        "gnewton_path": {
            "strict_ready": g_ready,
            "gaps": g_gaps,
            "current_bridge": {
                "bridge_observable_origin": s2101.get("bridge_observable_origin"),
                "strict_provenance_ready": s2101.get("strict_provenance_ready"),
            },
        },
        "t3t4_meta_gate": {
            "qw2104_verdict": d2104.get("verdict"),
            "qw2104_pass_count": d2104.get("pass_count"),
            "qw2104_total_flags": d2104.get("total_flags"),
        },
        "verdict": verdict,
        "required_next_step": (
            "COLLECT_STRICT_READY_HZ_AND_GNEWTON_EXTERNAL_INPUTS_THEN_RERUN_QW2099_QW2101_QW2102_QW2103_QW2090_QW2092_QW2104_QW2094"
            if verdict == "T3T4_STRICT_INPUT_GAPS_PRESENT"
            else "RERUN_T3T4_STRICT_CHAIN_AND_PROMOTE_CLOSURE"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2105: T3T4 STRICT INPUT GAP REPORT",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## H(z) Path",
        f"- strict_ready: `{hz_ready}`",
        f"- current n_nodes: `{m2099.get('n_nodes')}`",
        f"- current z_span: `{m2099.get('z_span')}`",
        f"- current e_span: `{m2099.get('e_span')}`",
        f"- current cond([E,1]): `{m2099.get('design_condition_number')}`",
        "- gaps:",
    ]
    if hz_gaps:
        lines.extend([f"  - {g}" for g in hz_gaps])
    else:
        lines.append("  - none")

    lines.extend(
        [
            "",
            "## G_newton Path",
            f"- strict_ready: `{g_ready}`",
            f"- bridge_observable_origin: `{s2101.get('bridge_observable_origin')}`",
            "- gaps:",
        ]
    )
    if g_gaps:
        lines.extend([f"  - {g}" for g in g_gaps])
    else:
        lines.append("  - none")

    lines.extend(
        [
            "",
            "## Meta",
            f"- QW-2104 verdict: `{d2104.get('verdict')}` ({d2104.get('pass_count')}/{d2104.get('total_flags')})",
            "",
            "## Artifact",
            f"- JSON: `{OUT_JSON.name}`",
            f"- required_next_step: `{out['required_next_step']}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2105] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2105] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2105] verdict={verdict} hz_ready={hz_ready} g_ready={g_ready}")


if __name__ == "__main__":
    main()
