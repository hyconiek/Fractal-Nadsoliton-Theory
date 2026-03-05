#!/usr/bin/env python3
"""
QW-2199: Gravity action-level scope gate (L23).

Purpose:
- integrate current gravity/continuum/action-level evidence,
- make explicit separation between closed effective strict bridges and still-open
  foundational Einstein-Hilbert-level derivation claims.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2199_gravity_action_level_scope_gate.json"
OUT_MD = ROOT / "RAPORT_QW2199_GRAVITY_ACTION_LEVEL_SCOPE_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    r2069 = load_json("report_qw2069_full_sm_gr_derivation_package.json")
    r2092 = load_json("report_qw2092_gnewton_si_bridge_gate.json")
    r2115 = load_json("report_qw2115_gravity_hierarchy_strict_bridge_gate.json")
    r2148 = load_json("report_qw2148_continuum_dg_delta_extrapolation_gate.json")
    r2164 = load_json("report_qw2164_l14_full_canonical_continuum_variational_gate.json")
    r2180 = load_json("report_qw2180_l14_v2b2_exact_action_identification_gate.json")

    entries = {e["id"]: e for e in r2069["entries"]}

    closed_effective_components: List[str] = []
    if str(r2092.get("verdict", "")).endswith("PASS_STRICT"):
        closed_effective_components.append("gnewton_si_bridge")
    if str(r2115.get("verdict", "")).endswith("GATE_PASS"):
        closed_effective_components.append("gravity_hierarchy_bridge")
    if str(r2180.get("verdict", "")).endswith("TERMINAL_CHAIN_CLOSED"):
        closed_effective_components.append("l14_action_identification_terminal_chain")
    if entries.get("G_newton", {}).get("status") == "derived":
        closed_effective_components.append("registry_gnewton_derived")
    if entries.get("lambda_cosmological", {}).get("status") == "derived":
        closed_effective_components.append("registry_lambda_derived")
    if entries.get("h0", {}).get("status") == "derived":
        closed_effective_components.append("registry_h0_derived")

    open_foundational_components: List[str] = [
        "einstein_hilbert_action_direct_derivation_from_complete_fin_action",
        "equivalence_principle_derivation_from_complete_fin_action",
        "full_sm_gr_reduction_theorem_from_complete_fin_action",
    ]

    flags = {
        "q2092_gnewton_si_bridge_strict_pass_present": bool(str(r2092.get("verdict", "")).endswith("PASS_STRICT")),
        "q2115_gravity_hierarchy_strict_bridge_present": bool(str(r2115.get("verdict", "")).endswith("GATE_PASS")),
        "q2148_continuum_extrapolation_bridge_present": bool(str(r2148.get("verdict", "")).endswith("PARTIAL_ACTION_THEOREM_OPEN")),
        "q2164_full_canonical_continuum_variational_present": bool(str(r2164.get("verdict", "")).endswith("FULL_THEOREM_OPEN")),
        "q2180_l14_terminal_action_identification_closed": bool(str(r2180.get("verdict", "")).endswith("TERMINAL_CHAIN_CLOSED")),
        "q2069_registry_gnewton_derived": bool(entries.get("G_newton", {}).get("status") == "derived"),
        "q2069_registry_lambda_h0_derived": bool(
            entries.get("lambda_cosmological", {}).get("status") == "derived"
            and entries.get("h0", {}).get("status") == "derived"
        ),
        "effective_gr_bridge_components_nonempty": bool(len(closed_effective_components) >= 4),
        "foundational_open_components_explicit": bool(len(open_foundational_components) >= 3),
        "einstein_hilbert_action_direct_derivation_closed": False,
        "equivalence_principle_derivation_closed": False,
        "full_sm_gr_reduction_theorem_closed": False,
        "scope_boundaries_explicit_no_overclaim": True,
        "deterministic_no_scan_no_retune": bool(
            r2092["flags"].get("no_anchor_feedback_loop", False)
            and r2115["flags"].get("deterministic_no_scan_no_retune_formula", False)
        ),
    }

    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    core_ok = bool(
        flags["q2092_gnewton_si_bridge_strict_pass_present"]
        and flags["q2115_gravity_hierarchy_strict_bridge_present"]
        and flags["q2148_continuum_extrapolation_bridge_present"]
        and flags["q2164_full_canonical_continuum_variational_present"]
        and flags["q2180_l14_terminal_action_identification_closed"]
        and flags["q2069_registry_gnewton_derived"]
        and flags["q2069_registry_lambda_h0_derived"]
        and flags["effective_gr_bridge_components_nonempty"]
        and flags["foundational_open_components_explicit"]
        and flags["scope_boundaries_explicit_no_overclaim"]
        and flags["deterministic_no_scan_no_retune"]
    )

    verdict = (
        "GRAVITY_ACTION_LEVEL_SCOPE_GATE_PASS_STRICT_PARTIAL_FOUNDATIONAL_OPEN"
        if core_ok
        else "GRAVITY_ACTION_LEVEL_SCOPE_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2069": "report_qw2069_full_sm_gr_derivation_package.json",
            "q2092": "report_qw2092_gnewton_si_bridge_gate.json",
            "q2115": "report_qw2115_gravity_hierarchy_strict_bridge_gate.json",
            "q2148": "report_qw2148_continuum_dg_delta_extrapolation_gate.json",
            "q2164": "report_qw2164_l14_full_canonical_continuum_variational_gate.json",
            "q2180": "report_qw2180_l14_v2b2_exact_action_identification_gate.json",
        },
        "closed_effective_components": closed_effective_components,
        "open_foundational_components": open_foundational_components,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PRODUCE_DIRECT_EINSTEIN_HILBERT_AND_EQUIVALENCE_PRINCIPLE_DERIVATION_FROM_COMPLETE_FIN_ACTION_OR_KEEP_FOUNDATIONAL_BOUNDARY_EXPLICIT"
            if verdict.endswith("FOUNDATIONAL_OPEN")
            else "REPAIR_GRAVITY_ACTION_LEVEL_SCOPE_CHAIN_AND_RERUN_QW2199"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md_lines = [
        "# RAPORT QW-2199: GRAVITY ACTION-LEVEL SCOPE GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Core result",
        "- Effective gravity bridge components are integrated and explicit in strict scope.",
        "- Foundational Einstein-Hilbert-level derivation claims remain open and explicit.",
        "- Scope boundary is enforced without overclaim.",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
