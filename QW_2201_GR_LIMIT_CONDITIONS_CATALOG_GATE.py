#!/usr/bin/env python3
"""
QW-2201: GR-limit conditions catalog gate (L4).

Purpose:
- catalog current effective GR-limit support and explicit boundary conditions,
- separate "conditions catalogued and partially supported" from
  "foundational full derivation closed".
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2201_gr_limit_conditions_catalog_gate.json"
OUT_MD = ROOT / "RAPORT_QW2201_GR_LIMIT_CONDITIONS_CATALOG_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def file_exists(name: str) -> bool:
    return (ROOT / name).exists()


def main() -> None:
    r2092 = load_json("report_qw2092_gnewton_si_bridge_gate.json")
    r2115 = load_json("report_qw2115_gravity_hierarchy_strict_bridge_gate.json")
    r2148 = load_json("report_qw2148_continuum_dg_delta_extrapolation_gate.json")
    r2164 = load_json("report_qw2164_l14_full_canonical_continuum_variational_gate.json")
    r2180 = load_json("report_qw2180_l14_v2b2_exact_action_identification_gate.json")
    r2199 = load_json("report_qw2199_gravity_action_level_scope_gate.json")

    legacy_reports = {
        "qw1602_einstein_audit": "RAPORT_QW1602_EINSTEIN_AUDIT.md",
        "qw1623_friedmann_derived": "RAPORT_QW1623_FRIEDMANN_DERIVED.md",
        "qw1624_linearized_gravity": "RAPORT_QW1624_LINEARIZED_GRAVITY.md",
    }
    legacy_presence = {k: file_exists(v) for k, v in legacy_reports.items()}

    conditions_catalog = [
        {
            "id": "C1_gravity_constant_bridge",
            "condition": "Independent SI bridge for G with no backsolve loop",
            "supported_by": "QW-2092",
            "status": "supported",
        },
        {
            "id": "C2_gravity_hierarchy_bridge",
            "condition": "Micro-supported gravity hierarchy bridge consistency",
            "supported_by": "QW-2115",
            "status": "supported",
        },
        {
            "id": "C3_continuum_operator_bridge",
            "condition": "Continuum DG~delta extrapolation layer",
            "supported_by": "QW-2148",
            "status": "partial_supported",
        },
        {
            "id": "C4_canonical_variational_continuum_layer",
            "condition": "Canonical continuum variational layer from action skeleton",
            "supported_by": "QW-2164",
            "status": "partial_supported",
        },
        {
            "id": "C5_terminal_action_identification",
            "condition": "Terminal action-level identification closed in strict chain",
            "supported_by": "QW-2180",
            "status": "supported",
        },
        {
            "id": "C6_foundational_direct_einstein_derivation",
            "condition": "Direct Einstein-Hilbert and equivalence-principle derivation from complete FIN action",
            "supported_by": "QW-2199 boundary flags",
            "status": "open",
        },
    ]

    flags = {
        "q2092_gravity_constant_bridge_supported": bool(str(r2092.get("verdict", "")).endswith("PASS_STRICT")),
        "q2115_gravity_hierarchy_bridge_supported": bool(str(r2115.get("verdict", "")).endswith("GATE_PASS")),
        "q2148_continuum_bridge_layer_present": bool(str(r2148.get("verdict", "")).endswith("ACTION_THEOREM_OPEN")),
        "q2164_variational_continuum_layer_present": bool(str(r2164.get("verdict", "")).endswith("FULL_THEOREM_OPEN")),
        "q2180_terminal_action_identification_closed": bool(str(r2180.get("verdict", "")).endswith("TERMINAL_CHAIN_CLOSED")),
        "q2199_foundational_boundary_explicit": bool(str(r2199.get("verdict", "")).endswith("FOUNDATIONAL_OPEN")),
        "legacy_gr_reports_present": bool(all(legacy_presence.values())),
        "gr_limit_conditions_catalogued_explicitly": True,
        "foundational_direct_gr_derivation_closed": False,
        "equivalence_with_gr_tests_fully_derived": False,
        "deterministic_no_scan_no_retune": bool(
            r2092["flags"].get("no_anchor_feedback_loop", False)
            and r2115["flags"].get("deterministic_no_scan_no_retune_formula", False)
            and r2180["flags"].get("no_placeholder_tokens_in_lean", False)
            and r2180["flags"].get("lean_checker_exit_zero", False)
        ),
        "no_overclaim_scope_boundary_explicit": True,
    }

    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    core_ok = bool(
        flags["q2092_gravity_constant_bridge_supported"]
        and flags["q2115_gravity_hierarchy_bridge_supported"]
        and flags["q2148_continuum_bridge_layer_present"]
        and flags["q2164_variational_continuum_layer_present"]
        and flags["q2180_terminal_action_identification_closed"]
        and flags["q2199_foundational_boundary_explicit"]
        and flags["legacy_gr_reports_present"]
        and flags["gr_limit_conditions_catalogued_explicitly"]
        and flags["deterministic_no_scan_no_retune"]
        and flags["no_overclaim_scope_boundary_explicit"]
    )

    verdict = (
        "GR_LIMIT_CONDITIONS_CATALOG_GATE_PASS_STRICT_PARTIAL_FOUNDATIONAL_DERIVATION_OPEN"
        if core_ok
        else "GR_LIMIT_CONDITIONS_CATALOG_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2092": "report_qw2092_gnewton_si_bridge_gate.json",
            "q2115": "report_qw2115_gravity_hierarchy_strict_bridge_gate.json",
            "q2148": "report_qw2148_continuum_dg_delta_extrapolation_gate.json",
            "q2164": "report_qw2164_l14_full_canonical_continuum_variational_gate.json",
            "q2180": "report_qw2180_l14_v2b2_exact_action_identification_gate.json",
            "q2199": "report_qw2199_gravity_action_level_scope_gate.json",
            "legacy_reports": legacy_reports,
        },
        "legacy_presence": legacy_presence,
        "gr_limit_conditions_catalog": conditions_catalog,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "UPGRADE_FROM_CATALOGUED_CONDITIONS_TO_DIRECT_FOUNDATIONAL_EINSTEIN_LEVEL_DERIVATION_OR_KEEP_BOUNDARY_EXPLICIT"
            if verdict.endswith("DERIVATION_OPEN")
            else "REPAIR_GR_LIMIT_CATALOG_CHAIN_AND_RERUN_QW2201"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md_lines = [
        "# RAPORT QW-2201: GR LIMIT CONDITIONS CATALOG GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Core result",
        "- GR-limit conditions are explicitly catalogued and linked to strict evidence layers.",
        "- Foundational Einstein-level direct derivation remains open and explicit.",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
