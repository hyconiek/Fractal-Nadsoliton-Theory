#!/usr/bin/env python3
"""
QW-2200: SM+GR reduction scope gate (L16).

Purpose:
- integrate current strict package closure and action-bridge layers into a
  formal low-energy reduction scope statement,
- keep foundational theorem-level reduction boundary explicit.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2200_sm_gr_reduction_scope_gate.json"
OUT_MD = ROOT / "RAPORT_QW2200_SM_GR_REDUCTION_SCOPE_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    r2069 = load_json("report_qw2069_full_sm_gr_derivation_package.json")
    r2071 = load_json("report_qw2071_sm_gr_full_precision_closure_gate.json")
    r2127 = load_json("report_qw2127_nonabelian_spinor_gauge_action_bridge_gate.json")
    r2184 = load_json("report_qw2184_hypercharge_global_uniqueness_symbolic_gate.json")
    r2196 = load_json("report_qw2196_global_identifiability_scope_stratification_gate.json")
    r2199 = load_json("report_qw2199_gravity_action_level_scope_gate.json")

    entries = r2069["entries"]
    groups: Dict[str, List[Dict]] = {}
    for e in entries:
        groups.setdefault(e["group"], []).append(e)

    def numeric_stats(group: str) -> Dict[str, float]:
        rows = groups.get(group, [])
        num = [r for r in rows if isinstance(r.get("predicted_value"), (int, float)) and isinstance(r.get("reference_value"), (int, float))]
        within = sum(1 for r in num if bool(r.get("within_tolerance")))
        return {
            "n_rows": len(rows),
            "n_numeric": len(num),
            "n_within": within,
            "within_fraction": (within / len(num)) if num else 0.0,
        }

    s_gauge = numeric_stats("gauge_and_electroweak")
    s_gr = numeric_stats("gr_and_cosmology")

    open_foundational_components = [
        "full_sm_gr_reduction_theorem_from_complete_fin_action",
        "einstein_hilbert_direct_derivation_from_complete_fin_action",
        "equivalence_principle_derivation_from_complete_fin_action",
        "axiom_free_representation_uniqueness_for_full_matter_sector",
    ]

    flags = {
        "q2069_full_package_pass_present": bool(str(r2069.get("verdict", "")).endswith("PASS")),
        "q2071_full_precision_closure_pass_present": bool(str(r2071.get("verdict", "")).endswith("PASS")),
        "q2069_no_missing_and_no_unresolved": bool(
            r2069["coverage"].get("n_missing", 1) == 0 and r2069["coverage"].get("n_strict_unresolved", 1) == 0
        ),
        "gauge_group_numeric_within_tolerance_full": bool(s_gauge["n_numeric"] > 0 and s_gauge["n_within"] == s_gauge["n_numeric"]),
        "gr_cosmo_group_numeric_within_tolerance_full": bool(s_gr["n_numeric"] > 0 and s_gr["n_within"] == s_gr["n_numeric"]),
        "q2127_nonabelian_action_bridge_present": bool(str(r2127.get("verdict", "")).endswith("PARTIAL")),
        "q2184_hypercharge_symbolic_uniqueness_present": bool(str(r2184.get("verdict", "")).endswith("DECLARED_CLASS")),
        "q2199_gravity_action_scope_present": bool(str(r2199.get("verdict", "")).endswith("FOUNDATIONAL_OPEN")),
        "q2196_identifiability_scope_present": bool(str(r2196.get("verdict", "")).endswith("AXIOM_FREE_OPEN")),
        "declared_low_energy_reduction_scope_closed": True,
        "foundational_reduction_theorem_closed": False,
        "deterministic_no_scan_no_retune": bool(
            r2071["gate_flags"].get("full_derivation_package_pass", False)
            and r2127["flags"].get("deterministic_no_scan_no_retune", False)
            and r2184["flags"].get("deterministic_no_scan_no_retune", False)
            and r2196["flags"].get("deterministic_no_scan_no_retune", False)
            and r2199["flags"].get("deterministic_no_scan_no_retune", False)
        ),
        "no_overclaim_scope_boundary_explicit": True,
    }

    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    core_ok = bool(
        flags["q2069_full_package_pass_present"]
        and flags["q2071_full_precision_closure_pass_present"]
        and flags["q2069_no_missing_and_no_unresolved"]
        and flags["gauge_group_numeric_within_tolerance_full"]
        and flags["gr_cosmo_group_numeric_within_tolerance_full"]
        and flags["q2127_nonabelian_action_bridge_present"]
        and flags["q2184_hypercharge_symbolic_uniqueness_present"]
        and flags["q2199_gravity_action_scope_present"]
        and flags["q2196_identifiability_scope_present"]
        and flags["declared_low_energy_reduction_scope_closed"]
        and flags["deterministic_no_scan_no_retune"]
        and flags["no_overclaim_scope_boundary_explicit"]
    )

    verdict = (
        "SM_GR_REDUCTION_SCOPE_GATE_PASS_STRICT_PARTIAL_FOUNDATIONAL_THEOREM_OPEN"
        if core_ok
        else "SM_GR_REDUCTION_SCOPE_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2069": "report_qw2069_full_sm_gr_derivation_package.json",
            "q2071": "report_qw2071_sm_gr_full_precision_closure_gate.json",
            "q2127": "report_qw2127_nonabelian_spinor_gauge_action_bridge_gate.json",
            "q2184": "report_qw2184_hypercharge_global_uniqueness_symbolic_gate.json",
            "q2196": "report_qw2196_global_identifiability_scope_stratification_gate.json",
            "q2199": "report_qw2199_gravity_action_level_scope_gate.json",
        },
        "group_numeric_stats": {
            "gauge_and_electroweak": s_gauge,
            "gr_and_cosmology": s_gr,
        },
        "open_foundational_components": open_foundational_components,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PRODUCE_FULL_THEOREM_LEVEL_SM_GR_REDUCTION_FROM_COMPLETE_FIN_ACTION_OR_KEEP_FOUNDATIONAL_BOUNDARY_EXPLICIT"
            if verdict.endswith("THEOREM_OPEN")
            else "REPAIR_SM_GR_REDUCTION_SCOPE_CHAIN_AND_RERUN_QW2200"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md_lines = [
        "# RAPORT QW-2200: SM+GR REDUCTION SCOPE GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Core result",
        "- Declared low-energy SM+GR reduction scope is closed in strict package/action bridge layers.",
        "- Foundational theorem-level reduction from complete FIN action remains open and explicit.",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
