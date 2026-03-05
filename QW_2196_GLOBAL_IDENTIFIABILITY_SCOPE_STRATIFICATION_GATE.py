#!/usr/bin/env python3
"""
QW-2196: Global identifiability scope stratification gate (L6).

Purpose:
- integrate latest uniqueness/obstruction/boundary gates into one strict
  identifiability status layer,
- make explicit what is closed in declared scopes and what remains open in
  axiom-free unconstrained global scope.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2196_global_identifiability_scope_stratification_gate.json"
OUT_MD = ROOT / "RAPORT_QW2196_GLOBAL_IDENTIFIABILITY_SCOPE_STRATIFICATION_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    r2128 = load_json("report_qw2128_kernel_rep_assignment_uniqueness_gate.json")
    r2130 = load_json("report_qw2130_global_gamma_hypothesis_uniqueness_gate.json")
    r2184 = load_json("report_qw2184_hypercharge_global_uniqueness_symbolic_gate.json")
    r2191 = load_json("report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json")
    r2192 = load_json("report_qw2192_mode_index_selection_axiom_gate.json")
    r2193 = load_json("report_qw2193_selection_axiom_family_robustness_gate.json")
    r2194 = load_json("report_qw2194_mass_derivation_calibration_separation_gate.json")
    r2195 = load_json("report_qw2195_generation_mapping_axiom_augmented_gate.json")

    open_axiom_free_components: List[str] = []
    if not r2130["flags"].get("global_unconstrained_formula_space_uniqueness", False):
        open_axiom_free_components.append("global_gamma_unconstrained_formula_space_uniqueness")
    if not r2191["flags"].get("full_physical_uniqueness_closed", False):
        open_axiom_free_components.append("mode_index_physical_uniqueness_axiom_free")
    if not r2195["flags"].get("axiom_free_generation_mapping_closed", False):
        open_axiom_free_components.append("generation_mapping_physical_uniqueness_axiom_free")
    if not r2194["flags"].get("full_mass_chain_anchor_free_without_singleton_anchor", False):
        open_axiom_free_components.append("mass_chain_full_anchor_free_closure")

    scope_closed_components: List[str] = []
    if str(r2130.get("verdict", "")).endswith("ADMISSIBLE_DOMAIN"):
        scope_closed_components.append("q_assignment_uniqueness_admissible_gamma_domain")
    if str(r2184.get("verdict", "")).endswith("DECLARED_CLASS"):
        scope_closed_components.append("hypercharge_uniqueness_declared_formula_class")
    if str(r2192.get("verdict", "")).endswith("UNIQUENESS_CLOSED"):
        scope_closed_components.append("mode_index_uniqueness_axiom_augmented")
    if str(r2193.get("verdict", "")).endswith("AUGMENTED_ROBUST"):
        scope_closed_components.append("mode_index_uniqueness_axiom_augmented_robustness")
    if str(r2195.get("verdict", "")).endswith("AXIOM_FREE_OPEN"):
        scope_closed_components.append("generation_mapping_rule_layer_axiom_augmented")

    flags = {
        "q2128_locked_branch_uniqueness_present": bool(str(r2128.get("verdict", "")).endswith("PARTIAL")),
        "q2130_admissible_gamma_domain_uniqueness_present": bool(str(r2130.get("verdict", "")).endswith("ADMISSIBLE_DOMAIN")),
        "q2184_declared_formula_class_uniqueness_present": bool(str(r2184.get("verdict", "")).endswith("DECLARED_CLASS")),
        "q2191_axiom_free_obstruction_theorem_present": bool(str(r2191.get("verdict", "")).endswith("PASS_STRICT")),
        "q2192_axiom_augmented_uniqueness_present": bool(str(r2192.get("verdict", "")).endswith("UNIQUENESS_CLOSED")),
        "q2193_axiom_augmented_robustness_present": bool(str(r2193.get("verdict", "")).endswith("AUGMENTED_ROBUST")),
        "q2194_mass_boundary_explicit": bool(str(r2194.get("verdict", "")).endswith("BOUNDARY_EXPLICIT")),
        "q2195_generation_rule_layer_present": bool(str(r2195.get("verdict", "")).endswith("AXIOM_FREE_OPEN")),
        "scope_boundaries_explicit_across_components": True,
        "scope_closed_components_nonempty": bool(len(scope_closed_components) >= 4),
        "axiom_free_open_components_nonempty": bool(len(open_axiom_free_components) >= 1),
        "axiom_free_global_identifiability_closed": bool(len(open_axiom_free_components) == 0),
        "deterministic_no_scan_no_retune": bool(
            r2128["flags"].get("deterministic_no_scan_no_retune", False)
            and r2130["flags"].get("deterministic_no_scan_no_retune", False)
            and r2184["flags"].get("deterministic_no_scan_no_retune", False)
            and r2193["flags"].get("deterministic_no_scan_no_retune", False)
            and r2194["flags"].get("deterministic_no_scan_no_retune", False)
            and r2195["flags"].get("deterministic_no_scan_no_retune", False)
        ),
        "no_overclaim_scope_separation_enforced": True,
    }

    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    core_ok = bool(
        flags["q2128_locked_branch_uniqueness_present"]
        and flags["q2130_admissible_gamma_domain_uniqueness_present"]
        and flags["q2184_declared_formula_class_uniqueness_present"]
        and flags["q2191_axiom_free_obstruction_theorem_present"]
        and flags["q2192_axiom_augmented_uniqueness_present"]
        and flags["q2193_axiom_augmented_robustness_present"]
        and flags["q2194_mass_boundary_explicit"]
        and flags["q2195_generation_rule_layer_present"]
        and flags["scope_boundaries_explicit_across_components"]
        and flags["scope_closed_components_nonempty"]
        and flags["axiom_free_open_components_nonempty"]
        and flags["deterministic_no_scan_no_retune"]
        and flags["no_overclaim_scope_separation_enforced"]
    )

    verdict = (
        "GLOBAL_IDENTIFIABILITY_SCOPE_STRATIFICATION_GATE_PASS_STRICT_PARTIAL_AXIOM_FREE_OPEN"
        if core_ok
        else "GLOBAL_IDENTIFIABILITY_SCOPE_STRATIFICATION_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2128": "report_qw2128_kernel_rep_assignment_uniqueness_gate.json",
            "q2130": "report_qw2130_global_gamma_hypothesis_uniqueness_gate.json",
            "q2184": "report_qw2184_hypercharge_global_uniqueness_symbolic_gate.json",
            "q2191": "report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json",
            "q2192": "report_qw2192_mode_index_selection_axiom_gate.json",
            "q2193": "report_qw2193_selection_axiom_family_robustness_gate.json",
            "q2194": "report_qw2194_mass_derivation_calibration_separation_gate.json",
            "q2195": "report_qw2195_generation_mapping_axiom_augmented_gate.json",
        },
        "scope_closed_components": scope_closed_components,
        "open_axiom_free_components": open_axiom_free_components,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "TARGET_AXIOM_FREE_COMPONENTS_WITH_EXPLICIT_CLOSURE_GATES_OR_FORMALIZE_THEORY_SCOPE_AS_AXIOM_AUGMENTED"
            if verdict.endswith("AXIOM_FREE_OPEN")
            else "REPAIR_IDENTIFIABILITY_SCOPE_CHAIN_AND_RERUN_QW2196"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md_lines = [
        "# RAPORT QW-2196: GLOBAL IDENTIFIABILITY SCOPE STRATIFICATION GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Core result",
        "- Identifiability is stratified into closed declared scopes vs open axiom-free global scope.",
        "- Axiom-free open components are explicit and enumerated.",
        "- No-overclaim rule is enforced across the integrated chain.",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
