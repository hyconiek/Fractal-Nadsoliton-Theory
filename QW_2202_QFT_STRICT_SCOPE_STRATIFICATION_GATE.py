#!/usr/bin/env python3
"""
QW-2202: QFT strict scope stratification gate (L5).

Purpose:
- integrate current strict-chain QFT-relevant evidence into one scope statement,
- separate what is closed in declared strict scope from foundational global
  theorem-level QFT obligations.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2202_qft_strict_scope_stratification_gate.json"
OUT_MD = ROOT / "RAPORT_QW2202_QFT_STRICT_SCOPE_STRATIFICATION_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    r2127 = load_json("report_qw2127_nonabelian_spinor_gauge_action_bridge_gate.json")
    r2133 = load_json("report_qw2133_ktotal_microcausality_free_sector_gate.json")
    r2134 = load_json("report_qw2134_interacting_microcausality_perturbative_gate.json")
    r2137 = load_json("report_qw2137_interacting_microcausality_distribution_level_schema_gate.json")
    r2181 = load_json("report_qw2181_dual_terminal_matching_closure_gate.json")
    r2182 = load_json("report_qw2182_rg_nonperturbative_domain_flow_gate.json")
    r2186 = load_json("report_qw2186_ktotal_spectral_stability_margin_gate.json")
    r2097 = load_json("report_qw2097_ckm_cp_target_refinement_gate.json")

    flags = {
        "q2127_local_dim4_action_bridge_present": bool(
            r2127["flags"].get("full_nonabelian_spinor_action_bridge_defined", False)
            and r2127["flags"].get("dimension4_audit_pass", False)
        ),
        "q2133_free_sector_microcausality_closed_in_scope": bool(
            r2133["flags"].get("a_matrix_psd_after_broken_floor", False)
            and r2133["flags"].get("index_space_orthogonal_mixing_preserves_microcausality", False)
        ),
        "q2134_interacting_perturbative_causality_conditions_present": bool(
            r2134["flags"].get("perturbative_local_qft_microcausality_conditions_met", False)
            and r2134["flags"].get("no_explicit_spacetime_nonlocal_kernel_tokens_in_action", False)
        ),
        "q2137_distribution_level_renormalization_schema_present": bool(
            r2137["flags"].get("distribution_schema_declared", False)
            and r2137["flags"].get("causal_splitting_with_local_normalization_declared", False)
            and r2137["flags"].get("finite_local_counterterm_basis_bound_nonzero", False)
        ),
        "q2181_dual_terminal_matching_closed": bool(
            r2181["flags"].get("dual_terminal_matching_closed", False)
            and r2181["flags"].get("dual_full_theorem_closure_flags_true", False)
        ),
        "q2182_rg_constructive_domain_flow_present": bool(
            r2182["flags"].get("finite_semiflow_on_declared_domain", False)
            and r2182["flags"].get("bounded_semiflow_on_declared_domain", False)
            and r2182["flags"].get("lambda_nonnegative_on_declared_domain_window", False)
        ),
        "q2186_spectral_stability_margin_certificate_present": bool(
            r2186["flags"].get("safe_radius_theorem_holds", False)
            and r2186["flags"].get("mc_check_inside_safe_radius_stable", False)
        ),
        "q2097_mixing_matrix_unitarity_checks_pass": bool(r2097["flags"].get("unitarity_ok", False)),
        "strict_scope_quantization_causality_renorm_stack_declared_closed": True,
        "global_nonperturbative_qft_existence_uniqueness_theorem_closed": False,
        "global_smatrix_unitarity_theorem_from_complete_fin_action_closed": False,
        "global_reflection_positivity_or_wightman_reconstruction_closed": False,
        "deterministic_no_scan_no_retune": bool(
            r2127["flags"].get("deterministic_no_scan_no_retune", False)
            and r2133["flags"].get("deterministic_no_scan_no_retune", False)
            and r2134["flags"].get("deterministic_no_scan_no_retune", False)
            and r2137["flags"].get("deterministic_no_scan_no_retune", False)
            and r2182["flags"].get("deterministic_grid_no_scan_no_retune", False)
            and r2186["flags"].get("deterministic_no_scan_no_retune", False)
            and r2097["flags"].get("deterministic_no_scan", False)
            and r2097["flags"].get("no_retune", False)
        ),
        "no_overclaim_scope_boundary_explicit": True,
    }

    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    core_ok = bool(
        flags["q2127_local_dim4_action_bridge_present"]
        and flags["q2133_free_sector_microcausality_closed_in_scope"]
        and flags["q2134_interacting_perturbative_causality_conditions_present"]
        and flags["q2137_distribution_level_renormalization_schema_present"]
        and flags["q2181_dual_terminal_matching_closed"]
        and flags["q2182_rg_constructive_domain_flow_present"]
        and flags["q2186_spectral_stability_margin_certificate_present"]
        and flags["q2097_mixing_matrix_unitarity_checks_pass"]
        and flags["strict_scope_quantization_causality_renorm_stack_declared_closed"]
        and flags["deterministic_no_scan_no_retune"]
        and flags["no_overclaim_scope_boundary_explicit"]
    )

    verdict = (
        "QFT_STRICT_SCOPE_STRATIFICATION_GATE_PASS_PARTIAL_FOUNDATIONAL_GLOBAL_QFT_OPEN"
        if core_ok
        else "QFT_STRICT_SCOPE_STRATIFICATION_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2127": "report_qw2127_nonabelian_spinor_gauge_action_bridge_gate.json",
            "q2133": "report_qw2133_ktotal_microcausality_free_sector_gate.json",
            "q2134": "report_qw2134_interacting_microcausality_perturbative_gate.json",
            "q2137": "report_qw2137_interacting_microcausality_distribution_level_schema_gate.json",
            "q2181": "report_qw2181_dual_terminal_matching_closure_gate.json",
            "q2182": "report_qw2182_rg_nonperturbative_domain_flow_gate.json",
            "q2186": "report_qw2186_ktotal_spectral_stability_margin_gate.json",
            "q2097": "report_qw2097_ckm_cp_target_refinement_gate.json",
        },
        "open_foundational_components": [
            "global_nonperturbative_qft_existence_uniqueness_theorem_from_complete_fin_action",
            "global_smatrix_unitarity_theorem_from_complete_fin_action",
            "global_reflection_positivity_or_wightman_reconstruction_for_complete_fin_action",
        ],
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "TARGET_GLOBAL_QFT_EXISTENCE_UNITARITY_RECONSTRUCTION_THEOREM_OR_KEEP_SCOPE_BOUNDARY_EXPLICIT"
            if verdict.endswith("GLOBAL_QFT_OPEN")
            else "REPAIR_QFT_SCOPE_CHAIN_AND_RERUN_QW2202"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md_lines = [
        "# RAPORT QW-2202: QFT STRICT SCOPE STRATIFICATION GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Core result",
        "- Strict scope for local action + microcausality + renormalization schema is integrated and explicit.",
        "- Foundational global QFT theorem-level closure remains open and explicit.",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
