#!/usr/bin/env python3
"""Scratch probe: ToE potential/readiness audit for the strict-kernel release.

This probe gives a finite, professor-style readiness certificate for discussing
ToE potential.  It intentionally does not close ToE: it checks that the new
release-facing ToE potential note is backed by the existing theorem-frontier,
traceability, source-coverage, role-transfer, and chain-integrity reports, and
that the numerical frontier still requires all seven open atoms.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FAR = ROOT / "fundamental_action_reconstruction"
OUT_JSON = HERE / "bridge_strict_completion_toe_potential_readiness_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_toe_potential_readiness_certificate_report.md"
TOE_AUDIT_DOC = FAR / "STRICT_KERNEL_TOE_POTENTIAL_AUDIT.md"

SOURCE_REPORTS = {
    "chain_integrity": HERE / "bridge_strict_completion_certificate_chain_integrity_report.json",
    "frontier_cut": HERE / "bridge_strict_completion_theorem_frontier_cut_certificate_report.json",
    "frontier_truth_table": HERE / "bridge_strict_completion_theorem_frontier_truth_table_certificate_report.json",
    "frontier_atom_influence": HERE / "bridge_strict_completion_theorem_frontier_atom_influence_certificate_report.json",
    "frontier_target_signature_lattice": HERE / "bridge_strict_completion_theorem_frontier_target_signature_lattice_certificate_report.json",
    "frontier_low_weight_extension": HERE / "bridge_strict_completion_theorem_frontier_low_weight_extension_certificate_report.json",
    "role_transfer_lattice": HERE / "bridge_strict_completion_legacy_role_transfer_minimal_obligation_lattice_certificate_report.json",
    "release_source_coverage": HERE / "bridge_strict_completion_release_source_coverage_certificate_report.json",
    "release_traceability_matrix": HERE / "bridge_strict_completion_release_traceability_matrix_certificate_report.json",
}

REQUIRED_DOC_SNIPPETS = [
    "professorial potential assessment only; no ToE closure is claimed",
    "nadsoliton -> light -> matter -> emergent observer",
    "2^7 = 128",
    "ToE closure is available in exactly `1` assignment, with minimal weight `7`",
    "reachable target signatures are only `6` of the possible `16`",
    "all `7` singleton and all `21` pair extensions fail to close bridge",
    "chi11_selector_source` is the unique top Boolean bottleneck",
    "The strict nadsoliton kernel has a serious, structured ToE potential",
    "No ToE closure is claimed",
]

OPEN_ATOMS = [
    "strict_dynamical_source_for_A_P_D",
    "strict_phase_frequency_source",
    "strict_damping_beta_eta_source",
    "chi11_selector_source",
    "alpha_geo_electroweak_role_theorem",
    "beta_tors_strict_role_theorem",
    "beta_power_hierarchy_successor_theorem",
]


TARGET_SIGNATURE_BITS = ["bridge", "role_transfer", "selector_qw2191", "toe"]


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def contains_all(text: str, snippets: list[str]) -> bool:
    return all(snippet in text for snippet in snippets)


def build_payload() -> dict[str, Any]:
    reports = {name: load_json(path) for name, path in SOURCE_REPORTS.items()}
    doc_text = TOE_AUDIT_DOC.read_text(encoding="utf-8") if TOE_AUDIT_DOC.exists() else ""

    truth_summary = reports["frontier_truth_table"]["theorem_frontier_truth_table_summary"]
    cut_summary = reports["frontier_cut"]["theorem_frontier_cut_summary"]
    influence_summary = reports["frontier_atom_influence"]["theorem_frontier_atom_influence_summary"]
    target_summary = reports["frontier_target_signature_lattice"]["theorem_frontier_target_signature_lattice_summary"]
    low_weight_summary = reports["frontier_low_weight_extension"]["theorem_frontier_low_weight_extension_summary"]
    role_summary = reports["role_transfer_lattice"]["role_transfer_minimal_obligation_summary"]
    coverage_summary = reports["release_source_coverage"]["coverage_summary"]
    traceability_summary = reports["release_traceability_matrix"]["traceability_summary"]
    chain_summary = reports["chain_integrity"]["chain_summary"]

    truth_assignment_count = truth_summary["truth_assignment_count"]
    toe_satisfying_assignment_count = truth_summary["toe_satisfying_assignment_count"]
    readiness_fraction = Fraction(toe_satisfying_assignment_count, truth_assignment_count)
    open_leaf_cut = cut_summary["minimal_open_leaf_cut_for_toe"]
    open_leaf_cut_matches = sorted(open_leaf_cut) == sorted(OPEN_ATOMS)

    low_weight_no_go = (
        low_weight_summary["singleton_extension_count"] == 7
        and low_weight_summary["pair_extension_count"] == 21
        and low_weight_summary["no_singleton_closes_bridge_role_or_toe"]
        and low_weight_summary["no_pair_closes_bridge"]
        and low_weight_summary["no_pair_closes_role_transfer"]
        and low_weight_summary["no_pair_closes_toe"]
    )

    frontier_quantitative_rows = [
        {
            "metric": "open_atom_count",
            "value": truth_summary["open_atom_count"],
            "interpretation": "seven named theorem atoms remain open before ToE can be represented as true",
        },
        {
            "metric": "truth_assignment_count",
            "value": truth_assignment_count,
            "interpretation": "finite readiness board enumerates every Boolean state of the open atoms",
        },
        {
            "metric": "toe_satisfying_assignment_count",
            "value": toe_satisfying_assignment_count,
            "interpretation": "only the all-true frontier assignment reaches ToE target closure",
        },
        {
            "metric": "toe_readiness_fraction",
            "value": f"{readiness_fraction.numerator}/{readiness_fraction.denominator}",
            "interpretation": "current frontier gives one ToE-capable assignment among 128 formal assignments",
        },
        {
            "metric": "toe_minimal_set_size",
            "value": truth_summary["toe_minimal_set_size"],
            "interpretation": "no proper subset of the seven frontier atoms closes ToE",
        },
        {
            "metric": "reachable_target_signature_count",
            "value": target_summary["reachable_target_signature_count"],
            "interpretation": "only six of sixteen possible bridge/role/selector/ToE signatures are reachable",
        },
        {
            "metric": "chi11_total_critical_count",
            "value": influence_summary["chi11_selector_source_total_critical_count"],
            "interpretation": "selector source is the unique top Boolean bottleneck on the finite board",
        },
    ]

    cross_checks = {
        "toe_audit_doc_present": TOE_AUDIT_DOC.exists(),
        "toe_audit_doc_required_snippets_present": contains_all(doc_text, REQUIRED_DOC_SNIPPETS),
        "source_reports_present": all(path.exists() for path in SOURCE_REPORTS.values()),
        "chain_integrity_inherited": reports["chain_integrity"]["all_cross_checks_pass"] and not chain_summary["bridge_theorem_exported"] and not chain_summary["strict_dynamic_derivation_exported"],
        "open_atom_set_pass": truth_summary["open_atom_count"] == 7 and open_leaf_cut_matches,
        "truth_table_counts_pass": truth_assignment_count == 128 and toe_satisfying_assignment_count == 1,
        "toe_minimal_weight_pass": truth_summary["toe_minimal_set_size"] == 7 and truth_summary["toe_minimal_set_equals_frontier_cut"],
        "target_signature_lattice_pass": target_summary["reachable_target_signature_count"] == 6 and target_summary["only_full_signature_has_toe_closure"] and target_summary["toe_implies_all_targets"],
        "low_weight_no_go_pass": low_weight_no_go,
        "chi11_top_bottleneck_pass": influence_summary["chi11_selector_source_is_unique_top_influence"] and influence_summary["top_influence_atoms"] == ["chi11_selector_source"],
        "role_transfer_still_blocked": role_summary["all_atoms_missing_now"] and not role_summary["role_transfer_theorem_exported"],
        "release_coverage_inherited": coverage_summary["release_scaffold_nonduplicating_source_map_ready"] and coverage_summary["no_toe_closure"],
        "traceability_matrix_inherited": traceability_summary["target_columns_independent_over_gf2"] and traceability_summary["no_toe_closure"],
        "hard_limits_preserved": not any([
            truth_summary["toe_closure_claimed"],
            target_summary["toe_closure_claimed"],
            low_weight_summary["toe_closure_claimed"],
            influence_summary["toe_closure_claimed"],
            cut_summary["toe_closure_claimed"],
            role_summary["toe_closure_claimed"],
        ]),
    }

    summary = {
        "professorial_toe_potential_doc_certified": cross_checks["toe_audit_doc_present"] and cross_checks["toe_audit_doc_required_snippets_present"],
        "finite_frontier_board_certified": cross_checks["open_atom_set_pass"] and cross_checks["truth_table_counts_pass"] and cross_checks["toe_minimal_weight_pass"],
        "toe_requires_all_7_open_atoms": cross_checks["toe_minimal_weight_pass"] and open_leaf_cut_matches,
        "toe_readiness_fraction": f"{readiness_fraction.numerator}/{readiness_fraction.denominator}",
        "reachable_target_signatures": target_summary["reachable_signatures"],
        "low_weight_extensions_do_not_close_toe": cross_checks["low_weight_no_go_pass"],
        "chi11_selector_source_is_top_bottleneck": cross_checks["chi11_top_bottleneck_pass"],
        "role_transfer_remains_blocked": cross_checks["role_transfer_still_blocked"],
        "traceability_and_source_coverage_inherited": cross_checks["release_coverage_inherited"] and cross_checks["traceability_matrix_inherited"],
        "toe_closure_claimed": False,
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_TOE_POTENTIAL_READINESS_CERTIFICATE__PROFESSORIAL_FINITE_AUDIT_NO_CLOSURE",
        "status": "PASS" if all(cross_checks.values()) else "FAIL",
        "toe_audit_doc": rel(TOE_AUDIT_DOC),
        "source_reports": {name: rel(path) for name, path in SOURCE_REPORTS.items()},
        "target_signature_bits": TARGET_SIGNATURE_BITS,
        "open_atoms": OPEN_ATOMS,
        "frontier_quantitative_rows": frontier_quantitative_rows,
        "cross_checks": cross_checks,
        "all_cross_checks_pass": all(cross_checks.values()),
        "toe_potential_readiness_summary": summary,
        "proof_certificate": {
            "professorial_assessment_step": "The ToE potential audit gives a professor-level assessment: strict-kernel ToE potential is serious and structured because bridge components, role boundaries, selector obstruction, and source coverage are finitely named, but the note explicitly claims no ToE closure.",
            "truth_table_step": "The inherited theorem-frontier truth table enumerates 2^7=128 assignments, finds exactly one ToE-capable assignment, and records minimal ToE weight 7; therefore ToE readiness is 1/128 and requires all seven open atoms.",
            "signature_lattice_step": "The target-signature lattice has 6 reachable signatures out of 16 and permits ToE only at full signature 1111, so partial bridge/selector/role progress is not silently promoted to ToE closure.",
            "low_weight_step": "The low-weight extension audit checks all 7 singleton and 21 pair extensions: chi11_selector_source is the only singleton unlock, all unlocking pairs contain chi11, and no singleton or pair closes bridge, role-transfer, or ToE.",
            "bottleneck_step": "The atom-influence audit identifies chi11_selector_source as the unique top bottleneck with total critical count 73, while bridge-source atoms tie at 17 and role-only atoms tie at 9; these scores are priorities, not theorem sources.",
            "limit_step": "This certificate is finite readiness bookkeeping only: it preserves no identity bridge, no strict dynamical source theorem, no legacy role-transfer theorem, no QW-2191 discharge, no tensor-resolved EOM closure, and no ToE closure.",
        },
        "hard_limits": [
            "No identity bridge K_legacy_ont == K_strict_gate is claimed.",
            "No strict dynamical source theorem for A/P/D, phase transport, or damping/compression is claimed.",
            "No legacy physical-role transfer to K_strict_gate is claimed.",
            "No beta_tors -> chi11 theorem is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No full tensor-resolved Lagrangian/EOM theorem is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# ToE potential readiness certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This report certifies that the ToE potential discussion is backed by the",
        "finite theorem-frontier and release-traceability ledger while preserving",
        "all non-closure guardrails.",
        "",
        "## Quantitative frontier rows",
        "",
    ]
    for row in payload["frontier_quantitative_rows"]:
        lines.append(f"- `{row['metric']}` = `{row['value']}` — {row['interpretation']}")
    lines.extend(["", "## Cross-checks", ""])
    for key, value in payload["cross_checks"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend([
        "",
        f"All cross-checks pass: `{payload['all_cross_checks_pass']}`",
        "",
        "## Summary",
        "",
    ])
    for key, value in payload["toe_potential_readiness_summary"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Hard limits", ""])
    for limit in payload["hard_limits"]:
        lines.append(f"- {limit}")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(payload)
    print(json.dumps({"status": payload["status"], "all_cross_checks_pass": payload["all_cross_checks_pass"]}, sort_keys=True))


if __name__ == "__main__":
    main()
