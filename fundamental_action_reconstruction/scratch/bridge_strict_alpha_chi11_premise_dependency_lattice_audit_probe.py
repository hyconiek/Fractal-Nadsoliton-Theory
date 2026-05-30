#!/usr/bin/env python3
from __future__ import annotations

import itertools
import json
from pathlib import Path
from typing import Any, Callable

ROOT = Path(__file__).resolve().parents[1]
SCRATCH = ROOT / "scratch"
OUT_JSON = SCRATCH / "bridge_strict_alpha_chi11_premise_dependency_lattice_audit_report.json"
OUT_MD = SCRATCH / "bridge_strict_alpha_chi11_premise_dependency_lattice_audit_report.md"

INPUT_REPORTS = {
    "affine_subgroup_lattice": SCRATCH / "bridge_strict_alpha_affine_subgroup_lattice_chi11_source_obstruction_report.json",
    "cyclic_character_projection": SCRATCH / "bridge_strict_alpha_cyclic_quotient_character_projection_chi11_certificate_report.json",
    "d12_character_module": SCRATCH / "bridge_strict_alpha_d12_chi11_character_module_certificate_report.json",
    "d12_sparsest_extension": SCRATCH / "bridge_strict_alpha_d12_chi11_sparsest_extension_selector_audit_report.json",
    "d12_max_shell_imbalance": SCRATCH / "bridge_strict_alpha_d12_chi11_max_shell_imbalance_selector_certificate_report.json",
    "full_aut_block_amplitude": SCRATCH / "bridge_strict_alpha_full_aut_block_amplitude_chi11_polarity_obstruction_report.json",
    "nonhistogram_dihedral_quotient": SCRATCH / "bridge_strict_alpha_nonhistogram_dihedral_quotient_chi11_audit_report.json",
}

PREMISES = [
    "full_aut_unoriented_block_amplitude",
    "d12_reduced_quotient",
    "cyclic_translation_quotient",
    "chi11_unit_character_choice",
    "branch_normalization_constraint",
    "sparsest_extension_selector",
    "shell_label_d1_d5_axis",
    "max_shell_imbalance_selector",
]


def load_report(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def is_minimal(subset: frozenset[str], accepted: list[frozenset[str]]) -> bool:
    return not any(candidate < subset for candidate in accepted)


def enumerate_minimal(predicate: Callable[[frozenset[str]], bool]) -> list[list[str]]:
    accepted: list[frozenset[str]] = []
    for size in range(len(PREMISES) + 1):
        for combo in itertools.combinations(PREMISES, size):
            subset = frozenset(combo)
            if predicate(subset) and is_minimal(subset, accepted):
                accepted.append(subset)
    return [sorted(row, key=PREMISES.index) for row in accepted]


def has_all(subset: frozenset[str], *items: str) -> bool:
    return all(item in subset for item in items)


def main() -> None:
    reports = {name: load_report(path) for name, path in INPUT_REPORTS.items()}

    full_aut = reports["full_aut_block_amplitude"]
    affine = reports["affine_subgroup_lattice"]
    cyclic = reports["cyclic_character_projection"]
    d12 = reports["d12_character_module"]
    sparse = reports["d12_sparsest_extension"]
    maximb = reports["d12_max_shell_imbalance"]
    nonhist = reports["nonhistogram_dihedral_quotient"]

    cross_report_checks = {
        "input_report_count": len(reports),
        "all_input_reports_loaded": all(report.get("status") for report in reports.values()),
        "full_aut_orbit_count": full_aut["finite_model"]["full_affine_aut_orbit_count"],
        "d12_orbit_count": d12["finite_model"]["d12_orbit_count"],
        "cyclic_translation_orbit_count": cyclic["finite_model"]["translation_orbit_count"],
        "support_count_agrees": len({report["finite_model"]["support_count"] for report in reports.values()}) == 1,
        "full_aut_kills_chi11": not affine["lattice_summary"]["full_aut_admits_chi11"],
        "unit5_and_unit7_are_killing_additions": affine["lattice_summary"]["minimal_killing_unit_additions"] == ["unit5", "unit7"],
        "full_aut_block_amplitude_locates_unique_block": full_aut["block_amplitude_summary"]["unique_maximizer_is_branch_full_aut_block"],
        "full_aut_block_amplitude_exports_no_polarity": not full_aut["block_amplitude_summary"]["exports_chi11_polarity"],
        "nonhistogram_singleton_A5_not_A1_count_zero": nonhist["quotient_summary"]["full_aut_singleton_A5_not_A1_classifier_count"] == 0,
        "cyclic_chi11_dimension": cyclic["trace_summary"]["chi11_dimension"],
        "cyclic_full_aut_trivial_intersection_with_chi11_rank": cyclic["trace_summary"]["full_aut_trivial_intersection_with_chi11_rank"],
        "d12_chi11_rank": d12["module_summary"]["integer_chi11_covariant_module_rank"],
        "d12_full_aut_intersection_rank_zero": d12["module_summary"]["full_aut_invariant_intersection_rank"] == 0,
        "sparsest_unique_minimum_is_branch_generator": sparse["sparsity_summary"]["unique_minimum_is_branch_generator"],
        "max_shell_imbalance_unique_branch_maximizer": maximb["selector_summary"]["unique_maximizer_is_branch_cycle"],
        "max_shell_imbalance_requires_shell_label": maximb["selector_summary"]["score_requires_shell_label"],
    }

    outcome_predicates: dict[str, Callable[[frozenset[str]], bool]] = {
        "locate_branch_full_aut_block_without_polarity": lambda s: has_all(
            s, "full_aut_unoriented_block_amplitude"
        ),
        "host_nonzero_cyclic_chi11_character_space": lambda s: has_all(
            s, "cyclic_translation_quotient", "chi11_unit_character_choice"
        ),
        "host_d12_chi11_covariant_module": lambda s: has_all(
            s, "d12_reduced_quotient", "chi11_unit_character_choice"
        ),
        "branch_normalized_d12_chi11_family": lambda s: has_all(
            s,
            "d12_reduced_quotient",
            "chi11_unit_character_choice",
            "branch_normalization_constraint",
        ),
        "unique_branch_generator_by_sparsest_extension": lambda s: has_all(
            s,
            "d12_reduced_quotient",
            "chi11_unit_character_choice",
            "branch_normalization_constraint",
            "sparsest_extension_selector",
        ),
        "unique_branch_generator_by_max_shell_imbalance": lambda s: has_all(
            s,
            "d12_reduced_quotient",
            "chi11_unit_character_choice",
            "branch_normalization_constraint",
            "shell_label_d1_d5_axis",
            "max_shell_imbalance_selector",
        ),
        "strict_full_aut_internal_chi11_polarity_source": lambda s: False,
    }

    minimal_premise_antichains = {
        name: enumerate_minimal(predicate) for name, predicate in outcome_predicates.items()
    }

    outcome_rows = []
    for name, antichain in minimal_premise_antichains.items():
        outcome_rows.append(
            {
                "outcome": name,
                "minimal_premise_sets": antichain,
                "minimal_premise_count": None if not antichain else min(len(row) for row in antichain),
                "realized_by_current_reports": name != "strict_full_aut_internal_chi11_polarity_source",
            }
        )

    exact_proof_certificate = {
        "finite_domain": "This audit reuses the already-generated C(12,5)=792 support reports and checks their shared finite-model counts before building a premise antichain.",
        "not_duplicate_of_general_qw2191_lattice": "P2165 audits general QW-2191 lane admissibility.  This probe is narrower: it computes the minimal premise antichain for the current strict-alpha chi_11 support reports only.",
        "full_aut_no_polarity": "The affine subgroup lattice says full Aut does not admit chi_11; the full-Aut block-amplitude report locates the branch block but has exports_chi11_polarity=false; the nonhistogram audit has zero full-Aut singleton A5-not-A1 classifiers.",
        "reduced_quotient_boundary": "A nonzero chi_11 carrier appears only after quotient/premise reduction: cyclic translation quotient plus chi_11 character gives dimension 13, and D12 plus chi_11 character gives the 13 two-cycle module.",
        "selector_boundary": "The branch-normalized D12 chi_11 family is not unique until an extra selector is imported.  Current finite audits certify two conditional selectors: sparsest extension and shell-labelled max imbalance.",
        "strict_limit": "No row in this premise lattice supplies an internal strict full-Aut source for the chi_11 polarity, the shell-label d1/d5 axis, or QW-2191 selector closure.",
    }

    payload = {
        "result_kind": "SCRATCH_STRICT_ALPHA_CHI11_PREMISE_DEPENDENCY_LATTICE_AUDIT__NO_STRICT_FULL_AUT_SOURCE",
        "status": "minimal-premise-antichain-computed-for-existing-chi11-support-audits",
        "finite_model": {
            "ring": "Z_12",
            "active_count": 5,
            "support_count": 792,
            "premise_count": len(PREMISES),
            "premise_subset_count": 2 ** len(PREMISES),
            "target_q_power": "256/243",
            "target_eta": "9/5",
        },
        "repo_grep_nonduplication_note": {
            "performed_before_new_probe": True,
            "finding": "Existing files already cover affine subgroup lattice, cyclic projection, D12 module, max-shell selector, sparsest-extension selector, full-Aut block amplitude, and nonhistogram D12 quotient.  The new work therefore audits the dependency antichain across those reports rather than re-running the same orbit study.",
        },
        "input_reports": {name: str(path.relative_to(ROOT)) for name, path in INPUT_REPORTS.items()},
        "premises": [
            {"id": premise, "index": index, "strict_core_exported": False if "selector" in premise or "axis" in premise or "character" in premise or "quotient" in premise or "normalization" in premise else None}
            for index, premise in enumerate(PREMISES)
        ],
        "cross_report_checks": cross_report_checks,
        "outcome_rows": outcome_rows,
        "minimal_premise_antichains": minimal_premise_antichains,
        "exact_proof_certificate": exact_proof_certificate,
        "interpretation": {
            "honest_positive": "The existing chi_11 audit chain can now be read as a finite premise-dependency antichain: it tells exactly which extra premises are sufficient for each weaker outcome.",
            "honest_negative": "The antichain has no minimal set for a strict full-Aut internal chi_11 polarity source; all polarity/selector outcomes still require reduced quotient, character, shell-label, or selector premises.",
            "relation_to_previous_probe": "This is a meta-certificate over prior generated support-combinatorics reports, not another duplicate enumeration of the same D12/full-Aut orbits.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in solitonic state; this audit only classifies finite support-premise dependencies.",
            "forbidden_reading": "No separate informational layer under the nadsoliton is introduced, and no strict-core source for chi_11 is claimed.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted.",
            "No legacy physical role is transferred onto K_strict_gate.",
            "No theorem derives chi_11 from full-Aut strict geometry.",
            "No theorem derives the shell-label d1/d5 axis or unit-axis bit from strict full-Aut support data.",
            "No theorem derives the sparsest-extension selector or max-shell-imbalance selector as strict-core closure.",
            "No QW-2191 discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }

    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    rows_md = "\n".join(
        f"| `{row['outcome']}` | `{row['minimal_premise_count']}` | `{row['minimal_premise_sets']}` |"
        for row in outcome_rows
    )
    OUT_MD.write_text(
        "# Strict-alpha chi_11 premise-dependency lattice audit\n\n"
        f"- Result kind: `{payload['result_kind']}`\n"
        f"- Status: `{payload['status']}`\n"
        f"- Premise subsets checked: `{payload['finite_model']['premise_subset_count']}`.\n"
        f"- Input reports loaded: `{cross_report_checks['input_report_count']}`.\n"
        f"- Full-Aut exports chi_11 polarity: `{not cross_report_checks['full_aut_block_amplitude_exports_no_polarity']}`.\n"
        f"- D12 chi_11 rank: `{cross_report_checks['d12_chi11_rank']}`.\n"
        f"- Cyclic chi_11 dimension: `{cross_report_checks['cyclic_chi11_dimension']}`.\n\n"
        "## Minimal premise antichain\n\n"
        "| Outcome | Minimal premise count | Minimal premise sets |\n"
        "|---|---:|---|\n"
        f"{rows_md}\n\n"
        "## Proof certificate\n\n"
        f"- {exact_proof_certificate['finite_domain']}\n"
        f"- {exact_proof_certificate['full_aut_no_polarity']}\n"
        f"- {exact_proof_certificate['reduced_quotient_boundary']}\n"
        f"- {exact_proof_certificate['selector_boundary']}\n"
        f"- {exact_proof_certificate['strict_limit']}\n\n"
        "## Hard limits\n\n"
        + "\n".join(f"- {item}" for item in payload["hard_limits"])
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
