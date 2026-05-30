#!/usr/bin/env python3
"""Scratch probe: completed-kernel bridge blocker dependency lattice.

The repo now contains many legacy->strict bridge-side probes.  This script is a
non-duplicating next step: it does not recompute the kernel values, chi_11
orbits, or Reynolds matrices.  Instead it loads the existing reports and builds
an exact finite dependency lattice for what is already certified versus what is
still a blocker for a theorem-level statement that the current strict kernel is
the completed legacy/nadsoliton-characteristic carrier.

No false pass: exact factorization and local/measure candidates are separated
from strict derivations, global Z12 transport, the beta_tors->chi_11 source,
chi_11 uniqueness, Reynolds-obstruction escape, QW-2191 discharge, and ToE
closure.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path
from typing import Any, Callable

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_completed_kernel_blocker_dependency_lattice_report.json"
OUT_MD = HERE / "bridge_completed_kernel_blocker_dependency_lattice_report.md"

REPORTS = {
    "completed_factorization": HERE / "bridge_completed_strict_kernel_factorization_certificate_report.json",
    "reactivation": HERE / "bridge_legacy_kernel_reactivation_from_diagrams_candidate_report.json",
    "torsion_opinion": HERE / "bridge_legacy_torsion_chi11_opinion_audit_report.json",
    "reynolds": HERE / "bridge_strict_alpha_reynolds_annihilator_chi11_matrix_certificate_report.json",
    "puiseux_map": HERE / "bridge_phase_normalized_puiseux_map_report.json",
    "puiseux_eta_plus_one": HERE / "bridge_phase_normalized_puiseux_eta_plus_one_report.json",
    "measure_transport": HERE / "bridge_phase_normalized_measure_transport_report.json",
    "inverse_branch_ode": HERE / "bridge_phase_normalized_inverse_branch_ode_report.json",
    "global_admissibility": HERE / "bridge_phase_normalized_global_admissibility_report.json",
}

PREMISES = [
    "diagrams_legacy_carrier",
    "exact_completion_factorization",
    "local_puiseux_match",
    "eta_plus_one_puiseux_match",
    "measure_transport_identity",
    "monotone_flow_output_matching",
    "strict_transport_derivation",
    "global_z12_map_derivation",
    "orientation_chi11_source",
    "chi11_uniqueness",
    "reynolds_obstruction_escape",
    "role_transfer_theorem",
    "qw2191_selector_discharge",
]

EXPORTED_PREMISES = {
    "diagrams_legacy_carrier",
    "exact_completion_factorization",
    "local_puiseux_match",
    "eta_plus_one_puiseux_match",
    "measure_transport_identity",
    "monotone_flow_output_matching",
}

OPEN_BLOCKERS = [premise for premise in PREMISES if premise not in EXPORTED_PREMISES]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def has_all(subset: frozenset[str], *items: str) -> bool:
    return all(item in subset for item in items)


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


def currently_realized(minimal_sets: list[list[str]]) -> bool:
    return any(set(row).issubset(EXPORTED_PREMISES) for row in minimal_sets)


def build_payload() -> dict[str, Any]:
    reports = {name: load_json(path) for name, path in REPORTS.items()}

    cross_report_checks = {
        "input_report_count": len(reports),
        "completed_factorization_residual_pass": reports["completed_factorization"]["factorization_summary"]["residual_tolerance_pass"],
        "completed_factorization_current_live_kernel_is_strict": "K_strict_gate" in reports["completed_factorization"]["kernel_reading"]["current_live_kernel"],
        "diagrams_reactivation_all_evidence_found": reports["reactivation"]["source_cross_checks"]["diagrams_evidence_all_found"],
        "one_bit_frontier_name": reports["reactivation"]["one_bit_frontier"]["frontier_name"],
        "torsion_opinion_overall_verdict": reports["torsion_opinion"]["opinion_audit"]["overall_verdict"],
        "reynolds_annihilator_zero": reports["reynolds"]["matrix_certificate"]["reynolds_times_chi11_numerator_is_zero_matrix"],
        "reynolds_chi11_rank": reports["reynolds"]["matrix_certificate"]["chi11_numerator_rank"],
        "local_puiseux_match_pass": reports["puiseux_map"]["coefficient_match_passes_to_order_2"],
        "eta_plus_one_match_pass": reports["puiseux_eta_plus_one"]["coefficient_match_passes_through_eta_plus_one"],
        "measure_transport_balance_residual": reports["measure_transport"]["measure_balance"]["max_abs_balance_residual"],
        "monotone_flow_output_matching_residual": reports["inverse_branch_ode"]["integration_check"]["max_abs_output_residual"],
        "cutoff_formal_inverse_global_admissible": reports["global_admissibility"]["global_admissible_candidate"],
        "global_admissibility_obstructed": reports["global_admissibility"]["restricted_obstruction"]["is_obstructed"],
    }

    premise_status_rows = [
        {
            "premise": premise,
            "status": "EXPORTED_CERTIFICATE" if premise in EXPORTED_PREMISES else "OPEN_BLOCKER",
            "meaning": {
                "diagrams_legacy_carrier": "DIAGRAMS-backed legacy carrier reactivation evidence exists.",
                "exact_completion_factorization": "Pointwise strict-over-legacy completion factors reconstruct K_strict on d=1..11.",
                "local_puiseux_match": "Local two-channel phase-normalized Puiseux coefficients match through d^2.",
                "eta_plus_one_puiseux_match": "The local Puiseux candidate is extended through d^(eta+1).",
                "measure_transport_identity": "A differentiated output-matched branch yields exact unit measure transport bookkeeping.",
                "monotone_flow_output_matching": "The monotone inverse branch ODE reproduces output matching numerically.",
                "strict_transport_derivation": "Still missing: derive the transport/ODE/factors from strict nadsoliton dynamics rather than output matching.",
                "global_z12_map_derivation": "Still missing: prove a global admissible Z12 distance/transport map; one truncated formal inverse is obstructed.",
                "orientation_chi11_source": "Still missing: prove beta_tors/K_tors/topology exports the chi_11 unit-axis source.",
                "chi11_uniqueness": "Still missing: select chi_11 uniquely over chi_5 and chi_7.",
                "reynolds_obstruction_escape": "Still missing: show the source is not killed by full-Aut Reynolds averaging.",
                "role_transfer_theorem": "Still missing: authorize any legacy-role use on the strict side explicitly.",
                "qw2191_selector_discharge": "Still missing: actual strict-core selector closure/QW-2191 discharge.",
            }[premise],
        }
        for premise in PREMISES
    ]

    outcome_predicates: dict[str, Callable[[frozenset[str]], bool]] = {
        "finite_exact_completion_certificate": lambda s: has_all(s, "exact_completion_factorization"),
        "guarded_completed_kernel_candidate": lambda s: has_all(s, "diagrams_legacy_carrier", "exact_completion_factorization"),
        "local_two_channel_completion_candidate": lambda s: has_all(s, "local_puiseux_match"),
        "eta_plus_one_local_completion_candidate": lambda s: has_all(s, "local_puiseux_match", "eta_plus_one_puiseux_match"),
        "measure_transport_bookkeeping_candidate": lambda s: has_all(s, "measure_transport_identity"),
        "monotone_output_matching_flow_candidate": lambda s: has_all(s, "monotone_flow_output_matching"),
        "theorem_level_completed_legacy_to_strict_bridge": lambda s: has_all(
            s,
            "diagrams_legacy_carrier",
            "exact_completion_factorization",
            "strict_transport_derivation",
            "global_z12_map_derivation",
            "orientation_chi11_source",
            "chi11_uniqueness",
            "reynolds_obstruction_escape",
            "role_transfer_theorem",
        ),
        "selector_closed_completed_kernel_toe_step": lambda s: has_all(
            s,
            "diagrams_legacy_carrier",
            "exact_completion_factorization",
            "strict_transport_derivation",
            "global_z12_map_derivation",
            "orientation_chi11_source",
            "chi11_uniqueness",
            "reynolds_obstruction_escape",
            "role_transfer_theorem",
            "qw2191_selector_discharge",
        ),
        "strict_full_aut_internal_chi11_source": lambda _s: False,
    }

    antichains = {name: enumerate_minimal(predicate) for name, predicate in outcome_predicates.items()}
    outcome_rows = []
    for name, minimal_sets in antichains.items():
        min_count = None if not minimal_sets else min(len(row) for row in minimal_sets)
        outcome_rows.append(
            {
                "outcome": name,
                "minimal_premise_sets": minimal_sets,
                "minimal_premise_count": min_count,
                "currently_realized_by_loaded_reports": currently_realized(minimal_sets),
                "missing_open_blockers_for_first_minimal_set": None
                if not minimal_sets
                else [premise for premise in minimal_sets[0] if premise in OPEN_BLOCKERS],
            }
        )

    bridge_minimal = antichains["theorem_level_completed_legacy_to_strict_bridge"][0]
    remaining_for_bridge = [premise for premise in bridge_minimal if premise in OPEN_BLOCKERS]

    return {
        "result_kind": "SCRATCH_COMPLETED_KERNEL_BLOCKER_DEPENDENCY_LATTICE__NO_FALSE_PASS",
        "status": "minimal-premise-lattice-computed-for-completed-strict-kernel-bridge-frontier",
        "repo_grep_nonduplication_note": "Before adding this probe, rg was used over bridge/K_legacy_ont/K_strict_gate/beta_tors/chi_11/Reynolds/Puiseux/compression blockers to avoid repeating prior orbit, factorization, and local-transport audits.",
        "finite_lattice": {
            "premise_count": len(PREMISES),
            "premise_subset_count": 2 ** len(PREMISES),
            "exported_certificate_count": len(EXPORTED_PREMISES),
            "open_blocker_count": len(OPEN_BLOCKERS),
        },
        "premises": PREMISES,
        "exported_premises": sorted(EXPORTED_PREMISES, key=PREMISES.index),
        "open_blockers": OPEN_BLOCKERS,
        "cross_report_checks": cross_report_checks,
        "premise_status_rows": premise_status_rows,
        "minimal_premise_antichains": antichains,
        "outcome_rows": outcome_rows,
        "current_frontier_summary": {
            "current_live_kernel": "K_strict_gate is the live full kernel; legacy is the historical/nadsoliton-characteristic carrier.",
            "already_certified": sorted(EXPORTED_PREMISES, key=PREMISES.index),
            "remaining_for_theorem_level_bridge": remaining_for_bridge,
            "main_one_bit_blocker": "orientation_chi11_source plus chi11_uniqueness and reynolds_obstruction_escape",
            "transport_blocker": "strict_transport_derivation plus global_z12_map_derivation",
            "role_and_selector_blocker": "role_transfer_theorem remains below any QW-2191 selector discharge",
        },
        "exact_proof_certificate": {
            "finite_domain": f"The dependency lattice enumerates all 2^{len(PREMISES)}={2 ** len(PREMISES)} subsets of named bridge premises.",
            "separation_result": "Exact completion factorization, local Puiseux matching, measure transport bookkeeping, and monotone output-flow candidates are already certificates, but they are not strict derivations.",
            "bridge_antichain": "The theorem-level completed legacy-to-strict bridge requires strict transport derivation, global Z12 map derivation, chi_11 source/uniqueness/Reynolds escape, and role-transfer theorem in addition to the existing carrier/factorization certificates.",
            "full_aut_limit": "The strict_full_aut_internal_chi11_source outcome has no minimal premise set in this lattice because current Reynolds evidence keeps full-Aut invariant data in the annihilating/orthogonal sector.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in a solitonic state; legacy kernel data is a historical internal characteristic carrier completed by the strict kernel candidate lane.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced.",
        },
        "hard_limits": [
            "No claim that the legacy kernel is the current live kernel; K_strict_gate is the current full form.",
            "No unqualified identity K_legacy_ont == K_strict_gate is asserted.",
            "No beta_tors -> chi_11 theorem is asserted.",
            "No strict derivation of the completion factors, transport ODE, or global Z12 map is claimed.",
            "No legacy physical-role transfer onto K_strict_gate is used without an explicit bridge theorem.",
            "No QW-2191 discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    lattice = payload["finite_lattice"]
    frontier = payload["current_frontier_summary"]
    lines = [
        "# Completed kernel blocker dependency lattice",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite lattice",
        "",
        f"- Premises: `{lattice['premise_count']}`",
        f"- Enumerated subsets: `{lattice['premise_subset_count']}`",
        f"- Exported certificates: `{lattice['exported_certificate_count']}`",
        f"- Open blockers: `{lattice['open_blocker_count']}`",
        "",
        "## Current frontier summary",
        "",
    ]
    for key, value in frontier.items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Outcome antichains", ""])
    for row in payload["outcome_rows"]:
        lines.extend(
            [
                f"### {row['outcome']}",
                f"- Minimal premise count: `{row['minimal_premise_count']}`",
                f"- Currently realized: `{row['currently_realized_by_loaded_reports']}`",
                f"- Missing blockers for first set: `{row['missing_open_blockers_for_first_minimal_set']}`",
                f"- Minimal sets: `{row['minimal_premise_sets']}`",
                "",
            ]
        )
    lines.extend(["## Proof certificate", ""])
    for key, value in payload["exact_proof_certificate"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Hard limits", ""])
    lines.extend(f"- {item}" for item in payload["hard_limits"])
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    OUT_MD.write_text(write_markdown(payload), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
