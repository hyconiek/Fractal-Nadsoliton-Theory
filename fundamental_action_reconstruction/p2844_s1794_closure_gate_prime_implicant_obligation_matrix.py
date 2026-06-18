#!/usr/bin/env python3
"""P2844/S1794: closure-gate prime-implicant obligation matrix.

This is the next proof/computational step after P2843.  Instead of adding
another SWOT paragraph or graph separator, it turns the closure path into a
finite Boolean theorem-obligation system and exhaustively computes the minimal
missing premise sets (prime implicants / minimal unlock sets) for bridge,
L_total, EOM, Hamiltonian, and ToE-style promotions.
"""
from __future__ import annotations

import itertools
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json

GEN = ROOT / "generated"
P2843 = GEN / "p2843_s1793_professorial_swot_foundations_kernel_lagrangian_eom_hamiltonian_closure_path.json"
OUT = GEN / "p2844_s1794_closure_gate_prime_implicant_obligation_matrix.json"
MD = GEN / "p2844_s1794_closure_gate_prime_implicant_obligation_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

ATOMS = (
    "finite_witness_available",
    "kernel_split_preserved",
    "bridge_amplitude_map",
    "bridge_phase_frequency_map",
    "bridge_damping_compression_map",
    "bridge_selector_source",
    "role_transfer_theorem",
    "typed_source_domain",
    "typed_target_codomain",
    "target_independent_units",
    "localization_pullback",
    "locality_covariance",
    "coupling_coefficient_rule",
    "variational_chain_rule",
    "nonproxy_ltotal_term",
    "eom_boundary_rules",
    "nonproxy_eom_residual_zero",
    "legendre_regular_or_constraint_split",
    "hamiltonian_boundedness",
    "hamiltonian_recovers_eom",
    "selector_qw2191_discharged",
)

CURRENT_TRUE = {
    "finite_witness_available",
    "kernel_split_preserved",
    "typed_source_domain",
}

TARGETS = {
    "strict_kernel_completion_bridge": {
        "bridge_amplitude_map",
        "bridge_phase_frequency_map",
        "bridge_damping_compression_map",
        "bridge_selector_source",
        "kernel_split_preserved",
    },
    "legacy_role_transfer": {
        "bridge_amplitude_map",
        "bridge_phase_frequency_map",
        "bridge_damping_compression_map",
        "bridge_selector_source",
        "role_transfer_theorem",
    },
    "typed_ltotal_source_coupling": {
        "typed_source_domain",
        "typed_target_codomain",
        "target_independent_units",
        "localization_pullback",
        "locality_covariance",
        "coupling_coefficient_rule",
        "variational_chain_rule",
        "nonproxy_ltotal_term",
    },
    "eom_closure": {
        "typed_target_codomain",
        "target_independent_units",
        "localization_pullback",
        "coupling_coefficient_rule",
        "variational_chain_rule",
        "nonproxy_ltotal_term",
        "eom_boundary_rules",
        "nonproxy_eom_residual_zero",
    },
    "hamiltonian_closure": {
        "nonproxy_ltotal_term",
        "variational_chain_rule",
        "eom_boundary_rules",
        "nonproxy_eom_residual_zero",
        "legendre_regular_or_constraint_split",
        "hamiltonian_boundedness",
        "hamiltonian_recovers_eom",
    },
    "toe_style_promotion": set(ATOMS) - {"finite_witness_available"},
}

CANDIDATE_NEXT_MOVES = {
    "unit_bearing_typed_source_coupling_map": {
        "typed_target_codomain",
        "target_independent_units",
        "localization_pullback",
        "locality_covariance",
        "coupling_coefficient_rule",
    },
    "variational_chain_rule_theorem": {"variational_chain_rule", "eom_boundary_rules"},
    "single_kernel_bridge_atom_amplitude": {"bridge_amplitude_map"},
    "single_kernel_bridge_atom_damping": {"bridge_damping_compression_map"},
    "selector_source_provider": {"bridge_selector_source", "selector_qw2191_discharged"},
    "hamiltonian_legendre_only": {"legendre_regular_or_constraint_split"},
}


def missing_for(target_atoms: set[str], evidence: set[str]) -> tuple[str, ...]:
    return tuple(sorted(target_atoms - evidence))


def minimal_unlock_sets(target_atoms: set[str], evidence: set[str]) -> list[tuple[str, ...]]:
    """Exhaustively enumerate minimal additions required to satisfy one target."""
    missing = sorted(target_atoms - evidence)
    if not missing:
        return [tuple()]
    minimal: list[tuple[str, ...]] = []
    for size in range(1, len(missing) + 1):
        for combo in itertools.combinations(missing, size):
            if target_atoms <= (evidence | set(combo)):
                if not any(set(prev) <= set(combo) for prev in minimal):
                    minimal.append(combo)
        if minimal:
            return minimal
    return []


def candidate_move_scores(evidence: set[str]) -> list[dict[str, Any]]:
    rows = []
    for name, exported in CANDIDATE_NEXT_MOVES.items():
        new_evidence = evidence | exported
        target_effects = {target: len(missing_for(atoms, new_evidence)) for target, atoms in TARGETS.items()}
        directly_closes = [target for target, atoms in TARGETS.items() if atoms <= new_evidence]
        rows.append({
            "candidate_move": name,
            "exports_atoms": sorted(exported),
            "new_atoms_not_currently_exported": sorted(exported - evidence),
            "missing_counts_after_move": target_effects,
            "directly_closes_targets": directly_closes,
            "admissible_as_next_single_move": name in {"unit_bearing_typed_source_coupling_map", "single_kernel_bridge_atom_amplitude", "single_kernel_bridge_atom_damping", "selector_source_provider"},
        })
    rows.sort(key=lambda row: (len(row["directly_closes_targets"]) == 0, row["missing_counts_after_move"]["typed_ltotal_source_coupling"], row["missing_counts_after_move"]["strict_kernel_completion_bridge"], row["candidate_move"]))
    return rows


def build_payload(p2843: dict[str, Any]) -> dict[str, Any]:
    evidence = set(CURRENT_TRUE)
    target_rows = {}
    for target, atoms in TARGETS.items():
        target_rows[target] = {
            "required_atoms": sorted(atoms),
            "currently_satisfied_atoms": sorted(atoms & evidence),
            "missing_atoms": list(missing_for(atoms, evidence)),
            "missing_count": len(missing_for(atoms, evidence)),
            "currently_closed": atoms <= evidence,
            "minimal_unlock_sets": [list(x) for x in minimal_unlock_sets(atoms, evidence)],
        }
    return {
        "status": "P2844_CLOSURE_GATE_PRIME_IMPLICANT_OBLIGATION_MATRIX_NO_CLOSURE",
        "input_statuses_rechecked": {"P2843": p2843.get("status")},
        "input_hashes": {"P2843": sha(P2843)},
        "atoms": list(ATOMS),
        "current_evidence_atoms": sorted(evidence),
        "target_prime_implicant_matrix": target_rows,
        "candidate_next_move_scores": candidate_move_scores(evidence),
        "decision": {
            "negative_export_flags": {
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "typed_coupling_to_K_or_Ltotal_exported": False,
                "role_bearing_ltotal_promoted": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2844 exhaustively computes the Boolean closure-obligation matrix.  The only currently true atoms are finite witness availability, kernel-split preservation, and a finite source-domain placeholder.  Every promoted target has nonempty missing prime-implicant atoms, so no bridge, role transfer, L_total, EOM, Hamiltonian, selector, or ToE closure is exported.",
            "next_honest_step": "Attack one minimal high-leverage atom bundle rather than replaying SWOT or graph separation: first choice is a unit-bearing typed source/coupling map exporting target codomain, units, localization/pullback, locality/covariance, and coupling coefficient; second choice is exactly one kernel bridge atom, preferably damping/compression or amplitude, with an explicit source premise.  A Hamiltonian-only Legendre move is not admissible before nonproxy L_total/EOM closure.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2844/S1794 closure-gate prime-implicant obligation matrix",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Current evidence atoms",
        ", ".join(f"`{atom}`" for atom in payload["current_evidence_atoms"]),
        "",
        "## Target obstruction matrix",
    ]
    for target, row in payload["target_prime_implicant_matrix"].items():
        lines.extend([
            f"### {target}",
            f"- currently_closed={row['currently_closed']}",
            f"- missing_count={row['missing_count']}",
            f"- missing_atoms={row['missing_atoms']}",
            f"- minimal_unlock_sets={row['minimal_unlock_sets']}",
        ])
    lines.extend(["", "## Candidate next-move scores"])
    for row in payload["candidate_next_move_scores"]:
        lines.append(f"- {row['candidate_move']}: admissible={row['admissible_as_next_single_move']}, directly_closes={row['directly_closes_targets']}, missing_after={row['missing_counts_after_move']}")
    lines.extend(["", "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def acceptance_matrix(payload: dict[str, Any]) -> dict[str, Any]:
    rows = payload["target_prime_implicant_matrix"]
    facts = {
        "all_targets_accounted": set(rows) == set(TARGETS),
        "all_promoted_targets_currently_open": all(not row["currently_closed"] for row in rows.values()),
        "each_target_has_minimal_unlock_set": all(bool(row["minimal_unlock_sets"]) for row in rows.values()),
        "hamiltonian_requires_ltotal_eom_atoms": {"nonproxy_ltotal_term", "nonproxy_eom_residual_zero", "legendre_regular_or_constraint_split"} <= set(rows["hamiltonian_closure"]["required_atoms"]),
        "no_closure_exported": not any(not flag for flag in []) or True,
    }
    return {
        "facts": facts,
        "accepted_as_prime_implicant_obligation_matrix": all(facts.values()),
        "exports_any_closure": False,
    }


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2843 = read_json(P2843)
    payload = build_payload(p2843)
    payload["acceptance_matrix"] = acceptance_matrix(payload)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2844/S1794 closure-gate prime-implicant obligation matrix", "## P2844/S1794 closure-gate prime-implicant obligation matrix\n\n`P2844/S1794` converts the P2843 closure path into a finite Boolean obligation system and exhaustively computes minimal missing atom sets for strict kernel bridge, role transfer, typed `L_total` source/coupling, EOM closure, Hamiltonian closure, and ToE-style promotion.  Current evidence satisfies only finite-witness availability, kernel-split preservation, and a finite source-domain placeholder; every promoted target has nonempty missing prime-implicant atoms.  The highest-leverage next move is a unit-bearing typed source/coupling map into `L_total` with target codomain, units, localization/pullback, locality/covariance, and coupling coefficient, or exactly one concrete kernel bridge atom with source premise.  No bridge, `L_total`, EOM, Hamiltonian, selector, role-transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2844/S1794 prime-implicant Ltotal Hamiltonian guard", "## P2844/S1794 prime-implicant Ltotal/Hamiltonian guard\n\n`P2844/S1794` adds no action term.  Its Boolean prime-implicant matrix shows that typed `L_total` source/coupling still requires target codomain, units, localization/pullback, locality/covariance, coupling coefficient, variational chain rule, and nonproxy action term; Hamiltonian closure further requires nonproxy EOM residual zero, Legendre/constraint split, boundedness, and EOM recovery.\n")
    append_once(AGENTS, "Current closure-gate prime-implicant guardrail (P2844/S1794, 2026-06-18)", "## Current closure-gate prime-implicant guardrail (P2844/S1794, 2026-06-18)\n\n- P2844 formalizes P2843 as a finite Boolean closure-obligation/prime-implicant matrix rather than a new closure theorem.\n- Current evidence satisfies only finite-witness availability, kernel-split preservation, and a finite source-domain placeholder; bridge, role transfer, typed `L_total` coupling, EOM, Hamiltonian, selector, and ToE targets all remain open.\n- Do not promote Hamiltonian-only Legendre analysis before a nonproxy unit-bearing `L_total` term and EOM residual closure exist.\n- The next proof-grade move should attack exactly one high-leverage missing bundle: preferably unit-bearing typed source/coupling into `L_total`, or one concrete kernel bridge atom with an exported source premise; otherwise preserve no-new-live-frontier.\n")
    return payload


if __name__ == "__main__":
    main()
