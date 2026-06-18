#!/usr/bin/env python3
"""P2852/S1802: kernel bridge obligation reconciliation matrix.

P2851 recommends a combined bridge-obligation reconciliation over the bridge
atoms that have now been tested separately: damping/compression (P2849), EML
syntax impact (P2850), and amplitude normalization (P2851).  P2852 is a finite
Boolean proof ledger: it recomputes which obligations are discharged, which are
blocked, and which single typed source atom would be an honest next input rather
than a replay of a closed lane.
"""
from __future__ import annotations

import json
from itertools import combinations
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json

GEN = ROOT / "generated"
P2849 = GEN / "p2849_s1799_damping_compression_kernel_bridge_atom_audit.json"
P2850 = GEN / "p2850_s1800_eml_single_operator_kernel_bridge_impact_audit.json"
P2851 = GEN / "p2851_s1801_amplitude_normalization_kernel_bridge_atom_audit.json"
OUT = GEN / "p2852_s1802_kernel_bridge_obligation_reconciliation_matrix.json"
MD = GEN / "p2852_s1802_kernel_bridge_obligation_reconciliation_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

ATOMS = (
    "legacy_kernel_defined",
    "strict_kernel_defined",
    "eml_syntax_basis_available",
    "damping_compression_atom_exported",
    "amplitude_normalization_atom_exported",
    "phase_frequency_bridge_atom_exported",
    "selector_topological_source_exported",
    "eta_source_law_exported",
    "target_independent_beta_source_exported",
    "alpha_geo_strict_source_law_exported",
    "completion_map_semantics_exported",
    "role_transfer_theorem_exported",
)

TARGETS = {
    "damping_compression_bridge_atom": {
        "damping_compression_atom_exported",
        "eta_source_law_exported",
        "target_independent_beta_source_exported",
    },
    "amplitude_normalization_bridge_atom": {
        "amplitude_normalization_atom_exported",
        "alpha_geo_strict_source_law_exported",
        "phase_frequency_bridge_atom_exported",
        "damping_compression_atom_exported",
    },
    "syntax_level_common_expression_basis": {
        "legacy_kernel_defined",
        "strict_kernel_defined",
        "eml_syntax_basis_available",
    },
    "full_kernel_completion_bridge": {
        "damping_compression_atom_exported",
        "amplitude_normalization_atom_exported",
        "phase_frequency_bridge_atom_exported",
        "selector_topological_source_exported",
        "eta_source_law_exported",
        "target_independent_beta_source_exported",
        "alpha_geo_strict_source_law_exported",
        "completion_map_semantics_exported",
    },
    "role_transfer": {
        "full_kernel_completion_bridge",
        "role_transfer_theorem_exported",
    },
}

CANDIDATE_NEXT_ATOMS = {
    "phase_frequency_bridge_atom": {
        "adds": {"phase_frequency_bridge_atom_exported"},
        "replay_risk": "low",
        "reason": "Not yet audited as its own bridge atom; required by amplitude passage and full completion.",
    },
    "eta_beta_strict_source_law": {
        "adds": {"eta_source_law_exported", "target_independent_beta_source_exported"},
        "replay_risk": "medium",
        "reason": "P2849 blocked legacy-linear sourcing, but a genuinely new strict source law would be admissible.",
    },
    "alpha_geo_strict_source_law": {
        "adds": {"alpha_geo_strict_source_law_exported"},
        "replay_risk": "medium",
        "reason": "P2851 blocked constant rescaling, but a new role-safe strict alpha source would be admissible.",
    },
    "eml_syntax_replay": {
        "adds": {"eml_syntax_basis_available"},
        "replay_risk": "high",
        "reason": "P2850 already records syntax availability; replay does not add typed source semantics.",
    },
}


def current_atoms(p2849: dict[str, Any], p2850: dict[str, Any], p2851: dict[str, Any]) -> dict[str, bool]:
    return {
        "legacy_kernel_defined": True,
        "strict_kernel_defined": True,
        "eml_syntax_basis_available": p2850.get("acceptance_matrix", {}).get("accepted_as_eml_external_information_impact_audit", False),
        "damping_compression_atom_exported": p2849.get("acceptance_matrix", {}).get("exports_damping_compression_bridge_atom", False),
        "amplitude_normalization_atom_exported": p2851.get("acceptance_matrix", {}).get("exports_amplitude_normalization_bridge_atom", False),
        "phase_frequency_bridge_atom_exported": False,
        "selector_topological_source_exported": False,
        "eta_source_law_exported": False,
        "target_independent_beta_source_exported": False,
        "alpha_geo_strict_source_law_exported": False,
        "completion_map_semantics_exported": False,
        "role_transfer_theorem_exported": False,
    }


def target_statuses(atoms: dict[str, bool]) -> dict[str, Any]:
    statuses: dict[str, Any] = {}
    for target, requirements in TARGETS.items():
        expanded_requirements = set(requirements)
        if target == "role_transfer":
            full_bridge_status = statuses.get("full_kernel_completion_bridge", {})
            full_bridge_ok = bool(full_bridge_status.get("satisfied", False))
            satisfied = full_bridge_ok and atoms["role_transfer_theorem_exported"]
            missing = []
            if not full_bridge_ok:
                missing.append("full_kernel_completion_bridge")
            if not atoms["role_transfer_theorem_exported"]:
                missing.append("role_transfer_theorem_exported")
        else:
            satisfied = all(atoms[atom] for atom in expanded_requirements)
            missing = sorted(atom for atom in expanded_requirements if not atoms[atom])
        statuses[target] = {"satisfied": satisfied, "missing": missing}
    return statuses


def minimal_unlock_sets(atoms: dict[str, bool], target: str, max_size: int = 3) -> list[list[str]]:
    if target == "role_transfer":
        requirements = set(TARGETS["full_kernel_completion_bridge"]) | {"role_transfer_theorem_exported"}
    else:
        requirements = set(TARGETS[target])
    missing = sorted(atom for atom in requirements if not atoms.get(atom, False))
    unlocks: list[list[str]] = []
    for size in range(1, min(max_size, len(missing)) + 1):
        for combo in combinations(missing, size):
            # The target unlocks exactly when the combo covers all missing requirements.
            if set(combo) == set(missing):
                unlocks.append(list(combo))
        if unlocks:
            break
    return unlocks


def score_candidates(atoms: dict[str, bool]) -> dict[str, Any]:
    base = target_statuses(atoms)
    base_missing_total = sum(len(row["missing"]) for row in base.values())
    rows = {}
    for name, candidate in CANDIDATE_NEXT_ATOMS.items():
        new_atoms = dict(atoms)
        for atom in candidate["adds"]:
            new_atoms[atom] = True
        new_status = target_statuses(new_atoms)
        new_missing_total = sum(len(row["missing"]) for row in new_status.values())
        rows[name] = {
            "adds": sorted(candidate["adds"]),
            "replay_risk": candidate["replay_risk"],
            "reason": candidate["reason"],
            "missing_obligation_reduction": base_missing_total - new_missing_total,
            "directly_unlocks_targets": sorted(target for target, row in new_status.items() if row["satisfied"] and not base[target]["satisfied"]),
            "admissible_without_new_source_premise": False if candidate["replay_risk"] != "high" else False,
        }
    return rows


def build_payload(p2849: dict[str, Any], p2850: dict[str, Any], p2851: dict[str, Any]) -> dict[str, Any]:
    atoms = current_atoms(p2849, p2850, p2851)
    statuses = target_statuses(atoms)
    unlocks = {target: minimal_unlock_sets(atoms, target) for target in TARGETS}
    candidates = score_candidates(atoms)
    facts = {
        "p2849_rechecked": p2849.get("status") == "P2849_DAMPING_COMPRESSION_KERNEL_BRIDGE_ATOM_AUDIT_NO_CLOSURE",
        "p2850_rechecked": p2850.get("status") == "P2850_EML_SINGLE_OPERATOR_KERNEL_BRIDGE_IMPACT_AUDIT_NO_CLOSURE",
        "p2851_rechecked": p2851.get("status") == "P2851_AMPLITUDE_NORMALIZATION_KERNEL_BRIDGE_ATOM_AUDIT_NO_CLOSURE",
        "syntax_basis_only_target_satisfied": statuses["syntax_level_common_expression_basis"]["satisfied"],
        "no_semantic_bridge_target_satisfied": not any(statuses[target]["satisfied"] for target in ("damping_compression_bridge_atom", "amplitude_normalization_bridge_atom", "full_kernel_completion_bridge", "role_transfer")),
        "full_bridge_has_missing_atoms": len(statuses["full_kernel_completion_bridge"]["missing"]) > 0,
        "role_transfer_still_downstream": "full_kernel_completion_bridge" in statuses["role_transfer"]["missing"],
    }
    return {
        "status": "P2852_KERNEL_BRIDGE_OBLIGATION_RECONCILIATION_MATRIX_NO_CLOSURE",
        "input_hashes": {"P2849": sha(P2849), "P2850": sha(P2850), "P2851": sha(P2851)},
        "kernel_bridge_obligation_reconciliation": {
            "input_statuses_rechecked": {"P2849": p2849.get("status"), "P2850": p2850.get("status"), "P2851": p2851.get("status")},
            "atoms": atoms,
            "target_requirements": {key: sorted(value) for key, value in TARGETS.items()},
            "target_statuses": statuses,
            "minimal_unlock_sets_size_leq_3": unlocks,
            "candidate_next_atom_scores": candidates,
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_bridge_reconciliation_obstruction_audit": all(facts.values()),
            "exports_full_kernel_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "damping_compression_bridge_atom_exported": False,
                "amplitude_normalization_bridge_atom_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "selector_closure_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2852 reconciles P2849-P2851.  Only the syntax-level common expression basis is satisfied; semantic bridge targets remain unsatisfied.  Damping/compression still lacks eta and target-independent beta sources, amplitude still lacks a constant residual-zero passage and alpha source law, and full bridge/role transfer remain downstream of multiple missing typed atoms.",
            "next_honest_step": "Do not replay damping, amplitude, EML syntax, density normalizers, role transfer, L_total, EOM, Hamiltonian, or ToE promotion.  The next admissible proof-grade move must introduce exactly one new typed bridge-source atom.  The cleanest non-replay candidate is a phase/frequency bridge-source audit for omega/phi/topological data; alternatively introduce a genuinely new strict eta/beta source law.  Without such a new source premise, preserve the no-new-live-frontier certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    rec = payload["kernel_bridge_obligation_reconciliation"]
    lines = [
        "# P2852/S1802 kernel bridge obligation reconciliation matrix",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Target statuses",
    ]
    for target, row in rec["target_statuses"].items():
        lines.append(f"- {target}: satisfied={row['satisfied']}; missing={row['missing']}")
    lines.extend(["", "## Candidate next atom scores"])
    for name, row in rec["candidate_next_atom_scores"].items():
        lines.append(f"- {name}: adds={row['adds']}; reduction={row['missing_obligation_reduction']}; replay_risk={row['replay_risk']}; unlocks={row['directly_unlocks_targets']}")
    lines.extend([
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload(read_json(P2849), read_json(P2850), read_json(P2851))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2852/S1802 kernel bridge obligation reconciliation matrix", "## P2852/S1802 kernel bridge obligation reconciliation matrix\n\n`P2852/S1802` reconciles P2849-P2851 in a finite Boolean bridge-obligation matrix.  Only the syntax-level common expression basis is satisfied; semantic bridge targets remain open.  Damping/compression lacks `eta` and target-independent `beta` sources, amplitude lacks constant residual-zero passage and role-safe `alpha_geo` source law, and full bridge/role transfer remain downstream of multiple missing typed atoms.  No damping bridge, amplitude bridge, full kernel bridge, role transfer, `L_total`, EOM, Hamiltonian, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2852/S1802 bridge reconciliation Ltotal guard", "## P2852/S1802 bridge reconciliation `L_total` guard\n\n`P2852/S1802` adds no action term.  Its bridge-obligation matrix confirms that current bridge atoms do not provide a unit-bearing source density, nonproxy `L_total` insertion, EOM, Hamiltonian, or role-bearing term.\n")
    append_once(AGENTS, "Current kernel bridge obligation reconciliation guardrail (P2852/S1802, 2026-06-18)", "## Current kernel bridge obligation reconciliation guardrail (P2852/S1802, 2026-06-18)\n\n- P2852 reconciles P2849-P2851: only syntax-level common expression basis is satisfied; semantic bridge targets remain unsatisfied.\n- Damping/compression lacks `eta` and target-independent `beta` sources; amplitude lacks constant residual-zero passage and role-safe `alpha_geo` source law; full bridge and role transfer remain downstream.\n- Do not replay damping, amplitude, EML syntax, density normalizers, role transfer, `L_total`, EOM, Hamiltonian, or ToE promotion as closure evidence.\n- The next admissible proof-grade move must introduce exactly one new typed bridge-source atom, preferably a phase/frequency bridge-source audit for `omega/phi`/topological data or a genuinely new strict `eta/beta` source law; otherwise preserve no-new-live-frontier.\n")
    return payload


if __name__ == "__main__":
    main()
