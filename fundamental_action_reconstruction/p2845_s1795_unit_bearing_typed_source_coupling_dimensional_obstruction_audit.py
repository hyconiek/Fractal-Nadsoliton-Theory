#!/usr/bin/env python3
"""P2845/S1795: unit-bearing typed source/coupling dimensional obstruction audit.

P2844 says the highest-leverage next move is a unit-bearing typed source/coupling
map into L_total.  P2845 attacks that exact bundle with a finite dimensional
obligation matrix: candidate source terms are dimension-checked, then audited for
typed codomain, localization/pullback, covariance, coefficient source, and
variational-chain-rule premises.  Formal unit balance alone is not promoted to a
source law.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json

GEN = ROOT / "generated"
P2844 = GEN / "p2844_s1794_closure_gate_prime_implicant_obligation_matrix.json"
OUT = GEN / "p2845_s1795_unit_bearing_typed_source_coupling_dimensional_obstruction_audit.json"
MD = GEN / "p2845_s1795_unit_bearing_typed_source_coupling_dimensional_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

SPACETIME_DIMENSION = 4
ACTION_DIMENSION = 0
LAGRANGIAN_DENSITY_DIMENSION = 4
FINITE_GRAPH_FUNCTIONAL_DIMENSION = 0

REQUIRED_PREMISES = (
    "typed_source_domain",
    "typed_target_codomain",
    "dimension_balanced_units",
    "target_independent_units",
    "localization_pullback",
    "locality_covariance",
    "coupling_coefficient_rule",
    "variational_chain_rule",
    "nonproxy_ltotal_term",
)


@dataclass(frozen=True)
class CandidateTerm:
    name: str
    schematic_term: str
    operator_dimension: Fraction
    graph_factor_dimension: Fraction
    has_typed_target_codomain: bool
    has_localization_pullback: bool
    has_locality_covariance: bool
    has_coupling_coefficient_rule: bool
    has_variational_chain_rule: bool
    has_nonproxy_ltotal_term: bool
    notes: str

    @property
    def required_coupling_dimension(self) -> Fraction:
        return Fraction(LAGRANGIAN_DENSITY_DIMENSION) - self.operator_dimension - self.graph_factor_dimension

    @property
    def dimension_balanced_units(self) -> bool:
        # In natural units, a formal coupling with this mass dimension can always
        # be assigned.  This is only a formal dimensional balance, not a source.
        return True

    @property
    def target_independent_units(self) -> bool:
        # P2836/P2844 leave the unit/source of the coupling unexported.
        return False

    def premise_payload(self) -> dict[str, Any]:
        premises = {
            "typed_source_domain": True,
            "typed_target_codomain": self.has_typed_target_codomain,
            "dimension_balanced_units": self.dimension_balanced_units,
            "target_independent_units": self.target_independent_units,
            "localization_pullback": self.has_localization_pullback,
            "locality_covariance": self.has_locality_covariance,
            "coupling_coefficient_rule": self.has_coupling_coefficient_rule,
            "variational_chain_rule": self.has_variational_chain_rule,
            "nonproxy_ltotal_term": self.has_nonproxy_ltotal_term,
        }
        return {
            "candidate": self.name,
            "schematic_term": self.schematic_term,
            "operator_dimension": frac_payload(self.operator_dimension),
            "graph_factor_dimension": frac_payload(self.graph_factor_dimension),
            "required_coupling_dimension": frac_payload(self.required_coupling_dimension),
            "premises": premises,
            "missing_premises": [key for key in REQUIRED_PREMISES if not premises[key]],
            "accepted_as_unit_bearing_typed_source_coupling": all(premises.values()),
            "notes": self.notes,
        }


def frac_payload(value: Fraction) -> dict[str, int | str]:
    return {"numerator": value.numerator, "denominator": value.denominator, "display": str(value)}


def candidate_terms() -> tuple[CandidateTerm, ...]:
    return (
        CandidateTerm(
            name="dimensionless_graph_times_dimension4_scalar_density",
            schematic_term="Delta L = lambda_0 * F_G * O_4(x)",
            operator_dimension=Fraction(4),
            graph_factor_dimension=Fraction(0),
            has_typed_target_codomain=True,
            has_localization_pullback=False,
            has_locality_covariance=True,
            has_coupling_coefficient_rule=False,
            has_variational_chain_rule=False,
            has_nonproxy_ltotal_term=False,
            notes="Formal dimensions allow lambda_0 dimension 0, but no strict source selects O_4(x), lambda_0, or graph-to-field pullback.",
        ),
        CandidateTerm(
            name="graph_mass_term_scalar_square",
            schematic_term="Delta L = lambda_2 * F_G * phi(x)^2",
            operator_dimension=Fraction(2),
            graph_factor_dimension=Fraction(0),
            has_typed_target_codomain=True,
            has_localization_pullback=False,
            has_locality_covariance=True,
            has_coupling_coefficient_rule=False,
            has_variational_chain_rule=False,
            has_nonproxy_ltotal_term=False,
            notes="Formal dimensions require lambda_2 with mass dimension 2; no target-independent mass/unit source is exported.",
        ),
        CandidateTerm(
            name="strict_kernel_modulated_density",
            schematic_term="Delta L = lambda_K * F_G * K_strict_gate(d) * O_4(x)",
            operator_dimension=Fraction(4),
            graph_factor_dimension=Fraction(0),
            has_typed_target_codomain=True,
            has_localization_pullback=False,
            has_locality_covariance=False,
            has_coupling_coefficient_rule=False,
            has_variational_chain_rule=False,
            has_nonproxy_ltotal_term=False,
            notes="K_strict_gate is dimensionless as written, but no strict d(x,y)/graph pullback or source coefficient is exported.",
        ),
        CandidateTerm(
            name="localized_delta_source",
            schematic_term="Delta L = lambda_delta * F_G * delta^4(x-x_G)",
            operator_dimension=Fraction(4),
            graph_factor_dimension=Fraction(0),
            has_typed_target_codomain=True,
            has_localization_pullback=False,
            has_locality_covariance=False,
            has_coupling_coefficient_rule=False,
            has_variational_chain_rule=False,
            has_nonproxy_ltotal_term=False,
            notes="A delta density can be dimension-balanced, but x_G is exactly the missing strict localization object.",
        ),
        CandidateTerm(
            name="hamiltonian_potential_placeholder",
            schematic_term="Delta H = c_H * F_G * q^2",
            operator_dimension=Fraction(2),
            graph_factor_dimension=Fraction(0),
            has_typed_target_codomain=False,
            has_localization_pullback=False,
            has_locality_covariance=False,
            has_coupling_coefficient_rule=False,
            has_variational_chain_rule=False,
            has_nonproxy_ltotal_term=False,
            notes="Hamiltonian placeholder is premature: no phase-space variables, Legendre map, constraints, or L_total source term exist.",
        ),
    )


def build_audit(p2844: dict[str, Any]) -> dict[str, Any]:
    rows = [candidate.premise_payload() for candidate in candidate_terms()]
    accepted = [row for row in rows if row["accepted_as_unit_bearing_typed_source_coupling"]]
    blocker_histogram: dict[str, int] = {key: 0 for key in REQUIRED_PREMISES}
    for row in rows:
        for missing in row["missing_premises"]:
            blocker_histogram[missing] += 1
    return {
        "input_statuses_rechecked": {"P2844": p2844.get("status")},
        "dimension_conventions": {
            "spacetime_dimension": SPACETIME_DIMENSION,
            "action_dimension": ACTION_DIMENSION,
            "lagrangian_density_mass_dimension": LAGRANGIAN_DENSITY_DIMENSION,
            "finite_graph_functional_dimension": FINITE_GRAPH_FUNCTIONAL_DIMENSION,
        },
        "required_premises": list(REQUIRED_PREMISES),
        "candidate_rows": rows,
        "accepted_candidate_count": len(accepted),
        "accepted_candidates": [row["candidate"] for row in accepted],
        "blocker_histogram": blocker_histogram,
        "formal_dimension_balance_count": sum(1 for row in rows if row["premises"]["dimension_balanced_units"]),
        "target_independent_units_count": sum(1 for row in rows if row["premises"]["target_independent_units"]),
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2844_rechecked": audit["input_statuses_rechecked"]["P2844"] == "P2844_CLOSURE_GATE_PRIME_IMPLICANT_OBLIGATION_MATRIX_NO_CLOSURE",
        "all_candidates_dimension_checked": audit["formal_dimension_balance_count"] == len(audit["candidate_rows"]),
        "no_candidate_exports_target_independent_units": audit["target_independent_units_count"] == 0,
        "no_candidate_accepted_as_source_coupling": audit["accepted_candidate_count"] == 0,
        "localization_pullback_still_blocked": audit["blocker_histogram"]["localization_pullback"] == len(audit["candidate_rows"]),
        "coupling_rule_still_blocked": audit["blocker_histogram"]["coupling_coefficient_rule"] == len(audit["candidate_rows"]),
        "hamiltonian_not_promoted": True,
    }
    return {
        "facts": facts,
        "accepted_as_unit_bearing_coupling_obstruction_audit": all(facts.values()),
        "exports_unit_bearing_typed_source_coupling": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["unit_bearing_typed_source_coupling_audit"]
    lines = [
        "# P2845/S1795 unit-bearing typed source/coupling dimensional obstruction audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Dimension conventions",
        f"- spacetime_dimension={audit['dimension_conventions']['spacetime_dimension']}",
        f"- lagrangian_density_mass_dimension={audit['dimension_conventions']['lagrangian_density_mass_dimension']}",
        f"- finite_graph_functional_dimension={audit['dimension_conventions']['finite_graph_functional_dimension']}",
        "",
        "## Candidate rows",
    ]
    for row in audit["candidate_rows"]:
        lines.extend([
            f"### {row['candidate']}",
            f"- term={row['schematic_term']}",
            f"- required_coupling_dimension={row['required_coupling_dimension']['display']}",
            f"- accepted={row['accepted_as_unit_bearing_typed_source_coupling']}",
            f"- missing_premises={row['missing_premises']}",
            f"- notes={row['notes']}",
        ])
    lines.extend([
        "",
        "## Summary",
        f"- accepted_candidate_count={audit['accepted_candidate_count']}",
        f"- formal_dimension_balance_count={audit['formal_dimension_balance_count']}",
        f"- target_independent_units_count={audit['target_independent_units_count']}",
        f"- blocker_histogram={audit['blocker_histogram']}",
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
    p2844 = read_json(P2844)
    audit = build_audit(p2844)
    payload: dict[str, Any] = {
        "status": "P2845_UNIT_BEARING_TYPED_SOURCE_COUPLING_DIMENSIONAL_OBSTRUCTION_NO_CLOSURE",
        "input_hashes": {"P2844": sha(P2844)},
        "unit_bearing_typed_source_coupling_audit": audit,
        "decision": {
            "negative_export_flags": {
                "unit_bearing_typed_source_coupling_exported": False,
                "target_independent_units_exported": False,
                "localization_pullback_exported": False,
                "coupling_coefficient_rule_exported": False,
                "nonproxy_ltotal_term_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "toe_closure_exported": False,
            },
            "reason": "P2845 attacks the P2844 high-leverage unit-bearing typed source/coupling bundle.  All five candidate terms can be formally dimension-balanced by assigning a coupling mass dimension, but none exports target-independent units, localization/pullback, a strict coupling coefficient rule, a variational chain rule, or a nonproxy L_total term.  Formal dimensional bookkeeping is therefore not a source law.",
            "next_honest_step": "Do not replay dimensional ansätze or Hamiltonian placeholders.  The next admissible proof-grade move should isolate exactly one missing source premise from P2845: either a strict localization/pullback object x_G or rho_G(x), or a target-independent coupling coefficient/unit source for one named density.  If neither is supplied, pivot to one concrete kernel bridge atom with an exported source premise; otherwise preserve no-new-live-frontier.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2845/S1795 unit-bearing typed source/coupling dimensional obstruction audit", "## P2845/S1795 unit-bearing typed source/coupling dimensional obstruction audit\n\n`P2845/S1795` attacks the P2844 high-leverage typed `L_total` source/coupling bundle with exact dimension bookkeeping for five candidate terms: `F_G O_4(x)`, `F_G phi^2`, `F_G K_strict_gate(d) O_4(x)`, `F_G delta^4(x-x_G)`, and a Hamiltonian placeholder.  All candidates can be formally dimension-balanced by assigning a coupling mass dimension, but none exports target-independent units, localization/pullback, a strict coupling coefficient rule, a variational chain rule, or a nonproxy `L_total` term.  No unit-bearing typed source/coupling, EOM, Hamiltonian, bridge, role-transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2845/S1795 dimensional ansatz Ltotal guard", "## P2845/S1795 dimensional ansatz Ltotal guard\n\n`P2845/S1795` adds no action term.  Formal unit balance for candidate graph-weighted densities is insufficient: current artifacts still lack target-independent units, graph-to-field localization/pullback, a coupling coefficient source, variational chain rule, and nonproxy `L_total` insertion.\n")
    append_once(AGENTS, "Current unit-bearing source/coupling dimensional obstruction guardrail (P2845/S1795, 2026-06-18)", "## Current unit-bearing source/coupling dimensional obstruction guardrail (P2845/S1795, 2026-06-18)\n\n- P2845 attacks the P2844 high-leverage typed `L_total` source/coupling bundle with exact dimension bookkeeping.\n- Formal dimension balance is available for candidate graph-weighted densities, but no candidate exports target-independent units, localization/pullback, strict coupling coefficient rule, variational chain rule, or nonproxy `L_total` term.\n- Do not promote dimensional ansätze or Hamiltonian placeholders to unit-bearing typed source/coupling, EOM, Hamiltonian, bridge, role-transfer, or ToE closure.\n- A next admissible move must isolate one missing source premise: either a strict localization/pullback object `x_G`/`rho_G(x)` or a target-independent coupling coefficient/unit source for one named density; otherwise pivot to one concrete kernel bridge atom with an exported source premise or preserve no-new-live-frontier.\n")
    return payload


if __name__ == "__main__":
    main()
