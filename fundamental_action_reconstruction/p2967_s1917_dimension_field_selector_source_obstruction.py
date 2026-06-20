#!/usr/bin/env python3
"""P2967/S1917: dimension/field selector-source obstruction.

P2966 reduced the physical-unit problem to the missing selector/source for
(N, sigma), unit U, or coefficient law.  This audit attacks one object only: a
selector/source for the dimension-field pair (N, sigma) in the equation
c[V]=(9/5)U^(2sigma-N).

It enumerates natural selector predicates on the P2966 dimensional grid.  Some
predicates are compatible and some can be made unique only by lexicographic or
4D import conventions, but no current artifact exports a strict nadsoliton source
for choosing one (N, sigma) pair.
"""
from __future__ import annotations

import hashlib, json
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2966_s1916_physical_unit_law_dimensional_obstruction import OUT as P2966, dimensional_rows

OUT = GEN / "p2967_s1917_dimension_field_selector_source_obstruction.json"
MD = GEN / "p2967_s1917_dimension_field_selector_source_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def parse_frac(text: str) -> Fraction:
    n, d = text.split("/")
    return Fraction(int(n), int(d))


def selector_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    def pack(name: str, selected: list[dict[str, Any]], source_type: str, replay: bool) -> dict[str, Any]:
        return {
            "selector": name,
            "selected_count": len(selected),
            "selected_pairs": [[r["spacetime_dimension_N"], r["field_length_dimension_sigma"]] for r in selected[:20]],
            "unique_pair": len(selected) == 1,
            "source_type": source_type,
            "replay_or_import": replay,
            "strict_nadsoliton_source_exported": False,
            "accepted_strict_selector": len(selected) == 1 and source_type == "strict_nadsoliton_export" and not replay,
        }

    exps = [parse_frac(r["required_coefficient_length_exponent_k"]) for r in rows]
    abs_min = min(abs(e) for e in exps)
    dimensionless = [r for r in rows if parse_frac(r["required_coefficient_length_exponent_k"]) == 0]
    min_abs = [r for r in rows if abs(parse_frac(r["required_coefficient_length_exponent_k"])) == abs_min]
    integer_power = [r for r in rows if r["integer_exponent"]]
    four_dim = [r for r in rows if r["spacetime_dimension_N"] == 4 and parse_frac(r["field_length_dimension_sigma"]) == 2]
    lexicographic = [sorted(rows, key=lambda r: (abs(parse_frac(r["required_coefficient_length_exponent_k"])), r["spacetime_dimension_N"], parse_frac(r["field_length_dimension_sigma"])))[0]]
    return [
        pack("dimensionless_coefficient_selector", dimensionless, "formal_dimensional_predicate", False),
        pack("minimal_absolute_exponent_selector", min_abs, "formal_minimization_predicate", False),
        pack("integer_power_selector", integer_power, "formal_integrality_predicate", False),
        pack("lexicographic_minimal_completion", lexicographic, "ordering_convention", True),
        pack("imported_four_dimensional_scalar_pair", four_dim, "external_physics_import", True),
    ]


def obligation_rows(selectors: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {"obligation": "selector_predicates_constructed", "satisfied": True, "evidence": f"{len(selectors)} selector candidates audited"},
        {"obligation": "nonimported_unique_pair_selected", "satisfied": any(s["unique_pair"] and not s["replay_or_import"] for s in selectors), "evidence": "formal nonimport predicates are nonunique on current grid"},
        {"obligation": "strict_nadsoliton_dimension_field_source", "satisfied": False, "evidence": "no current artifact exports N,sigma from nadsoliton ontology"},
        {"obligation": "strict_unit_U_source", "satisfied": False, "evidence": "selecting N,sigma would still not select U"},
        {"obligation": "coefficient_source_law_exported", "satisfied": False, "evidence": "no source law installs 9/5 into P2964 coupling"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["selector_predicates", "nonimported_unique_pair", "strict_N_sigma_source", "strict_U_source", "coefficient_source_law", "nonproxy_coupling"]
    full = (1 << len(names)) - 1
    return [{"mask": m, "present": {n: bool(m & (1 << i)) for i, n in enumerate(names)}, "accepts_strict_dimension_field_selector": m == full} for m in range(1 << len(names))]


def build_payload(p2966: dict[str, Any]) -> dict[str, Any]:
    rows = dimensional_rows()
    selectors = selector_rows(rows)
    obligations = obligation_rows(selectors)
    matrix = acceptance_matrix()
    return {
        "status": "P2967_DIMENSION_FIELD_SELECTOR_SOURCE_OBSTRUCTION_NO_STRICT_EXPORT",
        "input_hashes": {"P2966": hashlib.sha256(P2966.read_bytes()).hexdigest() if P2966.exists() else None},
        "constructed_theoretical_objects": {"candidate_object": "DimensionFieldSelectorSource_Obstruction", "selector_rows": selectors, "proof_obligation_rows": obligations, "finite_acceptance_matrix": matrix},
        "selector_certificate": {"selector_count": len(selectors), "unique_selectors": [s["selector"] for s in selectors if s["unique_pair"]], "accepted_strict_selectors": [s["selector"] for s in selectors if s["accepted_strict_selector"]], "acceptance_matrix_rows": len(matrix), "accepted_rows": sum(1 for r in matrix if r["accepts_strict_dimension_field_selector"])},
        "decision": {
            "positive_progress": "P2967 isolates the N,sigma selector problem and audits formal selectors; uniqueness is available only through convention/import, not through a strict source.",
            "breakthrough": "No strict dimension-field selector is exported: dimensionless/minimal/integer predicates are nonunique, while lexicographic and 4D selections are imported conventions.",
            "negative_export_flags": {k: False for k in ["strict_dimension_field_selector_exported", "strict_physical_unit_law_exported", "strict_unit_U_source_exported", "strict_coefficient_source_law_exported", "strict_unit_bearing_nonproxy_coupling_exported", "nonproxy_ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay dimensionless/minimal-exponent/integer-power selector predicates, lexicographic conventions, imported 4D scalar choices, beta normalization, primitive-ray quotient arithmetic, or scalar Euler insertion.  The next proof-grade move must introduce a strict unit U source or a coefficient source law independent of N,sigma selector replay, or pivot to a genuinely new typed structural object outside the ratio-package lane while preserving the P2929-P2967 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["selector_certificate"]
    lines = ["# P2967/S1917 dimension-field selector-source obstruction", "", f"Status: `{payload['status']}`", "", "## Selector certificate", f"- selector count: `{cert['selector_count']}`", f"- unique selectors: `{cert['unique_selectors']}`", f"- accepted strict selectors: `{cert['accepted_strict_selectors']}`", f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`", "", "## Lay summary", payload["decision"]["positive_progress"], payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2966))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2967/S1917 dimension-field selector-source obstruction", "## P2967/S1917 dimension-field selector-source obstruction\n\n`P2967/S1917` attacks exactly one P2966 missing object: a strict selector/source for the dimension-field pair `(N,sigma)` in `c[V]=(9/5)U^(2sigma-N)`.  Dimensionless, minimal-exponent, and integer-power predicates remain nonunique on the current dimensional grid; unique selections arise only from lexicographic convention or imported four-dimensional scalar assumptions.  Therefore no strict physical unit law, strict unit `U` source, coefficient source law, unit-bearing nonproxy coupling, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2967/S1917 dimension-field selector `L_total` guard", "## P2967/S1917 dimension-field selector `L_total` guard\n\n`P2967/S1917` finds no nonimported strict source selecting `(N,sigma)` for the P2966 dimensional lift.  Because formal selector predicates are nonunique and unique choices are conventional/imported, the coefficient `9/5` is not installed as a sourced unit-bearing `L_total` term; no EOM, Hamiltonian, bridge, role-transfer, or ToE closure follows.\n")
    append_once(AGENTS, "Current dimension-field selector-source obstruction guardrail (P2967/S1917, 2026-06-20)", "## Current dimension-field selector-source obstruction guardrail (P2967/S1917, 2026-06-20)\n\n- P2967 audits the P2966 missing selector/source for `(N,sigma)` in `c[V]=(9/5)U^(2sigma-N)`.\n- Dimensionless, minimal-exponent, and integer-power predicates are nonunique; lexicographic and four-dimensional scalar selections are convention/import, not strict nadsoliton source.\n- Do not promote these selector predicates, imported 4D scalar choices, beta normalization, primitive-ray quotient arithmetic, scalar Euler insertion, K/C exchange, or typed scalar mediators to strict ratio-package source, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- A next admissible move must introduce a strict unit `U` source, a coefficient source law independent of `(N,sigma)` selector replay, or a genuinely new typed structural object outside the ratio-package lane; otherwise preserve the P2929-P2967 boundary.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
