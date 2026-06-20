#!/usr/bin/env python3
"""P2966/S1916: physical unit law dimensional obstruction.

P2965 left the primitive-mean quotient as dimensionless arithmetic.  This audit
attacks exactly the next missing object: a strict physical unit law that would
make the quotient unit-bearing in the P2964 schema c[V] Phi^2 dVol_N.

The finite calculation enumerates a dimensional exponent grid.  If dVol_N has
length dimension N and Phi has length dimension -sigma, a dimensionless action
requires c[V] to have length exponent k=2*sigma-N.  The primitive mean 9/5 can
multiply every row, but the current artifacts do not select N, sigma, a length
unit U, or a strict coefficient source law.  Thus the quotient is compatible with
many unit assignments but does not source one.
"""
from __future__ import annotations

import hashlib, json
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2965_s1915_primitive_mean_scale_unit_quotient_candidate import OUT as P2965

OUT = GEN / "p2966_s1916_physical_unit_law_dimensional_obstruction.json"
MD = GEN / "p2966_s1916_physical_unit_law_dimensional_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

PRIMITIVE_MEAN = Fraction(9, 5)


def dimensional_rows() -> list[dict[str, Any]]:
    rows = []
    sigmas = [Fraction(n, 2) for n in range(0, 13)]
    for dimension in range(1, 13):
        for sigma in sigmas:
            exponent = 2 * sigma - dimension
            rows.append({
                "spacetime_dimension_N": dimension,
                "field_length_dimension_sigma": f"{sigma.numerator}/{sigma.denominator}",
                "required_coefficient_length_exponent_k": f"{exponent.numerator}/{exponent.denominator}",
                "formal_coefficient": f"(9/5) * U^({exponent.numerator}/{exponent.denominator})",
                "dimensionless_coefficient_case": exponent == 0,
                "integer_exponent": exponent.denominator == 1,
            })
    return rows


def candidate_unit_laws(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    zero_rows = [r for r in rows if r["dimensionless_coefficient_case"]]
    integer_rows = [r for r in rows if r["integer_exponent"]]
    laws = [
        {
            "candidate": "free_dimensional_completion_family",
            "description": "choose any N, sigma, and unit U, then set c=(9/5)U^(2sigma-N)",
            "nonempty": bool(rows),
            "unique_unit_law": False,
            "selects_N_sigma": False,
            "selects_length_unit_U": False,
            "strict_coefficient_source_law": False,
            "nonproxy_coupling_installed": False,
            "witness_count": len(rows),
        },
        {
            "candidate": "dimensionless_coefficient_subfamily",
            "description": "impose 2sigma=N so c=9/5 is dimensionless",
            "nonempty": bool(zero_rows),
            "unique_unit_law": len(zero_rows) == 1,
            "selects_N_sigma": False,
            "selects_length_unit_U": False,
            "strict_coefficient_source_law": False,
            "nonproxy_coupling_installed": False,
            "witness_count": len(zero_rows),
        },
        {
            "candidate": "integer_power_unit_subfamily",
            "description": "restrict coefficient to integer powers of U",
            "nonempty": bool(integer_rows),
            "unique_unit_law": len(integer_rows) == 1,
            "selects_N_sigma": False,
            "selects_length_unit_U": False,
            "strict_coefficient_source_law": False,
            "nonproxy_coupling_installed": False,
            "witness_count": len(integer_rows),
        },
        {
            "candidate": "completed_strict_physical_unit_law_schema",
            "description": "a future theorem selects N, sigma, U, and c_N[V] as a sourced coefficient",
            "nonempty": True,
            "unique_unit_law": True,
            "selects_N_sigma": True,
            "selects_length_unit_U": True,
            "strict_coefficient_source_law": True,
            "nonproxy_coupling_installed": True,
            "witness_count": 0,
        },
    ]
    for law in laws:
        current_artifact_available = law["candidate"] != "completed_strict_physical_unit_law_schema"
        law["current_artifact_available"] = current_artifact_available
        law["developmental_obstruction_witness"] = current_artifact_available and law["nonempty"] and not law["unique_unit_law"]
        law["strict_physical_unit_law_exported"] = current_artifact_available and all([
            law["unique_unit_law"],
            law["selects_N_sigma"],
            law["selects_length_unit_U"],
            law["strict_coefficient_source_law"],
            law["nonproxy_coupling_installed"],
        ])
    return laws


def obligation_rows(rows: list[dict[str, Any]], laws: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {"obligation": "dimensional_grid_constructed", "satisfied": True, "evidence": f"{len(rows)} (N,sigma) rows audited"},
        {"obligation": "primitive_mean_can_be_formally_unit_lifted", "satisfied": bool(rows), "evidence": "each row admits c=(9/5)U^k"},
        {"obligation": "unique_dimension_field_pair_selected", "satisfied": any(l["strict_physical_unit_law_exported"] for l in laws), "evidence": "current families have many rows and no selector for N,sigma"},
        {"obligation": "strict_length_unit_U_selected", "satisfied": False, "evidence": "U remains a free unit token in every current row"},
        {"obligation": "strict_coefficient_source_law_exported", "satisfied": False, "evidence": "9/5 is not promoted from primitive mean to sourced action coefficient"},
        {"obligation": "nonproxy_coupling_installed", "satisfied": False, "evidence": "P2964/P2965 installation remains pending"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["dimensional_grid", "formal_unit_lift", "unique_N_sigma", "strict_unit_U", "coefficient_source_law", "nonproxy_coupling"]
    full = (1 << len(names)) - 1
    return [{"mask": m, "present": {n: bool(m & (1 << i)) for i, n in enumerate(names)}, "accepts_strict_physical_unit_law": m == full} for m in range(1 << len(names))]


def build_payload(p2965: dict[str, Any]) -> dict[str, Any]:
    rows = dimensional_rows()
    laws = candidate_unit_laws(rows)
    obligations = obligation_rows(rows, laws)
    matrix = acceptance_matrix()
    exponents = sorted({r["required_coefficient_length_exponent_k"] for r in rows})
    return {
        "status": "P2966_PHYSICAL_UNIT_LAW_DIMENSIONAL_OBSTRUCTION_NO_STRICT_EXPORT",
        "input_hashes": {"P2965": hashlib.sha256(P2965.read_bytes()).hexdigest() if P2965.exists() else None},
        "constructed_theoretical_objects": {"candidate_object": "PhysicalUnitLaw_DimensionalExponentObstruction", "dimensional_exponent_rows": rows, "candidate_unit_law_rows": laws, "proof_obligation_rows": obligations, "finite_acceptance_matrix": matrix},
        "dimensional_certificate": {"row_count": len(rows), "distinct_exponent_count": len(exponents), "dimensionless_rows": sum(1 for r in rows if r["dimensionless_coefficient_case"]), "integer_exponent_rows": sum(1 for r in rows if r["integer_exponent"]), "strict_physical_unit_laws": [l["candidate"] for l in laws if l["strict_physical_unit_law_exported"]], "acceptance_matrix_rows": len(matrix), "accepted_rows": sum(1 for r in matrix if r["accepts_strict_physical_unit_law"])},
        "decision": {
            "positive_progress": "P2966 turns the missing physical-unit law into an explicit dimensional equation: c[V] Phi^2 dVol_N requires c[V]=(9/5)U^(2sigma-N).",
            "breakthrough": "No strict physical unit law is exported: many (N,sigma,k) completions are compatible, U remains free, and no source law selects the coefficient as a nonproxy action term.",
            "negative_export_flags": {k: False for k in ["strict_physical_unit_law_exported", "strict_coefficient_source_law_exported", "strict_unit_bearing_nonproxy_coupling_exported", "strict_ratio_package_source_exported", "damping_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay dimensional exponent grids, dimensionless subfamilies, beta length normalization, primitive-ray quotient arithmetic, scalar Euler insertion, K/C exchange, or typed scalar mediators.  The next proof-grade move must introduce exactly one strict selector/source for N and sigma, or a strict unit U source, or a coefficient source law installing 9/5 into the P2964 nonproxy coupling; otherwise pivot outside the ratio-package lane while preserving the P2929-P2966 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["dimensional_certificate"]
    lines = ["# P2966/S1916 physical unit law dimensional obstruction", "", f"Status: `{payload['status']}`", "", "## Dimensional certificate", f"- row count: `{cert['row_count']}`", f"- distinct exponent count: `{cert['distinct_exponent_count']}`", f"- dimensionless rows: `{cert['dimensionless_rows']}`", f"- integer exponent rows: `{cert['integer_exponent_rows']}`", f"- strict physical unit laws: `{cert['strict_physical_unit_laws']}`", f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`", "", "## Lay summary", payload["decision"]["positive_progress"], payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2965))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2966/S1916 physical unit law dimensional obstruction", "## P2966/S1916 physical unit law dimensional obstruction\n\n`P2966/S1916` attacks exactly one P2965 missing object: a strict physical unit law for the primitive-mean quotient.  For a candidate action density `c[V] Phi^2 dVol_N`, the dimensional equation gives `c[V]=(9/5)U^(2sigma-N)`.  The finite grid shows many compatible `(N,sigma,k)` completions, including dimensionless and integer-exponent subfamilies, but current artifacts do not select `N`, `sigma`, a length unit `U`, or a coefficient source law.  Therefore no strict unit-bearing nonproxy coupling, ratio-package source, damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2966/S1916 physical-unit dimensional `L_total` guard", "## P2966/S1916 physical-unit dimensional `L_total` guard\n\n`P2966/S1916` shows that the P2965 coefficient `9/5` can be formally dressed as `(9/5)U^(2sigma-N)` for many choices of dimension and field scaling.  Because the current theory does not source `N`, `sigma`, `U`, or the coefficient/source map, no sourced damping coefficient enters `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE from this dimensional lift.\n")
    append_once(AGENTS, "Current physical unit law dimensional obstruction guardrail (P2966/S1916, 2026-06-20)", "## Current physical unit law dimensional obstruction guardrail (P2966/S1916, 2026-06-20)\n\n- P2966 attacks the P2965 missing physical-unit law by deriving the dimensional equation for `c[V] Phi^2 dVol_N`: the primitive mean can be dressed as `c[V]=(9/5)U^(2sigma-N)`.\n- The finite grid contains many compatible `(N,sigma,k)` completions, so dimensional compatibility is not a strict source; current artifacts do not select `N`, `sigma`, `U`, or a coefficient source law.\n- Do not promote dimensional exponent grids, dimensionless subfamilies, beta normalization, primitive-ray quotient arithmetic, scalar Euler insertion, K/C exchange, or typed scalar mediators to strict ratio-package source, damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- A next admissible move must introduce exactly one strict selector/source for `N` and `sigma`, a strict unit `U` source, or a coefficient source law installing `9/5` into the P2964 nonproxy coupling; otherwise preserve the P2929-P2966 boundary.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
