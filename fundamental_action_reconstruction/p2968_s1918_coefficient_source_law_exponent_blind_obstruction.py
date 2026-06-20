#!/usr/bin/env python3
"""P2968/S1918: coefficient source-law exponent-blind obstruction.

P2967 left an honest route: a coefficient source law independent of replaying the
(N, sigma) selector.  This audit attacks that exact object.  For the P2966 lift,
each admissible dimensional row has coefficient c_k=(9/5)U^k.  A source law that
is independent of (N, sigma) must either be exponent-blind, select a unit U, or
fix k by another source.

The finite calculation shows 24 distinct k values on the audited grid.  The
primitive mean 9/5 is k-blind only by dropping units; U=1 is a convention; k=0
replays dimension/field selection.  Thus no strict coefficient source law is
exported.
"""
from __future__ import annotations

import hashlib, json
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2966_s1916_physical_unit_law_dimensional_obstruction import dimensional_rows
from p2967_s1917_dimension_field_selector_source_obstruction import OUT as P2967, parse_frac

OUT = GEN / "p2968_s1918_coefficient_source_law_exponent_blind_obstruction.json"
MD = GEN / "p2968_s1918_coefficient_source_law_exponent_blind_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def exponent_rows() -> list[dict[str, Any]]:
    counts: dict[Fraction, int] = {}
    for row in dimensional_rows():
        k = parse_frac(row["required_coefficient_length_exponent_k"])
        counts[k] = counts.get(k, 0) + 1
    return [{"exponent_k": f"{k.numerator}/{k.denominator}", "row_count": counts[k], "formal_coefficient": f"(9/5) * U^({k.numerator}/{k.denominator})"} for k in sorted(counts)]


def source_law_rows(exponents: list[dict[str, Any]]) -> list[dict[str, Any]]:
    k_count = len(exponents)
    rows = [
        {
            "candidate": "primitive_mean_exponent_blind_coefficient",
            "definition": "c=9/5 for all k",
            "independent_of_N_sigma_selector": True,
            "unit_bearing": False,
            "uses_unit_convention": False,
            "fixes_k_by_source": False,
            "strict_source_exported": False,
            "witness": f"drops {k_count} distinct unit exponents",
        },
        {
            "candidate": "unit_convention_U_equals_1",
            "definition": "c_k=(9/5)*1^k",
            "independent_of_N_sigma_selector": True,
            "unit_bearing": False,
            "uses_unit_convention": True,
            "fixes_k_by_source": False,
            "strict_source_exported": False,
            "witness": "collapses units by convention rather than source",
        },
        {
            "candidate": "dimensionless_k_zero_section",
            "definition": "select k=0 then c=9/5",
            "independent_of_N_sigma_selector": False,
            "unit_bearing": True,
            "uses_unit_convention": False,
            "fixes_k_by_source": False,
            "strict_source_exported": False,
            "witness": "replays P2967 dimension/field selector obstruction",
        },
        {
            "candidate": "completed_strict_coefficient_source_law_schema",
            "definition": "strict source chooses U and/or k and installs c_k into P2964 nonproxy action density",
            "independent_of_N_sigma_selector": True,
            "unit_bearing": True,
            "uses_unit_convention": False,
            "fixes_k_by_source": True,
            "strict_source_exported": True,
            "witness": "schema only; not currently exported",
        },
    ]
    for row in rows:
        row["current_artifact_available"] = row["candidate"] != "completed_strict_coefficient_source_law_schema"
        row["accepted_current_strict_coefficient_source"] = row["current_artifact_available"] and row["independent_of_N_sigma_selector"] and row["unit_bearing"] and not row["uses_unit_convention"] and row["fixes_k_by_source"] and row["strict_source_exported"]
    return rows


def obligation_rows(laws: list[dict[str, Any]], exponents: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {"obligation": "distinct_exponent_set_constructed", "satisfied": True, "evidence": f"{len(exponents)} distinct k exponents audited"},
        {"obligation": "coefficient_law_independent_of_N_sigma", "satisfied": any(r["independent_of_N_sigma_selector"] and r["current_artifact_available"] for r in laws), "evidence": "exponent-blind and U=1 candidates are independent but fail unit/source criteria"},
        {"obligation": "unit_bearing_without_unit_convention", "satisfied": any(r["unit_bearing"] and not r["uses_unit_convention"] and r["current_artifact_available"] for r in laws), "evidence": "only k=0 is unit-bearing but it replays selector obstruction"},
        {"obligation": "k_or_U_fixed_by_strict_source", "satisfied": any(r["fixes_k_by_source"] and r["current_artifact_available"] for r in laws), "evidence": "no current law fixes U or k by source"},
        {"obligation": "accepted_current_strict_coefficient_source", "satisfied": any(r["accepted_current_strict_coefficient_source"] for r in laws), "evidence": "completed law remains unavailable"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["exponent_set", "N_sigma_independent", "unit_bearing", "no_unit_convention", "k_or_U_source", "nonproxy_installation"]
    full = (1 << len(names)) - 1
    return [{"mask": m, "present": {n: bool(m & (1 << i)) for i, n in enumerate(names)}, "accepts_strict_coefficient_source_law": m == full} for m in range(1 << len(names))]


def build_payload(p2967: dict[str, Any]) -> dict[str, Any]:
    exponents = exponent_rows()
    laws = source_law_rows(exponents)
    obligations = obligation_rows(laws, exponents)
    matrix = acceptance_matrix()
    return {
        "status": "P2968_COEFFICIENT_SOURCE_LAW_EXPONENT_BLIND_OBSTRUCTION_NO_STRICT_EXPORT",
        "input_hashes": {"P2967": hashlib.sha256(P2967.read_bytes()).hexdigest() if P2967.exists() else None},
        "constructed_theoretical_objects": {"candidate_object": "CoefficientSourceLaw_ExponentBlindObstruction", "exponent_rows": exponents, "source_law_rows": laws, "proof_obligation_rows": obligations, "finite_acceptance_matrix": matrix},
        "coefficient_certificate": {"distinct_exponent_count": len(exponents), "source_law_count": len(laws), "accepted_current_strict_sources": [r["candidate"] for r in laws if r["accepted_current_strict_coefficient_source"]], "acceptance_matrix_rows": len(matrix), "accepted_rows": sum(1 for r in matrix if r["accepts_strict_coefficient_source_law"])},
        "decision": {
            "positive_progress": "P2968 reduces coefficient sourcing to the finite exponent obstruction c_k=(9/5)U^k over the P2966 exponent set.",
            "breakthrough": "No strict coefficient source law is exported: exponent-blind 9/5 drops units, U=1 is convention, and k=0 replays the P2967 selector obstruction.",
            "negative_export_flags": {k: False for k in ["strict_coefficient_source_law_exported", "strict_unit_U_source_exported", "strict_physical_unit_law_exported", "strict_unit_bearing_nonproxy_coupling_exported", "strict_ratio_package_source_exported", "nonproxy_ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay exponent-blind 9/5, U=1 unit convention, k=0 sections, N/sigma selectors, beta normalization, primitive-ray quotient arithmetic, or scalar Euler insertion.  The next proof-grade move must introduce a strict unit U source, construct a nonconventional nonproxy installation law for c_k that fixes U or k internally, or pivot to a genuinely new typed structural object outside the ratio-package lane while preserving the P2929-P2968 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["coefficient_certificate"]
    lines = ["# P2968/S1918 coefficient source-law exponent-blind obstruction", "", f"Status: `{payload['status']}`", "", "## Coefficient certificate", f"- distinct exponent count: `{cert['distinct_exponent_count']}`", f"- source-law count: `{cert['source_law_count']}`", f"- accepted current strict sources: `{cert['accepted_current_strict_sources']}`", f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`", "", "## Lay summary", payload["decision"]["positive_progress"], payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2967))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2968/S1918 coefficient source-law exponent-blind obstruction", "## P2968/S1918 coefficient source-law exponent-blind obstruction\n\n`P2968/S1918` attacks the P2967-admissible coefficient-source route independent of `(N,sigma)` selector replay.  Across the P2966 dimensional grid there are 24 distinct exponents `k` in `c_k=(9/5)U^k`.  Exponent-blind `9/5` drops units, `U=1` is a unit convention, and `k=0` replays the dimension/field selector obstruction.  Therefore no strict coefficient source law, strict unit `U` source, unit-bearing nonproxy coupling, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2968/S1918 coefficient-source `L_total` guard", "## P2968/S1918 coefficient-source `L_total` guard\n\n`P2968/S1918` shows that installing the P2965 coefficient as `c_k=(9/5)U^k` cannot be done by an exponent-blind rule, a `U=1` convention, or a `k=0` selector replay.  Since no strict source fixes `U` or `k`, no sourced damping coefficient enters `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current coefficient source-law exponent-blind obstruction guardrail (P2968/S1918, 2026-06-20)", "## Current coefficient source-law exponent-blind obstruction guardrail (P2968/S1918, 2026-06-20)\n\n- P2968 audits a coefficient source law independent of `(N,sigma)` selector replay for `c_k=(9/5)U^k`; the current exponent set has 24 distinct `k` values.\n- Exponent-blind `9/5` drops units, `U=1` is convention, and `k=0` replays P2967; none exports a strict coefficient source law.\n- Do not promote exponent-blind coefficients, unit conventions, dimensionless `k=0` sections, N/sigma selector replay, beta normalization, primitive-ray quotient arithmetic, or scalar Euler insertion to strict ratio-package source, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- A next admissible move must introduce a strict unit `U` source, a nonconventional nonproxy installation law fixing `U` or `k`, or a genuinely new typed structural object outside the ratio-package lane; otherwise preserve the P2929-P2968 boundary.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
