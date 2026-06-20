#!/usr/bin/env python3
"""P2989/S1939: annihilator-lattice named source-atom coupling obstruction.

P2988 left two annihilator-lattice routes.  This audit attacks exactly one:
named source-atom coupling for the six Z/12Z ideal-annihilator rows.  It does
not replay strict provenance, nonpremise source-localizer, action installation,
nilradical/CRT/zero-derivation lanes, selector closure, bridge completion, role
transfer, or L_total promotion.

The finite calculation constructs a four-atom coupling matrix against current
missing source atoms: selector/orientation sign, target-independent positive
beta/Z_beta damping, legacy-to-strict bridge-source, and unit-bearing
action-density source.  Each annihilator row receives exact algebraic invariants,
but every named atom is blocked by the missing source-localizer/provenance and by
atom-specific strict-source obligations.  Hence the matrix is a bounded no-go,
not an installation theorem.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2986_s1936_z12_annihilator_lattice_source_candidate_obstruction import annihilator_lattice_witness
from p2988_s1938_annihilator_lattice_nonpremise_source_localizer_obstruction import OUT as P2988

OUT = GEN / "p2989_s1939_annihilator_lattice_named_source_atom_coupling_obstruction.json"
MD = GEN / "p2989_s1939_annihilator_lattice_named_source_atom_coupling_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

NAMED_ATOMS = [
    "selector_orientation_sign",
    "target_independent_positive_beta_Z_beta",
    "legacy_to_strict_bridge_source",
    "unit_bearing_action_density_source",
]


def row_invariants(row: dict[str, Any]) -> dict[str, Any]:
    ideal = row["ideal"]
    ann = row["annihilator"]
    return {
        "ideal_size": row["size"],
        "annihilator_size": row["annihilator_size"],
        "zero_product_pairs": row["size"] * row["annihilator_size"],
        "annihilator_ratio_num": row["annihilator_size"],
        "annihilator_ratio_den": row["size"],
        "contains_unit": 1 in ideal,
        "annihilator_contains_unit": 1 in ann,
        "double_annihilator_exact": row["double_annihilator_returns_ideal"],
    }


def atom_failure(atom: str, inv: dict[str, Any]) -> str:
    if atom == "selector_orientation_sign":
        return "annihilator rows are unoriented ideal data and provide no +omega/-omega sign selector"
    if atom == "target_independent_positive_beta_Z_beta":
        return "row ratios are row-dependent cardinalities, not a target-independent positive beta/Z_beta source"
    if atom == "legacy_to_strict_bridge_source":
        return "ideal-annihilator bookkeeping supplies no amplitude, phase, or nonlinear damping completion map"
    if atom == "unit_bearing_action_density_source":
        return "zero-product support is finite algebraic support without unit-bearing measure or variational density"
    raise ValueError(atom)


def coupling_witness() -> dict[str, Any]:
    lattice = annihilator_lattice_witness()
    rows = []
    for row in lattice["annihilator_rows"]:
        inv = row_invariants(row)
        atom_rows = []
        for atom in NAMED_ATOMS:
            atom_rows.append({
                "named_atom": atom,
                "formal_receiver_available": True,
                "algebraic_invariants": inv,
                "strict_provenance_available": False,
                "nonpremise_source_localizer_available": False,
                "atom_specific_source_theorem": False,
                "accepted_named_atom_coupling": False,
                "failure": atom_failure(atom, inv),
            })
        rows.append({"ideal": row["ideal"], "annihilator": row["annihilator"], "atom_couplings": atom_rows})
    accepted = [
        {"ideal": r["ideal"], "named_atom": a["named_atom"]}
        for r in rows for a in r["atom_couplings"] if a["accepted_named_atom_coupling"]
    ]
    return {
        "named_atoms": NAMED_ATOMS,
        "row_count": len(rows),
        "coupling_row_count": sum(len(r["atom_couplings"]) for r in rows),
        "coupling_rows": rows,
        "accepted_named_atom_couplings": accepted,
    }


def obligation_rows(witness: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"obligation": "finite_annihilator_rows_present", "satisfied": witness["row_count"] == 6, "evidence": "all six ideal-annihilator rows are available"},
        {"obligation": "named_atom_matrix_constructed", "satisfied": witness["coupling_row_count"] == 24, "evidence": "six rows times four named source atoms are audited"},
        {"obligation": "formal_algebraic_receivers_available", "satisfied": True, "evidence": "each row has exact size, annihilator size, zero-product pair count, and double-annihilator status"},
        {"obligation": "strict_provenance_available", "satisfied": False, "evidence": "P2987 found no strict nadsoliton source map for annihilator rows"},
        {"obligation": "nonpremise_source_localizer_available", "satisfied": False, "evidence": "P2988 found only cardinality/gcd bookkeeping row labels"},
        {"obligation": "atom_specific_source_theorem", "satisfied": False, "evidence": "no selector sign, beta/Z_beta source, bridge map, or unit-bearing action-density theorem is exported"},
        {"obligation": "accepted_current_named_atom_coupling", "satisfied": bool(witness["accepted_named_atom_couplings"]), "evidence": "no current annihilator row satisfies the full named-atom coupling profile"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_row", "named_atom", "formal_receiver", "strict_provenance", "source_localizer", "atom_theorem", "unit_measure_or_sign", "nonproxy_export"]
    return [{"present": dict(zip(names, bits)), "accepts_named_atom_coupling": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p2988_path: Any) -> dict[str, Any]:
    witness = coupling_witness()
    matrix = acceptance_matrix()
    return {
        "status": "P2989_ANNIHILATOR_LATTICE_NAMED_SOURCE_ATOM_COUPLING_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P2988": hashlib.sha256(p2988_path.read_bytes()).hexdigest() if p2988_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "AnnihilatorLatticeNamedSourceAtomCoupling_ObstructionMatrix",
            "coupling_witness": witness,
            "proof_obligation_rows": obligation_rows(witness),
            "finite_acceptance_matrix": matrix,
        },
        "coupling_certificate": {
            "row_count": witness["row_count"],
            "named_atoms": witness["named_atoms"],
            "coupling_row_count": witness["coupling_row_count"],
            "accepted_named_atom_couplings": witness["accepted_named_atom_couplings"],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_named_atom_coupling"]),
        },
        "decision": {
            "positive_progress": "P2989 attacks exactly one remaining annihilator-lattice route after P2988: named source-atom coupling for the six ideal-annihilator rows.",
            "breakthrough": "Bounded no-go: exact finite algebraic receivers exist for all 24 row/atom tests, but no strict provenance, nonpremise source-localizer, selector sign, target-independent beta/Z_beta source, bridge-source map, or unit-bearing action-density theorem is exported.",
            "negative_export_flags": {k: False for k in ["annihilator_named_atom_coupling_exported", "strict_provenance_exported", "source_localizer_exported", "selector_closure_exported", "damping_source_exported", "bridge_closure_exported", "unit_bearing_density_exported", "nonproxy_ltotal_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay annihilator named-atom rows, algebraic row signatures, strict provenance, nilradical/CRT/zero-derivation lanes, selector replay, bridge maps, role transfer, or L_total placeholders as source coupling.  The only remaining annihilator-lattice route is action installation with a genuinely unit-bearing measure, named density theorem, boundary/integration map, and nonproxy continuum lift; otherwise introduce a genuinely new strict typed object/provider or preserve the P2929-P2989 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["coupling_certificate"]
    lines = [
        "# P2989/S1939 annihilator-lattice named source-atom coupling obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Coupling certificate",
        f"- row count: `{cert['row_count']}`",
        f"- named atoms: `{cert['named_atoms']}`",
        f"- coupling rows: `{cert['coupling_row_count']}`",
        f"- accepted named atom couplings: `{cert['accepted_named_atom_couplings']}`",
        f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`",
        "",
        "## Lay summary",
        payload["decision"]["positive_progress"],
        payload["decision"]["breakthrough"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P2988)
    payload = build_payload(P2988)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2989/S1939 annihilator-lattice named source-atom coupling obstruction", "## P2989/S1939 annihilator-lattice named source-atom coupling obstruction\n\n`P2989/S1939` attacks exactly one remaining annihilator-lattice route after P2988: named source-atom coupling for the six ideal-annihilator rows.  The finite matrix builds 24 row/atom receiver tests against selector/orientation sign, target-independent positive `beta/Z_beta`, legacy-to-strict bridge-source, and unit-bearing action-density source atoms.  Exact algebraic receivers exist, but strict provenance, nonpremise source-localizer, atom-specific source theorem, unit-bearing measure/sign, and nonproxy export are absent.  No named source-atom coupling, damping source, selector closure, bridge closure, unit-bearing action density, nonproxy `L_total`, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2989/S1939 annihilator named source-atom coupling `L_total` guard", "## P2989/S1939 annihilator named source-atom coupling `L_total` guard\n\n`P2989/S1939` does not add an annihilator named-atom coupling term to `L_total`.  The row/atom receiver matrix lacks strict field provenance, nonpremise source localization, a named unit-bearing density, measure/sign theorem, boundary/integration theorem, nonproxy variational chain, EOM/Hamiltonian term, bridge closure, role transfer, and ToE.\n")
    append_once(AGENTS, "Current annihilator-lattice named source-atom coupling obstruction guardrail (P2989/S1939, 2026-06-20)", "## Current annihilator-lattice named source-atom coupling obstruction guardrail (P2989/S1939, 2026-06-20)\n\n- P2989 attacks exactly one remaining annihilator-lattice route after P2988: named source-atom coupling for the six ideal-annihilator rows.\n- The finite matrix has 24 row/atom tests across selector/orientation sign, target-independent positive `beta/Z_beta`, legacy-to-strict bridge-source, and unit-bearing action-density source atoms; exact algebraic receivers are available.\n- The route is bounded no-go because strict provenance, nonpremise source-localizer, atom-specific source theorem, unit-bearing measure/sign, and nonproxy export are absent.\n- Do not promote annihilator named-atom rows to selector closure, damping source, bridge closure, action installation, role transfer, nonproxy `L_total`, or ToE.  The only remaining annihilator-lattice route is action installation with a genuinely unit-bearing measure, named density theorem, boundary/integration map, and nonproxy continuum lift; otherwise introduce a genuinely new strict typed object/provider while preserving the P2929-P2989 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
