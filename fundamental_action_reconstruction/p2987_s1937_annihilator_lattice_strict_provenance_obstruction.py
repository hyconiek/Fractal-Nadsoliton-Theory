#!/usr/bin/env python3
"""P2987/S1937: annihilator-lattice strict provenance obstruction.

P2986 left four possible theorem routes for the new Z12 annihilator-lattice
object.  This audit attacks exactly one route: strict nadsoliton provenance for
the annihilator lattice.  It does not replay nilradical, CRT, zero-derivation,
ratio-package, incidence, selector, bridge, role-transfer, or L_total lanes.

The finite side remains exact: the six ideals and their annihilators are
recomputed, double annihilators close, and order reversal holds.  The provenance
side remains blocked: the construction uses imported ring addition and
multiplication on Z/12Z, translation of proper ideals is not ideal-preserving,
and no current artifact exports a strict nadsoliton source map producing the
ideal lattice or a nonpremise localizer selecting its rows.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2986_s1936_z12_annihilator_lattice_source_candidate_obstruction import OUT as P2986, MODULUS, annihilator_lattice_witness

OUT = GEN / "p2987_s1937_annihilator_lattice_strict_provenance_obstruction.json"
MD = GEN / "p2987_s1937_annihilator_lattice_strict_provenance_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def provenance_witness() -> dict[str, Any]:
    lattice = annihilator_lattice_witness()
    rows = []
    for row in lattice["annihilator_rows"]:
        ideal = set(row["ideal"])
        translated_sets = {tuple(sorted(((x + t) % MODULUS for x in ideal))) for t in range(MODULUS)}
        ideal_translates = [list(s) for s in sorted(translated_sets) if any(set(s) == set(other["ideal"]) for other in lattice["annihilator_rows"])]
        rows.append({
            "ideal": row["ideal"],
            "generators": row["generators"],
            "annihilator": row["annihilator"],
            "finite_annihilator_exact": row["annihilator_product_zero"] and row["double_annihilator_returns_ideal"],
            "translation_orbit_size": len(translated_sets),
            "translation_ideal_translate_count": len(ideal_translates),
            "translation_preserves_ideal_row": len(ideal_translates) == len(translated_sets),
            "strict_nadsoliton_source_map_exported": False,
            "nonpremise_row_localizer_exported": False,
            "accepted_strict_provenance": False,
        })
    return {
        "modulus": MODULUS,
        "lattice": lattice,
        "provenance_rows": rows,
        "proper_rows_translation_preserved": all(r["translation_preserves_ideal_row"] for r in rows if len(r["ideal"]) not in (1, MODULUS)),
        "finite_annihilator_rows_exact": all(r["finite_annihilator_exact"] for r in rows),
    }


def obligation_rows(witness: dict[str, Any]) -> list[dict[str, Any]]:
    rows = witness["provenance_rows"]
    return [
        {"obligation": "finite_annihilator_lattice_recomputed", "satisfied": witness["finite_annihilator_rows_exact"], "evidence": "all six rows retain annihilator product-zero and double-annihilator exactness"},
        {"obligation": "order_reversing_lattice_retained", "satisfied": witness["lattice"]["all_order_reversing_checks_pass"], "evidence": "P2986 order-reversing inclusion checks are recomputed and pass"},
        {"obligation": "proper_ideal_translation_invariance", "satisfied": witness["proper_rows_translation_preserved"], "evidence": "proper ideals are not stable under all Z12 translations as ideal rows"},
        {"obligation": "strict_nadsoliton_source_map", "satisfied": any(r["strict_nadsoliton_source_map_exported"] for r in rows), "evidence": "no current strict source map exports ring ideals/annihilators from nadsoliton data"},
        {"obligation": "nonpremise_row_localizer", "satisfied": any(r["nonpremise_row_localizer_exported"] for r in rows), "evidence": "no row is selected by a nonpremise physical/source localizer"},
        {"obligation": "accepted_strict_provenance", "satisfied": any(r["accepted_strict_provenance"] for r in rows), "evidence": "all current rows fail strict provenance premises"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_lattice", "order_reversing", "translation_safety", "strict_source_map", "row_localizer", "source_atom_readiness"]
    return [{"present": dict(zip(names, bits)), "accepts_annihilator_strict_provenance": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p2986_path: Any) -> dict[str, Any]:
    witness = provenance_witness()
    obligations = obligation_rows(witness)
    matrix = acceptance_matrix()
    return {
        "status": "P2987_ANNIHILATOR_LATTICE_STRICT_PROVENANCE_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P2986": hashlib.sha256(p2986_path.read_bytes()).hexdigest() if p2986_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "AnnihilatorLatticeStrictProvenance_ObstructionMatrix",
            "strict_provenance_witness": witness,
            "proof_obligation_rows": obligations,
            "finite_acceptance_matrix": matrix,
        },
        "provenance_certificate": {
            "ideal_count": witness["lattice"]["ideal_count"],
            "finite_annihilator_rows_exact": witness["finite_annihilator_rows_exact"],
            "order_reversing_checks_pass": witness["lattice"]["all_order_reversing_checks_pass"],
            "proper_rows_translation_preserved": witness["proper_rows_translation_preserved"],
            "accepted_current_provenance_rows": [r["ideal"] for r in witness["provenance_rows"] if r["accepted_strict_provenance"]],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_annihilator_strict_provenance"]),
        },
        "decision": {
            "positive_progress": "P2987 attacks exactly one P2986 missing theorem route: strict nadsoliton provenance for the Z12 annihilator lattice.",
            "breakthrough": "Bounded no-go: the finite annihilator lattice remains exact, but provenance is blocked because the construction imports Z12 ring operations, proper ideals are not translation-stable source sectors, and no strict nadsoliton source map or nonpremise row localizer is exported.",
            "negative_export_flags": {k: False for k in ["annihilator_lattice_strict_provenance_exported", "source_localizer_exported", "named_atom_coupling_exported", "damping_source_exported", "unit_bearing_density_exported", "nonproxy_ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay annihilator strict-provenance rows, nilradical/CRT/zero-derivation lanes, selector replay, bridge maps, role transfer, or L_total placeholders.  A next admissible annihilator-lattice move may attack exactly one remaining route: nonpremise source-localizer, named source-atom coupling, or action installation; otherwise introduce a genuinely new strict typed object/provider or preserve the P2929-P2987 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["provenance_certificate"]
    lines = [
        "# P2987/S1937 annihilator-lattice strict provenance obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Provenance certificate",
        f"- ideal count: `{cert['ideal_count']}`",
        f"- finite annihilator rows exact: `{cert['finite_annihilator_rows_exact']}`",
        f"- order-reversing checks pass: `{cert['order_reversing_checks_pass']}`",
        f"- proper rows translation-preserved: `{cert['proper_rows_translation_preserved']}`",
        f"- accepted current provenance rows: `{cert['accepted_current_provenance_rows']}`",
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
    read_json(P2986)
    payload = build_payload(P2986)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2987/S1937 annihilator-lattice strict provenance obstruction", "## P2987/S1937 annihilator-lattice strict provenance obstruction\n\n`P2987/S1937` attacks exactly one P2986 missing theorem route: strict nadsoliton provenance for the `Z/12Z` annihilator lattice.  The finite side remains exact: the six ideal rows retain product-zero annihilators, `Ann(Ann(I))=I`, and order reversal.  The provenance side is blocked: the construction imports ring addition/multiplication on `Z/12Z`, proper ideals are not translation-stable source sectors, and no strict nadsoliton source map or nonpremise row localizer exports the ideal lattice.  No annihilator-lattice strict provenance, source-localizer, named source-atom coupling, unit-bearing action density, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2987/S1937 annihilator provenance `L_total` guard", "## P2987/S1937 annihilator provenance `L_total` guard\n\n`P2987/S1937` does not add an annihilator-provenance term to `L_total`.  It preserves the finite annihilator lattice but finds no strict nadsoliton source map, field provenance, unit-bearing measure/density, nonproxy variational chain, EOM/Hamiltonian term, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current annihilator-lattice strict provenance obstruction guardrail (P2987/S1937, 2026-06-20)", "## Current annihilator-lattice strict provenance obstruction guardrail (P2987/S1937, 2026-06-20)\n\n- P2987 attacks exactly one P2986 missing theorem route: strict nadsoliton provenance for the `Z/12Z` annihilator lattice.\n- Finite positives remain real: all six ideal rows retain product-zero annihilators, `Ann(Ann(I))=I`, and order reversal.\n- The current route is bounded no-go: the construction still imports `Z/12Z` ring operations, proper ideals are not translation-stable source sectors, and no strict nadsoliton source map or nonpremise row localizer is exported.\n- Do not promote annihilator strict-provenance rows to source-localizer, named source-atom coupling, action installation, selector closure, bridge closure, role transfer, nonproxy `L_total`, or ToE.  A next admissible annihilator-lattice move may attack exactly one remaining route (nonpremise source-localizer, named source-atom coupling, or action installation), or else introduce a genuinely new strict typed object/provider while preserving the P2929-P2987 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
