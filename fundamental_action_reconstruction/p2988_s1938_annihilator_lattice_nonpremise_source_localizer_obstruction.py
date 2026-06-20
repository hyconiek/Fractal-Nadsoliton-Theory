#!/usr/bin/env python3
"""P2988/S1938: annihilator-lattice nonpremise source-localizer obstruction.

P2987 left three annihilator-lattice routes.  This audit attacks exactly one:
nonpremise source-localizer for the Z12 annihilator-lattice rows.  It does not
replay strict provenance, named source-atom coupling, action installation,
nilradical/CRT/zero-derivation lanes, selector closure, bridge completion, role
transfer, or L_total promotion.

The finite calculation builds algebraic localizer signatures for every ideal row
using ideal size, annihilator size, generator gcd class, and double-annihilator
closure.  These signatures distinguish the six rows.  The obstruction is that
all such localizers are algebraic/cardinality labels imported from Z/12Z ring
bookkeeping; none is exported as a nonpremise physical/source sector or strict
nadsoliton-local row selector.
"""
from __future__ import annotations

import hashlib, json, math
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2987_s1937_annihilator_lattice_strict_provenance_obstruction import OUT as P2987
from p2986_s1936_z12_annihilator_lattice_source_candidate_obstruction import MODULUS, annihilator_lattice_witness

OUT = GEN / "p2988_s1938_annihilator_lattice_nonpremise_source_localizer_obstruction.json"
MD = GEN / "p2988_s1938_annihilator_lattice_nonpremise_source_localizer_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def gcd_signature(generators: list[int]) -> int:
    values = [math.gcd(g, MODULUS) for g in generators if g != 0]
    return min(values) if values else MODULUS


def localizer_witness() -> dict[str, Any]:
    lattice = annihilator_lattice_witness()
    rows = []
    for row in lattice["annihilator_rows"]:
        signature = {
            "ideal_size": row["size"],
            "annihilator_size": row["annihilator_size"],
            "generator_gcd_min": gcd_signature(row["generators"]),
            "double_annihilator_exact": row["double_annihilator_returns_ideal"],
        }
        rows.append({
            "ideal": row["ideal"],
            "annihilator": row["annihilator"],
            "algebraic_localizer_signature": signature,
            "signature_key": tuple(signature.items()),
            "row_distinguished_by_signature": True,
            "nonpremise_physical_sector": False,
            "strict_nadsoliton_row_selector": False,
            "source_localizer_exported": False,
        })
    seen = {}
    for row in rows:
        seen.setdefault(row["signature_key"], []).append(row["ideal"])
    for row in rows:
        row["row_distinguished_by_signature"] = len(seen[row["signature_key"]]) == 1
        row.pop("signature_key")
    return {
        "modulus": MODULUS,
        "row_count": len(rows),
        "localizer_rows": rows,
        "all_rows_algebraically_distinguished": all(r["row_distinguished_by_signature"] for r in rows),
        "exported_nonpremise_source_localizers": [r["ideal"] for r in rows if r["source_localizer_exported"]],
    }


def obligation_rows(witness: dict[str, Any]) -> list[dict[str, Any]]:
    rows = witness["localizer_rows"]
    return [
        {"obligation": "finite_annihilator_rows_present", "satisfied": witness["row_count"] == 6, "evidence": "all six P2986 annihilator-lattice rows are present"},
        {"obligation": "algebraic_signature_distinguishes_rows", "satisfied": witness["all_rows_algebraically_distinguished"], "evidence": "ideal size, annihilator size, generator gcd, and double-ann closure distinguish the six rows"},
        {"obligation": "nonpremise_physical_sector", "satisfied": any(r["nonpremise_physical_sector"] for r in rows), "evidence": "cardinality/gcd signatures are algebraic bookkeeping labels only"},
        {"obligation": "strict_nadsoliton_row_selector", "satisfied": any(r["strict_nadsoliton_row_selector"] for r in rows), "evidence": "P2987 found no strict source map exporting annihilator rows"},
        {"obligation": "source_localizer_exported", "satisfied": bool(witness["exported_nonpremise_source_localizers"]), "evidence": "no current row has a nonpremise source-localizer theorem"},
        {"obligation": "accepted_current_localizer", "satisfied": False, "evidence": "algebraic row separation does not satisfy nonpremise source-localizer premises"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_rows", "algebraic_distinction", "nonpremise_sector", "strict_row_selector", "source_theorem", "coupling_readiness"]
    return [{"present": dict(zip(names, bits)), "accepts_nonpremise_source_localizer": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p2987_path: Any) -> dict[str, Any]:
    witness = localizer_witness()
    obligations = obligation_rows(witness)
    matrix = acceptance_matrix()
    return {
        "status": "P2988_ANNIHILATOR_LATTICE_NONPREMISE_SOURCE_LOCALIZER_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P2987": hashlib.sha256(p2987_path.read_bytes()).hexdigest() if p2987_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "AnnihilatorLatticeNonpremiseSourceLocalizer_ObstructionMatrix",
            "localizer_witness": witness,
            "proof_obligation_rows": obligations,
            "finite_acceptance_matrix": matrix,
        },
        "localizer_certificate": {
            "row_count": witness["row_count"],
            "all_rows_algebraically_distinguished": witness["all_rows_algebraically_distinguished"],
            "exported_nonpremise_source_localizers": witness["exported_nonpremise_source_localizers"],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_nonpremise_source_localizer"]),
        },
        "decision": {
            "positive_progress": "P2988 attacks exactly one remaining annihilator-lattice route after P2987: a nonpremise source-localizer for the six ideal-annihilator rows.",
            "breakthrough": "Bounded no-go: algebraic signatures distinguish all six rows, but they are cardinality/gcd bookkeeping labels and no nonpremise physical sector, strict row selector, or source-localizer theorem is exported.",
            "negative_export_flags": {k: False for k in ["annihilator_source_localizer_exported", "strict_provenance_exported", "named_atom_coupling_exported", "damping_source_exported", "unit_bearing_density_exported", "nonproxy_ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay algebraic annihilator signatures, strict provenance, nilradical/CRT/zero-derivation lanes, selector replay, bridge maps, role transfer, or L_total placeholders as source localization.  A next admissible annihilator-lattice move may attack exactly one remaining route: named source-atom coupling or action installation; otherwise introduce a genuinely new strict typed object/provider or preserve the P2929-P2988 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["localizer_certificate"]
    lines = [
        "# P2988/S1938 annihilator-lattice nonpremise source-localizer obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Localizer certificate",
        f"- row count: `{cert['row_count']}`",
        f"- all rows algebraically distinguished: `{cert['all_rows_algebraically_distinguished']}`",
        f"- exported nonpremise source localizers: `{cert['exported_nonpremise_source_localizers']}`",
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
    read_json(P2987)
    payload = build_payload(P2987)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2988/S1938 annihilator-lattice nonpremise source-localizer obstruction", "## P2988/S1938 annihilator-lattice nonpremise source-localizer obstruction\n\n`P2988/S1938` attacks exactly one remaining annihilator-lattice route after P2987: a nonpremise source-localizer for the six ideal-annihilator rows.  The finite side is positive in a limited algebraic sense: ideal size, annihilator size, generator gcd class, and double-annihilator closure distinguish all six rows.  The source-localizer side remains blocked because these signatures are cardinality/gcd bookkeeping labels, not nonpremise physical sectors or strict nadsoliton row selectors.  No source-localizer, named source-atom coupling, unit-bearing action density, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2988/S1938 annihilator localizer `L_total` guard", "## P2988/S1938 annihilator localizer `L_total` guard\n\n`P2988/S1938` distinguishes annihilator-lattice rows algebraically but exports no nonpremise source-localizer and adds no annihilator-localizer term to `L_total`.  The row signatures supply no strict field provenance, named unit-bearing density, measure theorem, nonproxy variational chain, EOM/Hamiltonian term, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current annihilator-lattice nonpremise source-localizer obstruction guardrail (P2988/S1938, 2026-06-20)", "## Current annihilator-lattice nonpremise source-localizer obstruction guardrail (P2988/S1938, 2026-06-20)\n\n- P2988 attacks exactly one remaining annihilator-lattice route after P2987: a nonpremise source-localizer for the six ideal-annihilator rows.\n- Finite positives are algebraic only: ideal size, annihilator size, generator gcd class, and double-annihilator closure distinguish all six rows.\n- The current route is bounded no-go: these signatures are cardinality/gcd bookkeeping labels, not nonpremise physical sectors or strict nadsoliton row selectors.\n- Do not promote annihilator row signatures to strict provenance, source-localizer, named source-atom coupling, action installation, selector closure, bridge closure, role transfer, nonproxy `L_total`, or ToE.  A next admissible annihilator-lattice move may attack exactly one remaining route (named source-atom coupling or action installation), or else introduce a genuinely new strict typed object/provider while preserving the P2929-P2988 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
