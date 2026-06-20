#!/usr/bin/env python3
"""P2993/S1943: Fourier-character strict provenance obstruction.

P2992 left three Fourier-character routes.  This audit attacks exactly one:
strict character provenance for the Z/12Z additive Fourier-character table.  It
does not replay the P2992 frequency-localizer attempt, source-coupling, action
installation, annihilator/nilradical/CRT/zero-derivation lanes, selector closure,
bridge completion, role transfer, or L_total promotion.

The finite calculation verifies that every character row is an exact additive
homomorphism into exponent addition mod 12 and retains the P2991 orthogonality
certificate.  The obstruction is provenance: the rows are still imported
Fourier labels on Z/12Z, not exported as a strict nadsoliton source map or as a
nonpremise internal character of the primordial self-coupled information state.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2991_s1941_z12_additive_fourier_character_source_candidate_obstruction import MODULUS, fourier_character_witness
from p2992_s1942_fourier_character_nonpremise_frequency_source_localizer_obstruction import OUT as P2992

OUT = GEN / "p2993_s1943_fourier_character_strict_provenance_obstruction.json"
MD = GEN / "p2993_s1943_fourier_character_strict_provenance_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def homomorphism_defects(k: int) -> list[dict[str, int]]:
    defects = []
    for x in range(MODULUS):
        for y in range(MODULUS):
            lhs = (k * ((x + y) % MODULUS)) % MODULUS
            rhs = ((k * x) + (k * y)) % MODULUS
            if lhs != rhs:
                defects.append({"x": x, "y": y, "lhs": lhs, "rhs": rhs})
    return defects


def provenance_witness() -> dict[str, Any]:
    fourier = fourier_character_witness()
    rows = []
    for row in fourier["character_rows"]:
        defects = homomorphism_defects(row["k"])
        rows.append({
            "k": row["k"],
            "conductor": row["conductor"],
            "homomorphism_defect_count": len(defects),
            "exact_additive_character_row": len(defects) == 0,
            "p2992_frequency_localizer_available": False,
            "strict_nadsoliton_source_map_exported": False,
            "nonpremise_internal_character_provenance": False,
            "accepted_strict_character_provenance": False,
        })
    return {
        "modulus": MODULUS,
        "row_count": len(rows),
        "provenance_rows": rows,
        "all_rows_exact_additive_characters": all(r["exact_additive_character_row"] for r in rows),
        "p2991_orthogonality_retained": fourier["all_orthogonality_checks_pass"],
        "accepted_strict_character_provenance_rows": [r["k"] for r in rows if r["accepted_strict_character_provenance"]],
    }


def obligation_rows(witness: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"obligation": "finite_character_rows_present", "satisfied": witness["row_count"] == MODULUS, "evidence": "all 12 Fourier-character rows are present"},
        {"obligation": "exact_additive_homomorphism_rows", "satisfied": witness["all_rows_exact_additive_characters"], "evidence": "all 12 rows pass all 144 x,y homomorphism checks"},
        {"obligation": "orthogonality_certificate_retained", "satisfied": witness["p2991_orthogonality_retained"], "evidence": "the P2991 12x12 orthogonality certificate is retained"},
        {"obligation": "frequency_localizer_available", "satisfied": any(r["p2992_frequency_localizer_available"] for r in witness["provenance_rows"]), "evidence": "P2992 found no nonpremise frequency/source localizer"},
        {"obligation": "strict_nadsoliton_source_map_exported", "satisfied": any(r["strict_nadsoliton_source_map_exported"] for r in witness["provenance_rows"]), "evidence": "no source map exports a Fourier row from the primordial nadsoliton state"},
        {"obligation": "nonpremise_internal_character_provenance", "satisfied": any(r["nonpremise_internal_character_provenance"] for r in witness["provenance_rows"]), "evidence": "character rows remain imported Z/12Z Fourier labels, not sourced internal characters"},
        {"obligation": "accepted_current_strict_character_provenance", "satisfied": bool(witness["accepted_strict_character_provenance_rows"]), "evidence": "no current row satisfies the full strict provenance profile"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_rows", "homomorphism", "orthogonality", "frequency_localizer", "strict_source_map", "internal_provenance", "nonproxy_export"]
    return [{"present": dict(zip(names, bits)), "accepts_strict_character_provenance": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p2992_path: Any) -> dict[str, Any]:
    witness = provenance_witness()
    matrix = acceptance_matrix()
    return {
        "status": "P2993_FOURIER_CHARACTER_STRICT_PROVENANCE_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P2992": hashlib.sha256(p2992_path.read_bytes()).hexdigest() if p2992_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "FourierCharacterStrictProvenance_ObstructionMatrix",
            "provenance_witness": witness,
            "proof_obligation_rows": obligation_rows(witness),
            "finite_acceptance_matrix": matrix,
        },
        "provenance_certificate": {
            "row_count": witness["row_count"],
            "all_rows_exact_additive_characters": witness["all_rows_exact_additive_characters"],
            "p2991_orthogonality_retained": witness["p2991_orthogonality_retained"],
            "accepted_strict_character_provenance_rows": witness["accepted_strict_character_provenance_rows"],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_strict_character_provenance"]),
        },
        "decision": {
            "positive_progress": "P2993 attacks exactly one remaining Fourier-character route after P2992: strict character provenance.",
            "breakthrough": "Bounded no-go: every row is an exact additive character and the orthogonality certificate is retained, but no frequency localizer, strict nadsoliton source map, or nonpremise internal character provenance is exported.",
            "negative_export_flags": {k: False for k in ["strict_character_provenance_exported", "frequency_source_localizer_exported", "source_coupling_exported", "unit_bearing_density_exported", "self_learning_update_law_exported", "bridge_closure_exported", "nonproxy_ltotal_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay exact character homomorphism rows, orthogonality, P2992 localizer labels, selector closure, bridge maps, role transfer, or L_total placeholders as provenance.  A next admissible Fourier-character move may attack exactly one remaining route: source-coupling theorem or unit-bearing action installation; otherwise introduce a genuinely new strict typed object/provider or preserve the P2929-P2993 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["provenance_certificate"]
    lines = [
        "# P2993/S1943 Fourier-character strict provenance obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Provenance certificate",
        f"- row count: `{cert['row_count']}`",
        f"- all rows exact additive characters: `{cert['all_rows_exact_additive_characters']}`",
        f"- P2991 orthogonality retained: `{cert['p2991_orthogonality_retained']}`",
        f"- accepted strict character provenance rows: `{cert['accepted_strict_character_provenance_rows']}`",
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
    read_json(P2992)
    payload = build_payload(P2992)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2993/S1943 Fourier-character strict provenance obstruction", "## P2993/S1943 Fourier-character strict provenance obstruction\n\n`P2993/S1943` attacks exactly one remaining Fourier-character route after P2992: strict character provenance.  The finite side is exact: every `Z/12Z` row is an additive character under all `144` `x,y` homomorphism checks, and the P2991 orthogonality certificate is retained.  The provenance side remains blocked because P2992 exported no nonpremise frequency/source localizer, no strict nadsoliton source map exports a Fourier row from the primordial information state, and the rows remain imported spectral labels rather than internally sourced characters.  No strict character provenance, self-learning update law, source-coupling theorem, unit-bearing density, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2993/S1943 Fourier-character provenance `L_total` guard", "## P2993/S1943 Fourier-character provenance `L_total` guard\n\n`P2993/S1943` adds no Fourier-character provenance term to `L_total`.  Exact homomorphism and orthogonality checks are spectral consistency data only; they do not supply strict field provenance, self-learning update functional, named unit-bearing density, variational chain, EOM/Hamiltonian term, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current Fourier-character strict provenance obstruction guardrail (P2993/S1943, 2026-06-20)", "## Current Fourier-character strict provenance obstruction guardrail (P2993/S1943, 2026-06-20)\n\n- P2993 attacks exactly one remaining Fourier-character route after P2992: strict character provenance.\n- Finite positives are exact but spectral only: every row is an additive character under all `144` `x,y` homomorphism checks, and the P2991 orthogonality certificate is retained.\n- The current route is bounded no-go: P2992 exported no nonpremise frequency/source localizer, no strict nadsoliton source map exports a Fourier row from the primordial information state, and the rows remain imported spectral labels.\n- Do not promote exact character homomorphism, orthogonality, or provenance labels to self-learning update law, source coupling, selector closure, bridge closure, role transfer, nonproxy `L_total`, or ToE.\n- A next admissible Fourier-character move may attack exactly one remaining route (source-coupling theorem or unit-bearing action installation), or else introduce a genuinely new strict typed object/provider while preserving the P2929-P2993 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
