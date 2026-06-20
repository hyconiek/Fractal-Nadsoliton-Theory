#!/usr/bin/env python3
"""P2977/S1927: nilradical strict-nadsoliton provenance obstruction.

P2976 made the Z12 nilradical filtration a concrete finite typed object and
then named strict nadsoliton provenance as one exact missing theorem.  This file
attacks that theorem route only.  It does not replay ratio arithmetic, Gamma
Jacobians, selector/QW-2191 closure, generic bridge completion, incidence
formalism, unit conventions, or L_total promotion.

The finite result is deliberately narrow: the nilradical {0, 6} is canonical as
a ring ideal and fixed by all units, but current artifacts do not export a
strict nadsoliton source map that chooses the Z12 ring multiplication/zero
section as a physical source rather than a formal algebraic carrier.  Translation
of the two-point coset {0,6} gives six indistinguishable antipodal cosets, so the
ring-canonical witness is not a non-premise source-localizer or selector.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from math import gcd
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2976_s1926_z12_nilradical_filtration_source_candidate import OUT as P2976

OUT = GEN / "p2977_s1927_nilradical_strict_nadsoliton_provenance_obstruction.json"
MD = GEN / "p2977_s1927_nilradical_strict_nadsoliton_provenance_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MODULUS = 12
NILRADICAL = [0, 6]
UNITS = [x for x in range(MODULUS) if gcd(x, MODULUS) == 1]


def translate_set(subset: list[int], shift: int) -> tuple[int, ...]:
    return tuple(sorted(((x + shift) % MODULUS for x in subset)))


def finite_provenance_witness() -> dict[str, Any]:
    cosets = sorted({translate_set(NILRADICAL, t) for t in range(MODULUS)})
    unit_images = {str(u): tuple(sorted((u * x) % MODULUS for x in NILRADICAL)) for u in UNITS}
    zero_section_candidates = [x for x in range(MODULUS) if x == 0]
    return {
        "modulus": MODULUS,
        "nilradical": NILRADICAL,
        "units": UNITS,
        "unit_images_of_nilradical": unit_images,
        "all_units_fix_nilradical_set": all(tuple(NILRADICAL) == image for image in unit_images.values()),
        "translated_antipodal_cosets": [list(c) for c in cosets],
        "translated_antipodal_coset_count": len(cosets),
        "nilradical_coset_translation_stabilizer": [t for t in range(MODULUS) if translate_set(NILRADICAL, t) == tuple(NILRADICAL)],
        "source_origin_is_ring_zero_only": True,
        "ring_zero_exported_as_strict_nadsoliton_source": False,
        "zero_section_candidates_debug_count": len(zero_section_candidates),
    }


def candidate_rows(witness: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "candidate": "ring_canonical_nilradical_provenance",
            "finite_ring_canonicity": True,
            "unit_invariance": witness["all_units_fix_nilradical_set"],
            "translation_orbit_degeneracy_removed": False,
            "strict_nadsoliton_source_map_exported": False,
            "nonpremise_zero_section_source_exported": False,
            "couples_to_named_missing_source_atom": False,
            "accepted_strict_provenance_theorem": False,
            "witness": "nilradical is unique as a ring ideal, but this only proves algebraic canonicity",
        },
        {
            "candidate": "unit_fixed_nilpotent_6_provenance",
            "finite_ring_canonicity": True,
            "unit_invariance": True,
            "translation_orbit_degeneracy_removed": False,
            "strict_nadsoliton_source_map_exported": False,
            "nonpremise_zero_section_source_exported": False,
            "couples_to_named_missing_source_atom": False,
            "accepted_strict_provenance_theorem": False,
            "witness": "6 is fixed by U(12), but unit-fixity is orientation/source-blind and not a source map",
        },
        {
            "candidate": "translated_antipodal_coset_family",
            "finite_ring_canonicity": False,
            "unit_invariance": False,
            "translation_orbit_degeneracy_removed": False,
            "strict_nadsoliton_source_map_exported": False,
            "nonpremise_zero_section_source_exported": False,
            "couples_to_named_missing_source_atom": False,
            "accepted_strict_provenance_theorem": False,
            "witness": f"translation produces {witness['translated_antipodal_coset_count']} antipodal cosets with no strict localizer",
        },
        {
            "candidate": "completed_strict_nilradical_provenance_schema",
            "finite_ring_canonicity": True,
            "unit_invariance": True,
            "translation_orbit_degeneracy_removed": True,
            "strict_nadsoliton_source_map_exported": True,
            "nonpremise_zero_section_source_exported": True,
            "couples_to_named_missing_source_atom": True,
            "accepted_strict_provenance_theorem": False,
            "witness": "schema row records the missing theorem shape; no current artifact supplies it",
        },
    ]


def obligation_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    current = [r for r in rows if r["candidate"] != "completed_strict_nilradical_provenance_schema"]
    return [
        {"obligation": "finite_ring_canonicity", "satisfied": any(r["finite_ring_canonicity"] for r in current), "evidence": "the nilradical is computed intrinsically from multiplication in Z/12Z"},
        {"obligation": "unit_invariance", "satisfied": any(r["unit_invariance"] for r in current), "evidence": "all units fix the nilradical set and the element 6"},
        {"obligation": "translation_orbit_degeneracy_removed", "satisfied": any(r["translation_orbit_degeneracy_removed"] for r in current), "evidence": "six translated antipodal cosets remain if ring zero is not strictly sourced"},
        {"obligation": "strict_nadsoliton_source_map_exported", "satisfied": any(r["strict_nadsoliton_source_map_exported"] for r in current), "evidence": "no theorem maps the nadsoliton ontology to this nilradical filtration"},
        {"obligation": "nonpremise_zero_section_source_exported", "satisfied": any(r["nonpremise_zero_section_source_exported"] for r in current), "evidence": "the zero section is algebraic bookkeeping, not a strict source-localizer"},
        {"obligation": "couples_to_named_missing_source_atom", "satisfied": any(r["couples_to_named_missing_source_atom"] for r in current), "evidence": "no selector, damping, bridge-source, or action-density atom is coupled"},
        {"obligation": "accepted_strict_provenance_theorem", "satisfied": any(r["accepted_strict_provenance_theorem"] for r in current), "evidence": "current rows are developmental/provenance-obstruction rows only"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_canonicity", "unit_invariance", "translation_degeneracy_removed", "strict_source_map", "zero_section_source", "named_coupling"]
    rows = []
    for bits in product([False, True], repeat=len(names)):
        present = dict(zip(names, bits))
        rows.append({"present": present, "accepts_nilradical_strict_provenance": all(bits)})
    return rows


def build_payload(p2976_path: Any) -> dict[str, Any]:
    witness = finite_provenance_witness()
    rows = candidate_rows(witness)
    obligations = obligation_rows(rows)
    matrix = acceptance_matrix()
    return {
        "status": "P2977_NILRADICAL_STRICT_NADSOLITON_PROVENANCE_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P2976": hashlib.sha256(p2976_path.read_bytes()).hexdigest() if p2976_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "NilradicalStrictNadsolitonProvenance_ObstructionMatrix",
            "finite_provenance_witness": witness,
            "candidate_rows": rows,
            "proof_obligation_rows": obligations,
            "finite_acceptance_matrix": matrix,
        },
        "provenance_certificate": {
            "candidate_count": len(rows),
            "translated_antipodal_coset_count": witness["translated_antipodal_coset_count"],
            "translation_stabilizer": witness["nilradical_coset_translation_stabilizer"],
            "all_units_fix_nilradical_set": witness["all_units_fix_nilradical_set"],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_nilradical_strict_provenance"]),
            "accepted_current_strict_provenance_theorems": [r["candidate"] for r in rows if r["accepted_strict_provenance_theorem"]],
        },
        "decision": {
            "positive_progress": "P2977 attacks exactly one P2976 missing theorem: strict nadsoliton provenance for the Z12 nilradical filtration.  It separates finite ring canonicity/unit-invariance from strict sourcehood.",
            "breakthrough": "Bounded no-go: the nilradical is canonical inside the ring and unit-fixed, but current artifacts do not source the Z12 ring zero/multiplication from the nadsoliton, remove the translated antipodal-coset degeneracy nonpremise, or couple the object to a named source atom.",
            "negative_export_flags": {k: False for k in ["strict_nilradical_provenance_exported", "strict_nilradical_source_exported", "selector_or_orientation_exported", "damping_source_exported", "nonproxy_ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay nilradical canonicity, unit-fixity of 6, translated-coset scans, ratio/Gamma/incidence/selector/bridge lanes, or formal L_total placeholders as provenance.  The next proof-grade move may attack exactly one remaining nilradical theorem route: a named source-atom coupling or an action/variational installation; otherwise pivot to a genuinely new typed object and preserve the P2929-P2977 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["provenance_certificate"]
    lines = [
        "# P2977/S1927 nilradical strict-nadsoliton provenance obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite provenance certificate",
        f"- candidate rows: `{cert['candidate_count']}`",
        f"- translated antipodal cosets: `{cert['translated_antipodal_coset_count']}`",
        f"- translation stabilizer of `{NILRADICAL}`: `{cert['translation_stabilizer']}`",
        f"- all units fix nilradical set: `{cert['all_units_fix_nilradical_set']}`",
        f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`",
        f"- accepted current strict provenance theorems: `{cert['accepted_current_strict_provenance_theorems']}`",
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
    read_json(P2976)
    payload = build_payload(P2976)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2977/S1927 nilradical strict-nadsoliton provenance obstruction", "## P2977/S1927 nilradical strict-nadsoliton provenance obstruction\n\n`P2977/S1927` attacks exactly one P2976 missing theorem route: strict nadsoliton provenance for the `Z/12Z` nilradical filtration.  The finite matrix confirms the positive algebraic side (`{0,6}` is ring-canonical and fixed by all `4` units) but also the obstruction: without a strict source for the ring zero/multiplication, translation gives `6` antipodal cosets and no nonpremise theorem maps the nadsoliton ontology to the nilradical object.  No named source-atom coupling, damping source, selector/orientation source, action-density installation, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2977/S1927 nilradical provenance `L_total` guard", "## P2977/S1927 nilradical provenance `L_total` guard\n\n`P2977/S1927` does not add a nilradical term to `L_total`.  Ring canonicity and unit-fixity of `{0,6}`/`6` remain finite provenance-readiness only; no strict nadsoliton source map, nonpremise zero-section source, named coupling, action-density theorem, EOM, Hamiltonian, bridge closure, role transfer, or ToE is exported.\n")
    append_once(AGENTS, "Current nilradical strict-nadsoliton provenance obstruction guardrail (P2977/S1927, 2026-06-20)", "## Current nilradical strict-nadsoliton provenance obstruction guardrail (P2977/S1927, 2026-06-20)\n\n- P2977 attacks exactly one P2976 missing theorem: strict nadsoliton provenance for the `Z/12Z` nilradical filtration.\n- Finite positives remain real: `{0,6}` is ring-canonical, `6^2=0`, all `4` units fix the nilradical set, and the acceptance matrix has `64` profiles with only the full profile accepting.\n- The current route is bounded no-go: no strict source exports the ring zero/multiplication from the nadsoliton, translation leaves `6` antipodal cosets without a nonpremise localizer, and no named source-atom coupling is supplied.\n- Do not promote nilradical canonicity, unit-fixity of `6`, translated-coset scans, or formal source placeholders to selector closure, damping source, bridge completion, role transfer, nonproxy `L_total`, or ToE.  A next admissible move may attack one remaining nilradical route (named source-atom coupling or action/variational installation) or pivot to a genuinely new typed object while preserving the P2929-P2977 no-strict-export boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
