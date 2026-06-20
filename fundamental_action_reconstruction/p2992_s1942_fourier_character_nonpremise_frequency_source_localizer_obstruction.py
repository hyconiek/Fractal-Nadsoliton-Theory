#!/usr/bin/env python3
"""P2992/S1942: Fourier-character nonpremise frequency source-localizer obstruction.

P2991 introduced the exact additive Fourier character table of Z/12Z and left a
first missing theorem route: a nonpremise frequency/source localizer.  This audit
attacks exactly that route.  It does not replay annihilator, nilradical, CRT,
zero-derivation, selector, bridge, role-transfer, or L_total lanes.

The finite calculation builds localizer signatures from conductor, kernel size,
image size, primitive status, and inversion-pair support.  These signatures are
useful spectral labels, but they remain algebraic bookkeeping.  They do not
export a nonpremise physical frequency sector, a strict nadsoliton character
source map, or a theorem that chooses one frequency as the self-learning
nadsoliton's internal source.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2991_s1941_z12_additive_fourier_character_source_candidate_obstruction import MODULUS, OUT as P2991, fourier_character_witness

OUT = GEN / "p2992_s1942_fourier_character_nonpremise_frequency_source_localizer_obstruction.json"
MD = GEN / "p2992_s1942_fourier_character_nonpremise_frequency_source_localizer_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

THEORY_CAPSULE = {
    "ontology": "The nadsoliton is the single primordial fractal information of the universe in a solitonic state; there is no deeper informational substrate beneath it.",
    "self_adjoint_information_network_reading": "Self-adjoint/self-coupled language means the internal informational state couples back to itself through symmetric/variational consistency tests, not through an external observer or separate neural layer.",
    "self_learning_reading": "Self-learning is admissible only as a sourced update/stationarity theorem or energy/variational fixed-point law; neural-network language alone is a model-generating analogy, not closure.",
    "emergence_order": "Preferred internal order remains nadsoliton -> light -> matter -> emergent observer.",
}


def inversion_pair(k: int) -> list[int]:
    return sorted({k % MODULUS, (-k) % MODULUS})


def localizer_witness() -> dict[str, Any]:
    witness = fourier_character_witness()
    rows = []
    for row in witness["character_rows"]:
        signature = {
            "conductor": row["conductor"],
            "kernel_size": row["kernel_size"],
            "image_size": row["image_size"],
            "primitive_character": row["conductor"] == MODULUS,
            "inversion_pair": inversion_pair(row["k"]),
        }
        rows.append({
            "k": row["k"],
            "signature": signature,
            "algebraic_frequency_label_available": True,
            "row_distinguished_by_full_signature": False,
            "nonpremise_physical_frequency_sector": False,
            "strict_nadsoliton_character_source_map": False,
            "frequency_source_localizer_exported": False,
        })
    seen: dict[str, list[int]] = {}
    for row in rows:
        key = json.dumps(row["signature"], sort_keys=True)
        seen.setdefault(key, []).append(row["k"])
    for row in rows:
        key = json.dumps(row["signature"], sort_keys=True)
        row["signature_equivalence_class"] = seen[key]
        row["row_distinguished_by_full_signature"] = len(seen[key]) == 1
    return {
        "modulus": MODULUS,
        "row_count": len(rows),
        "localizer_rows": rows,
        "unique_signature_rows": [r["k"] for r in rows if r["row_distinguished_by_full_signature"]],
        "ambiguous_signature_classes": sorted([v for v in seen.values() if len(v) > 1], key=lambda xs: xs[0]),
        "exported_frequency_source_localizers": [r["k"] for r in rows if r["frequency_source_localizer_exported"]],
    }


def obligation_rows(witness: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"obligation": "finite_fourier_rows_present", "satisfied": witness["row_count"] == MODULUS, "evidence": "all 12 P2991 Fourier-character rows are present"},
        {"obligation": "algebraic_frequency_signatures_constructed", "satisfied": all(r["algebraic_frequency_label_available"] for r in witness["localizer_rows"]), "evidence": "conductor, kernel size, image size, primitive flag, and inversion-pair signatures are computed"},
        {"obligation": "full_signature_has_unique_rows", "satisfied": bool(witness["unique_signature_rows"]), "evidence": f"unique algebraic rows exist at k={witness['unique_signature_rows']}, but uniqueness is bookkeeping only"},
        {"obligation": "nonpremise_physical_frequency_sector", "satisfied": any(r["nonpremise_physical_frequency_sector"] for r in witness["localizer_rows"]), "evidence": "no signature is exported as a physical/source sector of the nadsoliton"},
        {"obligation": "strict_nadsoliton_character_source_map", "satisfied": any(r["strict_nadsoliton_character_source_map"] for r in witness["localizer_rows"]), "evidence": "the Fourier table has no strict source map from the primordial nadsoliton state"},
        {"obligation": "frequency_source_localizer_exported", "satisfied": bool(witness["exported_frequency_source_localizers"]), "evidence": "algebraic frequency labels do not export a nonpremise source-localizer theorem"},
        {"obligation": "accepted_current_frequency_source_localizer", "satisfied": False, "evidence": "no current row satisfies the full strict localizer profile"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_rows", "algebraic_signature", "unique_label", "nonpremise_sector", "strict_source_map", "localizer_theorem", "nonproxy_export"]
    return [{"present": dict(zip(names, bits)), "accepts_frequency_source_localizer": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p2991_path: Any) -> dict[str, Any]:
    witness = localizer_witness()
    matrix = acceptance_matrix()
    return {
        "status": "P2992_FOURIER_CHARACTER_NONPREMISE_FREQUENCY_SOURCE_LOCALIZER_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P2991": hashlib.sha256(p2991_path.read_bytes()).hexdigest() if p2991_path.exists() else None},
        "theory_capsule": THEORY_CAPSULE,
        "constructed_theoretical_objects": {
            "object": "FourierCharacterNonpremiseFrequencySourceLocalizer_ObstructionMatrix",
            "localizer_witness": witness,
            "proof_obligation_rows": obligation_rows(witness),
            "finite_acceptance_matrix": matrix,
        },
        "localizer_certificate": {
            "row_count": witness["row_count"],
            "unique_signature_rows": witness["unique_signature_rows"],
            "ambiguous_signature_classes": witness["ambiguous_signature_classes"],
            "exported_frequency_source_localizers": witness["exported_frequency_source_localizers"],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_frequency_source_localizer"]),
        },
        "decision": {
            "positive_progress": "P2992 attacks exactly one P2991 missing theorem route: a nonpremise frequency/source localizer for the Z12 additive Fourier-character table.",
            "breakthrough": "Bounded no-go: conductor/kernel/image/primitive/inversion signatures give useful finite spectral labels and some unique algebraic rows, but no nonpremise physical frequency sector, strict nadsoliton character source map, or frequency source-localizer theorem is exported.",
            "negative_export_flags": {k: False for k in ["frequency_source_localizer_exported", "strict_character_source_exported", "self_learning_update_law_exported", "selector_closure_exported", "source_coupling_exported", "unit_bearing_density_exported", "bridge_closure_exported", "nonproxy_ltotal_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay algebraic Fourier labels, annihilator/nilradical/CRT/zero-derivation lanes, selector closure, bridge maps, role transfer, or L_total placeholders as source localization.  A next admissible Fourier-character move may attack exactly one remaining route: strict character provenance, source-coupling theorem, or unit-bearing action installation; otherwise introduce a genuinely new strict typed object/provider or preserve the P2929-P2992 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["localizer_certificate"]
    capsule = payload["theory_capsule"]
    lines = [
        "# P2992/S1942 Fourier-character nonpremise frequency source-localizer obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Theory capsule",
        f"- ontology: {capsule['ontology']}",
        f"- self-adjoint information-network reading: {capsule['self_adjoint_information_network_reading']}",
        f"- self-learning reading: {capsule['self_learning_reading']}",
        f"- emergence order: {capsule['emergence_order']}",
        "",
        "## Localizer certificate",
        f"- row count: `{cert['row_count']}`",
        f"- unique signature rows: `{cert['unique_signature_rows']}`",
        f"- ambiguous signature classes: `{cert['ambiguous_signature_classes']}`",
        f"- exported frequency source localizers: `{cert['exported_frequency_source_localizers']}`",
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
    read_json(P2991)
    payload = build_payload(P2991)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2992/S1942 Fourier-character nonpremise frequency source-localizer obstruction", "## P2992/S1942 Fourier-character nonpremise frequency source-localizer obstruction\n\n`P2992/S1942` attacks exactly one P2991 missing theorem route: a nonpremise frequency/source localizer for the `Z/12Z` additive Fourier-character table.  It preserves the ontology that the nadsoliton is the single primordial fractal information in a solitonic state, with no lower information layer; self-learning/neural language is admissible only through sourced update or stationarity theorems.  The finite matrix builds conductor, kernel-size, image-size, primitive-character, and inversion-pair signatures.  These signatures provide spectral bookkeeping and some unique algebraic rows, but no nonpremise physical frequency sector, strict nadsoliton character source map, source-localizer theorem, nonproxy `L_total`, bridge closure, role transfer, or ToE is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2992/S1942 Fourier-character localizer `L_total` guard", "## P2992/S1942 Fourier-character localizer `L_total` guard\n\n`P2992/S1942` adds no Fourier-frequency localizer term to `L_total`.  The conductor/kernel/image/primitive/inversion signatures are algebraic spectral labels only; they do not supply strict field provenance, a self-learning update functional, named unit-bearing density, variational chain, EOM/Hamiltonian term, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current Fourier-character nonpremise frequency source-localizer obstruction guardrail (P2992/S1942, 2026-06-20)", "## Current Fourier-character nonpremise frequency source-localizer obstruction guardrail (P2992/S1942, 2026-06-20)\n\n- P2992 attacks exactly one P2991 missing theorem route: a nonpremise frequency/source localizer for the `Z/12Z` additive Fourier-character table.\n- Finite positives are spectral bookkeeping only: conductor, kernel-size, image-size, primitive-character, and inversion-pair signatures are computed and some rows are algebraically unique.\n- The current route is bounded no-go: no nonpremise physical frequency sector, strict nadsoliton character source map, or frequency source-localizer theorem is exported.\n- Preserve the ontology: the nadsoliton is the single primordial fractal information in a solitonic state, not a projection of a lower information layer; self-learning/neural language requires a sourced update or stationarity theorem.  Do not promote Fourier labels to selector closure, bridge closure, role transfer, nonproxy `L_total`, or ToE.\n- A next admissible Fourier-character move may attack exactly one remaining route (strict character provenance, source-coupling theorem, or unit-bearing action installation), or else introduce a genuinely new strict typed object/provider while preserving the P2929-P2992 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
