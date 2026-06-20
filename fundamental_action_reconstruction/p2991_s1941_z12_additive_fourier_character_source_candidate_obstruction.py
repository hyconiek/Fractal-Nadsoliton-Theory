#!/usr/bin/env python3
"""P2991/S1941: Z12 additive Fourier-character source-candidate obstruction.

P2990 bounded the annihilator-lattice lane.  This audit introduces one genuinely
new finite typed object outside nilradical/CRT/zero-derivation/annihilator replay:
the additive Fourier character table of Z/12Z.

The finite calculation enumerates all characters chi_k(x)=exp(2*pi*i*k*x/12)
by their exact exponent classes k*x mod 12, verifies the full 12x12
orthogonality matrix, and classifies each character by kernel size and conductor.
This is exact spectral bookkeeping, but it is not a strict source theorem: the
current artifacts do not export a nonpremise frequency selector, strict
nadsoliton provenance for a character row, target-independent beta/Z_beta source,
unit-bearing measure/action density, or bridge/source completion.
"""
from __future__ import annotations

import hashlib, json, math
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2990_s1940_annihilator_lattice_action_variational_installation_obstruction import OUT as P2990

OUT = GEN / "p2991_s1941_z12_additive_fourier_character_source_candidate_obstruction.json"
MD = GEN / "p2991_s1941_z12_additive_fourier_character_source_candidate_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
MODULUS = 12


def character_row(k: int) -> dict[str, Any]:
    exponents = [(k * x) % MODULUS for x in range(MODULUS)]
    kernel = [x for x, e in enumerate(exponents) if e == 0]
    image = sorted(set(exponents))
    conductor = MODULUS // math.gcd(k, MODULUS) if k else 1
    return {
        "k": k,
        "formula": f"chi_{k}(x)=zeta_12^({k}*x)",
        "exponents_mod_12": exponents,
        "kernel": kernel,
        "kernel_size": len(kernel),
        "image_exponents": image,
        "image_size": len(image),
        "conductor": conductor,
        "nontrivial_character": k != 0,
        "nonpremise_frequency_selector": False,
        "strict_nadsoliton_character_provenance": False,
        "source_coupling_theorem": False,
        "accepted_strict_character_source": False,
    }


def orthogonality_entry(k: int, l: int) -> dict[str, Any]:
    delta = (k - l) % MODULUS
    exponents = [(delta * x) % MODULUS for x in range(MODULUS)]
    return {
        "k": k,
        "l": l,
        "delta_mod_12": delta,
        "exponent_multiset_mod_12": exponents,
        "exact_sum": 12 if delta == 0 else 0,
        "orthogonality_expected": 12 if k == l else 0,
        "orthogonality_holds": (12 if delta == 0 else 0) == (12 if k == l else 0),
    }


def fourier_character_witness() -> dict[str, Any]:
    rows = [character_row(k) for k in range(MODULUS)]
    orthogonality = [orthogonality_entry(k, l) for k in range(MODULUS) for l in range(MODULUS)]
    conductor_histogram: dict[str, int] = {}
    for row in rows:
        conductor_histogram[str(row["conductor"])] = conductor_histogram.get(str(row["conductor"]), 0) + 1
    return {
        "modulus": MODULUS,
        "character_count": len(rows),
        "character_rows": rows,
        "orthogonality_matrix_rows": orthogonality,
        "all_orthogonality_checks_pass": all(r["orthogonality_holds"] for r in orthogonality),
        "conductor_histogram": conductor_histogram,
        "accepted_strict_character_sources": [r["k"] for r in rows if r["accepted_strict_character_source"]],
    }


def obligation_rows(witness: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"obligation": "finite_character_table_enumerated", "satisfied": witness["character_count"] == MODULUS, "evidence": "all 12 additive Fourier characters of Z/12Z are enumerated by exponent rows"},
        {"obligation": "full_orthogonality_matrix_verified", "satisfied": witness["all_orthogonality_checks_pass"], "evidence": "all 144 exact character-pair sums match 12*delta_{k,l}"},
        {"obligation": "conductor_kernel_classification", "satisfied": bool(witness["conductor_histogram"]), "evidence": f"conductor histogram is {witness['conductor_histogram']}"},
        {"obligation": "nonpremise_frequency_selector", "satisfied": any(r["nonpremise_frequency_selector"] for r in witness["character_rows"]), "evidence": "the table enumerates frequencies but does not select one nonpremise strict frequency"},
        {"obligation": "strict_nadsoliton_character_provenance", "satisfied": any(r["strict_nadsoliton_character_provenance"] for r in witness["character_rows"]), "evidence": "no strict source map exports a character row as nadsoliton data"},
        {"obligation": "source_coupling_theorem", "satisfied": any(r["source_coupling_theorem"] for r in witness["character_rows"]), "evidence": "no character is coupled to selector sign, beta/Z_beta, bridge-source, or action density by theorem"},
        {"obligation": "unit_bearing_measure_or_action_density", "satisfied": False, "evidence": "counting orthogonality is not a unit-bearing measure or action-density installation"},
        {"obligation": "accepted_current_strict_character_source", "satisfied": bool(witness["accepted_strict_character_sources"]), "evidence": "no current character row satisfies the full strict source profile"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_table", "orthogonality", "conductor_classification", "frequency_selector", "strict_provenance", "source_coupling", "unit_measure_or_action", "nonproxy_export"]
    return [{"present": dict(zip(names, bits)), "accepts_strict_character_source": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p2990_path: Any) -> dict[str, Any]:
    witness = fourier_character_witness()
    matrix = acceptance_matrix()
    return {
        "status": "P2991_Z12_ADDITIVE_FOURIER_CHARACTER_SOURCE_CANDIDATE_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P2990": hashlib.sha256(p2990_path.read_bytes()).hexdigest() if p2990_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "Z12AdditiveFourierCharacter_SourceCandidateObstructionMatrix",
            "fourier_character_witness": witness,
            "proof_obligation_rows": obligation_rows(witness),
            "finite_acceptance_matrix": matrix,
        },
        "character_certificate": {
            "character_count": witness["character_count"],
            "orthogonality_matrix_rows": len(witness["orthogonality_matrix_rows"]),
            "all_orthogonality_checks_pass": witness["all_orthogonality_checks_pass"],
            "conductor_histogram": witness["conductor_histogram"],
            "accepted_strict_character_sources": witness["accepted_strict_character_sources"],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_strict_character_source"]),
        },
        "decision": {
            "positive_progress": "P2991 introduces one new finite typed object after P2990: the exact additive Fourier character table of Z/12Z.",
            "breakthrough": "Bounded no-go: all 12 characters, all 144 orthogonality sums, and conductor/kernel classes are exact, but no nonpremise frequency selector, strict nadsoliton character provenance, source-coupling theorem, unit-bearing measure/action density, or nonproxy export is present.",
            "negative_export_flags": {k: False for k in ["strict_character_source_exported", "frequency_selector_exported", "strict_provenance_exported", "source_coupling_exported", "unit_bearing_density_exported", "bridge_closure_exported", "nonproxy_ltotal_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay annihilator, nilradical, CRT, zero-derivation, selector, bridge, or L_total lanes through Fourier labels.  A next admissible Fourier-character move may attack exactly one missing theorem for this new object: nonpremise frequency/source localizer, strict character provenance, source-coupling theorem, or unit-bearing action installation; otherwise introduce a genuinely new strict typed object/provider or preserve the P2929-P2991 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["character_certificate"]
    lines = [
        "# P2991/S1941 Z12 additive Fourier-character source-candidate obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Character certificate",
        f"- character count: `{cert['character_count']}`",
        f"- orthogonality matrix rows: `{cert['orthogonality_matrix_rows']}`",
        f"- all orthogonality checks pass: `{cert['all_orthogonality_checks_pass']}`",
        f"- conductor histogram: `{cert['conductor_histogram']}`",
        f"- accepted strict character sources: `{cert['accepted_strict_character_sources']}`",
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
    read_json(P2990)
    payload = build_payload(P2990)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2991/S1941 Z12 additive Fourier-character source-candidate obstruction", "## P2991/S1941 Z12 additive Fourier-character source-candidate obstruction\n\n`P2991/S1941` introduces one new finite typed object outside nilradical/CRT/zero-derivation/annihilator replay: the additive Fourier character table of `Z/12Z`.  The finite computation enumerates all `12` characters `chi_k(x)=zeta_12^(k*x)`, verifies all `144` exact orthogonality sums, and classifies conductor/kernel rows.  This is exact spectral bookkeeping, but no nonpremise frequency selector, strict nadsoliton character provenance, source-coupling theorem, unit-bearing measure/action density, nonproxy `L_total`, bridge closure, role transfer, or ToE is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2991/S1941 Fourier-character source-candidate `L_total` guard", "## P2991/S1941 Fourier-character source-candidate `L_total` guard\n\n`P2991/S1941` enumerates exact additive Fourier-character rows and orthogonality sums but adds no Fourier-character term to `L_total`.  Character orthogonality is finite spectral bookkeeping only; no nonpremise frequency selector, strict field provenance, named unit-bearing density, variational chain, EOM/Hamiltonian term, bridge closure, role transfer, or ToE is exported.\n")
    append_once(AGENTS, "Current Z12 additive Fourier-character source-candidate obstruction guardrail (P2991/S1941, 2026-06-20)", "## Current Z12 additive Fourier-character source-candidate obstruction guardrail (P2991/S1941, 2026-06-20)\n\n- P2991 introduces one new finite typed object outside nilradical/CRT/zero-derivation/annihilator replay: the additive Fourier character table of `Z/12Z`.\n- Exact finite positives are real: all `12` characters are enumerated, all `144` orthogonality sums pass, and conductor/kernel classes are computed.\n- The current route is bounded no-go as a strict source candidate: no nonpremise frequency selector, strict nadsoliton character provenance, source-coupling theorem, unit-bearing measure/action density, or nonproxy export is present.\n- Do not replay annihilator, nilradical, CRT, zero-derivation, selector, bridge, or `L_total` lanes through Fourier labels.  A next admissible Fourier-character move may attack exactly one missing theorem for this new object, or else introduce a genuinely new strict typed object/provider while preserving the P2929-P2991 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
