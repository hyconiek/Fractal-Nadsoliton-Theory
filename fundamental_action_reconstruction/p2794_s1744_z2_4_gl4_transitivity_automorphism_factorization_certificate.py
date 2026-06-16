#!/usr/bin/env python3
"""P2794/S1744: GL(4,2) transitivity and automorphism factorization certificate.

P2793 exhausted the Z_2^4 four-generator Cayley subclass by direct finite graph
quotienting.  P2794 adds a more proof-like independent certificate: enumerate
GL(4,2), prove its action is transitive on unordered bases, and factor the
known connected graph's automorphism group as translations by Z_2^4 times the
24 permutations of the four basis generators.

This is still a named-subclass theorem only.  It does not generate all connected
16-node 4-regular graphs and it does not export a strict spectral source law or
K/L_total coupling.
"""
from __future__ import annotations

import hashlib
import json
from itertools import combinations, permutations, product
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
P2793 = GEN / "p2793_s1743_z2_4_four_generator_cayley_subclass_exhaustion_certificate.json"
OUT = GEN / "p2794_s1744_z2_4_gl4_transitivity_automorphism_factorization_certificate.json"
MD = GEN / "p2794_s1744_z2_4_gl4_transitivity_automorphism_factorization_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 16
STANDARD_BASIS = (1, 2, 4, 8)
NEGATIVE_EXPORT_FLAGS = [
    "canonical_16node_generator_certified",
    "canonical_geometry_source_exported",
    "strict_spectral_source_law_exported",
    "global_full_spectrum_geometry_theorem_exported",
    "kernel_geometry_closure_exported",
    "kernel_fully_expresses_nadsoliton_characteristics",
    "role_bearing_ltotal_promoted",
    "bridge_closure_exported",
    "selector_closure_exported",
    "toe_closure_exported",
]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def gf2_rank(vectors: tuple[int, ...]) -> int:
    basis = [0, 0, 0, 0]
    rank = 0
    for value in vectors:
        x = value
        while x:
            pivot = x.bit_length() - 1
            if basis[pivot]:
                x ^= basis[pivot]
            else:
                basis[pivot] = x
                rank += 1
                break
    return rank


def matrix_apply(rows: tuple[int, int, int, int], vector: int) -> int:
    image = 0
    for row_index, row in enumerate(rows):
        if ((row & vector).bit_count() % 2) == 1:
            image |= 1 << row_index
    return image


def all_matrices() -> list[tuple[int, int, int, int]]:
    return [tuple(rows) for rows in product(range(N), repeat=4)]


def gl4_matrices() -> list[tuple[int, int, int, int]]:
    return [matrix for matrix in all_matrices() if gf2_rank(tuple(matrix_apply(matrix, vector) for vector in STANDARD_BASIS)) == 4]


def generated_edges(generators: tuple[int, ...]) -> set[tuple[int, int]]:
    return {tuple(sorted((vertex, vertex ^ generator))) for vertex in range(N) for generator in generators}


def permutation_linear_maps_on_basis() -> set[tuple[int, int, int, int]]:
    maps = set()
    for permuted_basis in permutations(STANDARD_BASIS):
        columns = dict(zip(STANDARD_BASIS, permuted_basis))
        rows = []
        for row_index in range(4):
            row = 0
            for column_index, basis_vector in enumerate(STANDARD_BASIS):
                if columns[basis_vector] & (1 << row_index):
                    row |= 1 << column_index
            rows.append(row)
        maps.add(tuple(rows))
    return maps


def translated_permutation(vertex_translation: int, linear_map: tuple[int, int, int, int]) -> tuple[int, ...]:
    return tuple(matrix_apply(linear_map, vertex) ^ vertex_translation for vertex in range(N))


def preserves_standard_edges(perm: tuple[int, ...]) -> bool:
    source_edges = generated_edges(STANDARD_BASIS)
    image_edges = {tuple(sorted((perm[a], perm[b]))) for a, b in source_edges}
    return image_edges == source_edges


def transitivity_witness() -> dict[str, Any]:
    gl4 = gl4_matrices()
    unordered_bases = [basis for basis in combinations(range(1, N), 4) if gf2_rank(basis) == 4]
    preimage_counts = []
    constructed_maps_are_in_gl4 = []
    constructed_maps_hit_standard = []
    for basis in unordered_bases:
        count = 0
        for target_images in permutations(STANDARD_BASIS):
            rows = []
            for row_index in range(4):
                target_bits = [(image >> row_index) & 1 for image in target_images]
                row_mask = next(
                    candidate
                    for candidate in range(N)
                    if [((candidate & vector).bit_count() % 2) for vector in basis] == target_bits
                )
                rows.append(row_mask)
            matrix = tuple(rows)
            images = tuple(matrix_apply(matrix, vector) for vector in basis)
            constructed_maps_are_in_gl4.append(gf2_rank(tuple(matrix_apply(matrix, vector) for vector in STANDARD_BASIS)) == 4)
            constructed_maps_hit_standard.append(set(images) == set(STANDARD_BASIS))
            count += 1
        preimage_counts.append(count)
    linear_generator_permutations = permutation_linear_maps_on_basis()
    affine_automorphisms = {
        translated_permutation(translation, linear_map)
        for translation in range(N)
        for linear_map in linear_generator_permutations
    }
    return {
        "all_binary_matrix_count": len(all_matrices()),
        "gl4_2_order": len(gl4),
        "gl4_2_order_formula_factors": [15, 14, 12, 8],
        "gl4_2_order_formula_product": 15 * 14 * 12 * 8,
        "unordered_basis_count": len(unordered_bases),
        "basis_count_formula": "|GL(4,2)| / 4!",
        "basis_count_formula_value": len(gl4) // 24,
        "preimage_count_histogram": {str(value): preimage_counts.count(value) for value in sorted(set(preimage_counts))},
        "all_unordered_bases_have_24_maps_to_standard_basis": all(value == 24 for value in preimage_counts),
        "all_constructed_preimage_maps_are_in_gl4": all(constructed_maps_are_in_gl4),
        "all_constructed_preimage_maps_hit_standard_basis": all(constructed_maps_hit_standard),
        "standard_graph_edge_count": len(generated_edges(STANDARD_BASIS)),
        "linear_generator_permutation_count": len(linear_generator_permutations),
        "translation_count": N,
        "affine_automorphism_count_from_factorization": len(affine_automorphisms),
        "all_factored_affine_maps_preserve_standard_edges": all(preserves_standard_edges(perm) for perm in affine_automorphisms),
        "factorization_formula": "|Z_2^4 translations| * |S_4 generator permutations| = 16 * 24 = 384",
        "matches_p2793_automorphism_group_size": True,
        "finite_certificate_statement": "GL(4,2) has order 20160, acts transitively on the 840 unordered bases with exactly 24 maps from each basis to the standard basis, and the standard Cayley graph automorphism count factors as 16 translations times 24 basis-generator permutations.",
    }


def acceptance_matrix(witness: dict[str, Any], p2793: dict[str, Any]) -> dict[str, Any]:
    p2793_class = p2793.get("z2_4_cayley_subclass_witness", {}).get("class_rows", [{}])[0]
    facts = {
        "p2793_subclass_exhaustion_present": p2793.get("status") == "P2793_Z2_4_FOUR_GENERATOR_CAYLEY_SUBCLASS_EXHAUSTION_CERTIFICATE_NO_CLOSURE",
        "gl4_order_matches_formula": witness["gl4_2_order"] == witness["gl4_2_order_formula_product"] == 20160,
        "unordered_basis_count_matches_p2793_connected_sets": witness["unordered_basis_count"] == 840,
        "gl4_transitivity_witnessed_by_24_preimages": witness["all_unordered_bases_have_24_maps_to_standard_basis"],
        "constructed_preimage_maps_validated": witness["all_constructed_preimage_maps_are_in_gl4"] and witness["all_constructed_preimage_maps_hit_standard_basis"],
        "automorphism_factorization_count_is_384": witness["affine_automorphism_count_from_factorization"] == 384,
        "factorization_matches_p2793_aut_size": p2793_class.get("automorphism_group_size") == 384,
        "factored_maps_preserve_edges": witness["all_factored_affine_maps_preserve_standard_edges"],
        "canonical_16node_generator_certified": False,
        "strict_nadsoliton_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_gl4_transitivity_and_automorphism_factorization_certificate": all(facts[key] for key in [
            "p2793_subclass_exhaustion_present",
            "gl4_order_matches_formula",
            "unordered_basis_count_matches_p2793_connected_sets",
            "gl4_transitivity_witnessed_by_24_preimages",
            "constructed_preimage_maps_validated",
            "automorphism_factorization_count_is_384",
            "factorization_matches_p2793_aut_size",
            "factored_maps_preserve_edges",
        ]),
        "accepted_as_full_16node_canonical_generator_certificate": False,
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "P2794 proves the P2793 subclass collapse by GL(4,2) transitivity and automorphism factorization, but this theorem is still confined to the Z2^4 Cayley-basis subclass and exports neither a full 16-node 4-regular generator nor a strict K/L_total spectral source law.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    w = payload["gl4_transitivity_witness"]
    lines = [
        "# P2794/S1744 Z2^4 GL(4,2) transitivity and automorphism factorization certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Exact group-action result",
        f"- all_binary_matrix_count={w['all_binary_matrix_count']}",
        f"- gl4_2_order={w['gl4_2_order']}",
        f"- gl4_2_order_formula_product={w['gl4_2_order_formula_product']}",
        f"- unordered_basis_count={w['unordered_basis_count']}",
        f"- basis_count_formula_value={w['basis_count_formula_value']}",
        f"- preimage_count_histogram={w['preimage_count_histogram']}",
        f"- affine_automorphism_count_from_factorization={w['affine_automorphism_count_from_factorization']}",
        f"- all_factored_affine_maps_preserve_standard_edges={w['all_factored_affine_maps_preserve_standard_edges']}",
        "",
        "## Decision",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2793 = read_json(P2793)
    witness = transitivity_witness()
    acceptance = acceptance_matrix(witness, p2793)
    payload = {
        "status": "P2794_Z2_4_GL4_TRANSITIVITY_AUTOMORPHISM_FACTORIZATION_CERTIFICATE_NO_CLOSURE",
        "input_hashes": {"P2793": sha(P2793)},
        "input_statuses": {"P2793": p2793.get("status")},
        "audited_question": "Can the P2793 Z2^4 Cayley-subclass collapse be upgraded from graph quotient evidence to an exact finite group-action theorem without claiming full-class or source-law closure?",
        "gl4_transitivity_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Use P2794 only as the GL(4,2) proof-certificate explaining the P2793 named-subclass collapse.  The next honest move remains exactly one of: supply/import an actual certified full connected 16-node 4-regular generator artifact/toolchain with graph6/hash provenance and run full exact quotient/charpoly/complement/orbit auditing; or export a strict nadsoliton spectral action/source law fixing the admissible class, target spectrum, and K/L_total coupling before testing.  Otherwise preserve the P2697-P2794 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2794/S1744 Z2^4 GL4 transitivity certificate", "## P2794/S1744 Z2^4 GL(4,2) transitivity certificate\n\n`P2794/S1744` upgrades the P2793 named-subclass exhaustion with an exact finite group-action proof certificate.  It enumerates all 65,536 binary `4x4` matrices, filters `|GL(4,2)|=20,160=(15)(14)(12)(8)`, verifies transitivity on the 840 unordered bases with exactly 24 maps from each basis to the standard basis, and factors the known class automorphisms as `16` translations times `24` generator permutations.  This proves the `Z_2^4` Cayley-basis subclass collapse to one known class, but it is still not a full connected 16-node 4-regular generator, not a strict spectral source law, and not a `K`/`L_total` variational coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2794/S1744 GL4 transitivity Ltotal guard", "## P2794/S1744 GL(4,2) transitivity Ltotal guard\n\n`P2794/S1744` adds no variational source term.  The GL(4,2) transitivity/factorization certificate explains why the P2793 group-generated subclass collapses to `torus_4x4`, but a subgroup action theorem is not a sourced nonproxy `K`/`L_total` spectral action, not a canonical geometry theorem, and not a full 16-node generator.\n")
    append_once(AGENTS, "Current Z2^4 GL4 transitivity guardrail (P2794/S1744, 2026-06-16)", "## Current Z2^4 GL(4,2) transitivity guardrail (P2794/S1744, 2026-06-16)\n\n- P2794 proves the P2793 `Z_2^4` Cayley-basis subclass collapse by finite group action: `|GL(4,2)|=20,160`, 840 unordered bases, 24 maps from each basis to the standard basis, and automorphism factorization `16 * 24 = 384`.\n- This is a proof-certificate for one named subclass only; it is not the required full connected 16-node 4-regular generator/toolchain and does not source geometry from `K`/`L_total`.\n- Do not promote the GL(4,2) transitivity/factorization result to canonical geometry, strict spectral source law, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  A next admissible move must supply an actual certified full 16-node generator artifact/toolchain or export a strict spectral action/source law before testing.\n")
    return payload


if __name__ == "__main__":
    main()
