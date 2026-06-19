#!/usr/bin/env python3
"""P2897/S1847: circulant relation basepoint obstruction.

P2896 exhausted source-neutral invariant scalar scores on the P2895 free Z12
torsor.  P2897 pivots to a genuinely different finite typed object: arbitrary
translation-invariant binary relations on the same torsor, equivalently directed
circulant relation structures indexed by subsets of Z12 differences.

The computation exhausts all 2^12 such relations.  Every relation has the full
translation group as automorphisms, so every vertex remains in one orbit of size
12; row/degree profiles are translated copies and never isolate a unique
basepoint.  Relation structure can provide invariant geometry, but not a
nonimported source/origin without a translation-breaking law.
"""
from __future__ import annotations

import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

N = 12
P2896 = GEN / "p2896_s1846_invariant_scalar_basepoint_law_exhaustion.json"
OUT = GEN / "p2897_s1847_circulant_relation_basepoint_obstruction.json"
MD = GEN / "p2897_s1847_circulant_relation_basepoint_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def translate(point: int, shift: int) -> int:
    return (point + shift) % N


def difference_subset(mask: int) -> set[int]:
    return {d for d in range(N) if mask & (1 << d)}


def relation_row(diff_set: set[int], source: int) -> tuple[int, ...]:
    return tuple(1 if (target - source) % N in diff_set else 0 for target in range(N))


def is_translation_automorphism(diff_set: set[int], shift: int) -> bool:
    for source in range(N):
        for target in range(N):
            before = (target - source) % N in diff_set
            after = (translate(target, shift) - translate(source, shift)) % N in diff_set
            if before != after:
                return False
    return True


def orbit_under_translations(point: int) -> list[int]:
    return sorted({translate(point, shift) for shift in range(N)})


def relation_summary(mask: int) -> dict[str, Any]:
    diff_set = difference_subset(mask)
    rows = [relation_row(diff_set, source) for source in range(N)]
    out_degrees = [sum(row) for row in rows]
    in_degrees = [sum(1 for source in range(N) if (target - source) % N in diff_set) for target in range(N)]
    local_rows = [tuple(1 if d in diff_set else 0 for d in range(N)) for _ in range(N)]
    local_profile_multiplicities: dict[str, int] = {}
    for row in local_rows:
        key = "".join(str(bit) for bit in row)
        local_profile_multiplicities[key] = local_profile_multiplicities.get(key, 0) + 1
    translations_are_automorphisms = all(is_translation_automorphism(diff_set, shift) for shift in range(N))
    vertex_orbit = orbit_under_translations(0)
    return {
        "mask": mask,
        "difference_subset": sorted(diff_set),
        "difference_subset_size": len(diff_set),
        "translations_are_automorphisms": translations_are_automorphisms,
        "translation_vertex_orbit_size": len(vertex_orbit),
        "out_degree_values": sorted(set(out_degrees)),
        "in_degree_values": sorted(set(in_degrees)),
        "constant_in_out_degree": len(set(out_degrees)) == 1 and len(set(in_degrees)) == 1,
        "unique_degree_vertex_count": sum(1 for v in range(N) if out_degrees.count(out_degrees[v]) == 1 and in_degrees.count(in_degrees[v]) == 1),
        "unique_local_row_profile_count": sum(1 for row in local_rows if local_profile_multiplicities["".join(str(bit) for bit in row)] == 1),
        "can_select_unique_vertex_by_translation_invariant_relation": translations_are_automorphisms and len(vertex_orbit) == 1,
    }


def circulant_relation_exhaustion() -> dict[str, Any]:
    summaries = [relation_summary(mask) for mask in range(1 << N)]
    undirected_loopless = []
    for row in summaries:
        ds = set(row["difference_subset"])
        if 0 in ds:
            continue
        if all(((-d) % N) in ds for d in ds):
            undirected_loopless.append(row)
    return {
        "torsor_size": N,
        "directed_circulant_relation_count": len(summaries),
        "undirected_loopless_circulant_graph_count": len(undirected_loopless),
        "all_relations_have_translation_automorphisms": all(row["translations_are_automorphisms"] for row in summaries),
        "translation_orbit_size_histogram": {str(size): sum(1 for row in summaries if row["translation_vertex_orbit_size"] == size) for size in sorted({row["translation_vertex_orbit_size"] for row in summaries})},
        "relations_selecting_unique_vertex_count": sum(1 for row in summaries if row["can_select_unique_vertex_by_translation_invariant_relation"]),
        "relations_with_unique_degree_vertex_count": sum(1 for row in summaries if row["unique_degree_vertex_count"] > 0),
        "relations_with_unique_local_row_profile_count": sum(1 for row in summaries if row["unique_local_row_profile_count"] > 0),
        "sample_relations": [summaries[i] for i in (0, 1, 3, 5, 2730, 4095)],
        "proof_certificate": {
            "circulant_relation_fact": "A translation-invariant binary relation on a free Z12 torsor is determined by a subset D of Z12 differences.",
            "automorphism_fact": "For every such D, every translation t -> t+s preserves all differences, so Z12 acts by relation automorphisms.",
            "selection_obstruction": "Since the translation automorphism group is transitive on the 12 vertices, no relation-internal invariant can define one vertex without an added translation-breaking source.",
        },
    }


def build_payload(p2896: dict[str, Any]) -> dict[str, Any]:
    audit = circulant_relation_exhaustion()
    facts = {
        "p2896_rechecked": p2896.get("status") == "P2896_INVARIANT_SCALAR_BASEPOINT_LAW_EXHAUSTION_NO_CLOSURE",
        "relation_object_is_new_after_scalar_scores": True,
        "all_relations_translation_symmetric": audit["all_relations_have_translation_automorphisms"],
        "all_vertex_orbits_size_12": audit["translation_orbit_size_histogram"] == {"12": 1 << N},
        "no_relation_selects_unique_vertex": audit["relations_selecting_unique_vertex_count"] == 0,
        "no_degree_or_local_row_profile_basepoint": audit["relations_with_unique_degree_vertex_count"] == 0 and audit["relations_with_unique_local_row_profile_count"] == 0,
    }
    return {
        "status": "P2897_CIRCULANT_RELATION_BASEPOINT_OBSTRUCTION_NO_CLOSURE",
        "input_hashes": {"P2896": sha(P2896)},
        "circulant_relation_basepoint_obstruction": {
            "input_status_rechecked": p2896.get("status"),
            "candidate_class": "arbitrary translation-invariant binary/circulant relation structures on the free Z12 torsor, tested as non-scalar internal basepoint localizers",
            "relation_exhaustion": audit,
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_nonimported_basepoint_law": False,
            "exports_translation_breaking_relation_source": False,
            "exports_unique_coupling_to_9_over_5_carrier": False,
            "exports_nonimported_9_over_5_variational_density": False,
        },
        "decision": {
            "negative_export_flags": {
                "nonimported_basepoint_or_polarity_law_exported": False,
                "translation_breaking_relation_source_exported": False,
                "strict_free_12_torsor_source_law_exported": False,
                "strict_phase_origin_source_artifact_exported": False,
                "unique_coupling_to_9_over_5_carrier_exported": False,
                "coupling_to_9_over_5_variational_density_exported": False,
                "nonimported_9_over_5_variational_chain_rule_exported": False,
                "localized_action_density_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2897 exhausts all 4096 translation-invariant binary relations on the P2895/P2896 free Z12 torsor.  Every relation is a circulant difference relation preserved by all 12 translations, so every vertex remains in one orbit of size 12.  Degree profiles and translation-local row profiles never isolate a unique vertex.  Therefore relation geometry does not create a nonimported basepoint/polarity law or source the 9/5 carrier origin.",
            "next_honest_step": "Do not promote circulant relations, distance/adjacency profiles, invariant relation geometry, graph row profiles, degree profiles, scalar scores, canonical zeros, or unpointed free-torsor clocks to strict phase/origin sourcehood.  A next proof-grade move must either supply an explicit non-circulant strict translation-breaking source law with computed basepoint/polarity and coupling theorem to the 9/5 variational density, or pivot to a genuinely different typed object outside torsor/basepoint/scalar-score/relation/support/orbit/Fourier/inventory constructions; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["circulant_relation_basepoint_obstruction"]["relation_exhaustion"]
    lines = [
        "# P2897/S1847 circulant relation basepoint obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite relation audit",
        f"- torsor size: `{audit['torsor_size']}`",
        f"- directed circulant relation count: `{audit['directed_circulant_relation_count']}`",
        f"- undirected loopless circulant graph count: `{audit['undirected_loopless_circulant_graph_count']}`",
        f"- all relations have translation automorphisms: `{audit['all_relations_have_translation_automorphisms']}`",
        f"- translation orbit size histogram: `{audit['translation_orbit_size_histogram']}`",
        f"- relations selecting a unique vertex: `{audit['relations_selecting_unique_vertex_count']}`",
        f"- relations with a unique degree-profile vertex: `{audit['relations_with_unique_degree_vertex_count']}`",
        f"- relations with a unique translation-local row-profile vertex: `{audit['relations_with_unique_local_row_profile_count']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2896))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2897/S1847 circulant relation basepoint obstruction",
        "## P2897/S1847 circulant relation basepoint obstruction\n\n"
        "`P2897/S1847` exhausts all `4096` translation-invariant binary/circulant relations on the `P2895/P2896` free `Z12` torsor.  Every such relation is preserved by all translations, the vertex orbit size remains `12`, and neither degree profiles nor translation-local relation-row profiles isolate a unique vertex.  Thus invariant relation geometry does not export a nonimported basepoint/polarity law, `9/5` variational density, localized action density, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2897/S1847 circulant relation basepoint `L_total` guard",
        "## P2897/S1847 circulant relation basepoint `L_total` guard\n\n"
        "`P2897/S1847` is a finite relation-geometry obstruction, not a strict action construction.  It adds no non-circulant translation-breaking source, no localized unit-bearing density, no coupling theorem to the `9/5` carrier, and no variational chain rule into nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current circulant relation basepoint obstruction guardrail (P2897/S1847, 2026-06-19)",
        "## Current circulant relation basepoint obstruction guardrail (P2897/S1847, 2026-06-19)\n\n"
        "- P2897 exhausts all `4096` translation-invariant binary/circulant relations on the P2895/P2896 free `Z12` torsor as a non-scalar internal basepoint-localizer candidate.\n"
        "- Every relation is preserved by all `12` translations, every vertex remains in one orbit of size `12`, and neither degree profiles nor relation-row profiles isolate a unique basepoint.\n"
        "- Do not promote circulant relations, distance/adjacency profiles, invariant relation geometry, graph row profiles, degree profiles, scalar scores, canonical zero choices, unpointed free-torsor clocks, support/orbit/Fourier data, or inventory hits to strict phase/origin sourcehood, strict damping/compression bridge, selector closure, role transfer, `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must supply an explicit non-circulant strict translation-breaking source law with computed basepoint/polarity and coupling to the `9/5` variational density, pivot to a genuinely different typed object outside torsor/basepoint/scalar-score/relation/support/orbit/Fourier/inventory constructions, or preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
