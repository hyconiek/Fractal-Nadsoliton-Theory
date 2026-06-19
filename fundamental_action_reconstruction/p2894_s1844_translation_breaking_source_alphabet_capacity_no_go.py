#!/usr/bin/env python3
"""P2894/S1844: translation-breaking source alphabet capacity no-go.

P2893 proved that a source-neutral quotient section cannot choose an embedded
representative in the free Z12 orbit of 9/5 carriers.  The next honest question
is how much translation-breaking structure a new source would need.  P2894 tests
the finite Z12-set capacity boundary: no scalar, sign, low-period, or any finite
source alphabet of size < 12 can map equivariantly onto a free 12-origin carrier
orbit.  A successful translation-breaking source must itself contain a free
12-phase torsor, and even then a polarity/basepoint law is still required.
"""
from __future__ import annotations

import json
from functools import lru_cache
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

N = 12
P2893 = GEN / "p2893_s1843_free_orbit_section_obstruction_theorem.json"
OUT = GEN / "p2894_s1844_translation_breaking_source_alphabet_capacity_no_go.json"
MD = GEN / "p2894_s1844_translation_breaking_source_alphabet_capacity_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
ORBIT_SIZES = (1, 2, 3, 4, 6, 12)
LOW_CAPACITY_ORBIT_SIZES = (1, 2, 3, 4, 6)


@lru_cache(maxsize=None)
def orbit_partitions(total: int, max_orbit_size_index: int = len(ORBIT_SIZES) - 1) -> tuple[tuple[int, ...], ...]:
    """Unordered decompositions of a finite Z12-set into transitive orbit sizes."""
    if total == 0:
        return ((),)
    if total < 0 or max_orbit_size_index < 0:
        return ()
    size = ORBIT_SIZES[max_orbit_size_index]
    rows: list[tuple[int, ...]] = []
    for count in range(total // size + 1):
        remainder = total - count * size
        for tail in orbit_partitions(remainder, max_orbit_size_index - 1):
            rows.append(tuple(sorted((size,) * count + tail)))
    return tuple(sorted(set(rows)))


def equivariant_maps_from_orbit_to_free_target(source_orbit_size: int) -> int:
    """Number of Z12-equivariant maps from one transitive orbit to a free orbit.

    A map Z12/H -> Z12 exists iff H is contained in the target stabilizer.  The
    target stabilizer is trivial, so only source_orbit_size=12 works; then the
    image of one point may be any of the 12 target points.
    """
    return N if source_orbit_size == N else 0


def equivariant_maps_from_zset_to_free_target(orbit_sizes: tuple[int, ...]) -> int:
    total = 1
    for size in orbit_sizes:
        total *= equivariant_maps_from_orbit_to_free_target(size)
    return total


def capacity_audit() -> dict[str, Any]:
    low_rows = []
    total_low_zsets = 0
    low_with_map = 0
    for total_size in range(1, N):
        partitions = tuple(part for part in orbit_partitions(total_size) if all(size in LOW_CAPACITY_ORBIT_SIZES for size in part))
        total_low_zsets += len(partitions)
        with_map = [part for part in partitions if equivariant_maps_from_zset_to_free_target(part) > 0]
        low_with_map += len(with_map)
        low_rows.append({"source_alphabet_size": total_size, "z12_set_type_count": len(partitions), "types_with_equivariant_map_to_free_12_orbit": len(with_map), "sample_types": [list(part) for part in partitions[:8]]})
    size_12_partitions = orbit_partitions(N)
    size_12_rows = []
    for part in size_12_partitions:
        map_count = equivariant_maps_from_zset_to_free_target(part)
        if map_count:
            size_12_rows.append({"orbit_decomposition": list(part), "equivariant_map_count": map_count, "boundary": "maps exist only because the source already contains a free 12-torsor; choosing one map still requires basepoint/polarity law"})
    named_sources = [
        {"name": "scalar/trivial", "orbit_decomposition": [1]},
        {"name": "binary sign", "orbit_decomposition": [2]},
        {"name": "ternary phase", "orbit_decomposition": [3]},
        {"name": "quarter phase", "orbit_decomposition": [4]},
        {"name": "half-cycle phase", "orbit_decomposition": [6]},
        {"name": "full Z12 phase torsor", "orbit_decomposition": [12]},
    ]
    for row in named_sources:
        row["equivariant_map_count_to_free_12_orbit"] = equivariant_maps_from_zset_to_free_target(tuple(row["orbit_decomposition"]))
    return {
        "target_orbit_size": N,
        "low_capacity_source_alphabet_sizes_tested": list(range(1, N)),
        "low_capacity_z12_set_type_count": total_low_zsets,
        "low_capacity_types_with_equivariant_map_to_free_12_orbit": low_with_map,
        "low_capacity_rows": low_rows,
        "size_12_z12_set_type_count": len(size_12_partitions),
        "size_12_types_with_map_rows": size_12_rows,
        "named_source_capacity_rows": named_sources,
        "minimum_source_orbit_size_required": N,
        "free_12_torsor_necessary_for_equivariant_targeting": True,
    }


def build_payload(p2893: dict[str, Any]) -> dict[str, Any]:
    audit = capacity_audit()
    facts = {
        "p2893_rechecked": p2893.get("status") == "P2893_FREE_ORBIT_SECTION_OBSTRUCTION_THEOREM_NO_CLOSURE",
        "low_capacity_types_exhausted": audit["low_capacity_z12_set_type_count"] > 0,
        "no_low_capacity_source_maps_to_free_12_orbit": audit["low_capacity_types_with_equivariant_map_to_free_12_orbit"] == 0,
        "free_12_torsor_required": audit["free_12_torsor_necessary_for_equivariant_targeting"],
    }
    return {
        "status": "P2894_TRANSLATION_BREAKING_SOURCE_ALPHABET_CAPACITY_NO_GO_NO_CLOSURE",
        "input_hashes": {"P2893": sha(P2893)},
        "translation_breaking_source_alphabet_capacity_no_go": {
            "input_status_rechecked": p2893.get("status"),
            "candidate_class": "finite Z12-set source alphabets attempting to target a free 12-origin 9/5 carrier orbit equivariantly",
            "capacity_audit": audit,
            "proof_certificate": {
                "orbit_stabilizer_fact": "A Z12-equivariant map from a transitive source orbit Z12/H to the free target orbit Z12 exists only when H is trivial, equivalently when the source orbit has size 12.",
                "finite_result": "All Z12-set decompositions of total size 1..11 have zero equivariant maps to a free 12-origin target orbit.  Scalar, sign, 3/4/6-period phase sources are too small.",
                "free_torsor_boundary": "A full free 12-torsor is necessary but not sufficient for strict sourcehood: it supplies possible equivariant maps only after a basepoint/polarity law fixes which map is sourced.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_translation_breaking_strict_source": False,
            "exports_free_12_torsor_source_law": False,
            "exports_basepoint_or_polarity_law": False,
            "exports_nonimported_9_over_5_variational_density": False,
        },
        "decision": {
            "negative_export_flags": {
                "translation_breaking_strict_source_exported": False,
                "free_12_torsor_source_law_exported": False,
                "basepoint_or_polarity_law_exported": False,
                "strict_phase_origin_source_artifact_exported": False,
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
            "reason": "P2894 tests the finite source-alphabet capacity required after P2893.  Every Z12-set source alphabet of total size below 12, including scalar/sign/low-period phase sources, has zero equivariant maps to a free 12-origin carrier orbit.  A free 12-torsor is necessary, but it is not a strict source law without an exported basepoint/polarity and coupling theorem.",
            "next_honest_step": "Do not continue with scalar, binary sign, 3/4/6-period phase, low-capacity alphabet, canonical representative, or quotient-section conventions as translation-breaking sourcehood.  A next proof-grade move must either construct an explicit strict free-12-torsor source with a nonimported basepoint/polarity law and coupling theorem to the 9/5 variational density, or pivot to a genuinely different typed object outside the quotient-section/source-alphabet/orbit/Fourier/inventory family; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["translation_breaking_source_alphabet_capacity_no_go"]["capacity_audit"]
    lines = [
        "# P2894/S1844 translation-breaking source alphabet capacity no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite capacity audit",
        f"- low-capacity alphabet sizes tested: `{audit['low_capacity_source_alphabet_sizes_tested']}`",
        f"- low-capacity Z12-set types tested: `{audit['low_capacity_z12_set_type_count']}`",
        f"- low-capacity types with equivariant maps to free 12-orbit: `{audit['low_capacity_types_with_equivariant_map_to_free_12_orbit']}`",
        f"- minimum source orbit size required: `{audit['minimum_source_orbit_size_required']}`",
        "",
        "## Named source classes",
    ]
    for row in audit["named_source_capacity_rows"]:
        lines.append(f"- `{row['name']}` {row['orbit_decomposition']} -> maps `{row['equivariant_map_count_to_free_12_orbit']}`")
    lines.extend(["", "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2893))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2894/S1844 translation-breaking source alphabet capacity no-go",
        "## P2894/S1844 translation-breaking source alphabet capacity no-go\n\n"
        "`P2894/S1844` audits the finite source-alphabet capacity needed after `P2893`.  All `Z12`-set source alphabets of total size `1..11` have `0` equivariant maps to a free `12`-origin `9/5` carrier orbit; scalar, binary sign, and `3/4/6`-period phase sources are too small.  A free `12`-torsor is necessary, but it still does not export strict sourcehood without a nonimported basepoint/polarity law and coupling theorem to the `9/5` variational density.  No localized action density, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2894/S1844 translation-breaking source alphabet capacity `L_total` guard",
        "## P2894/S1844 translation-breaking source alphabet capacity `L_total` guard\n\n"
        "`P2894/S1844` is a finite source-capacity obstruction, not a strict action construction.  It provides no free-`12`-torsor source law, no basepoint/polarity law, no localized unit-bearing action density, and no variational chain rule into nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current translation-breaking source alphabet capacity guardrail (P2894/S1844, 2026-06-19)",
        "## Current translation-breaking source alphabet capacity guardrail (P2894/S1844, 2026-06-19)\n\n"
        "- P2894 audits finite `Z12`-set source alphabets after P2893 and proves that all total sizes `1..11` have `0` equivariant maps to a free `12`-origin `9/5` carrier orbit.\n"
        "- Scalar, binary sign, and `3/4/6`-period phase sources are too small; a free `12`-torsor is necessary but still not sufficient without a nonimported basepoint/polarity law and coupling theorem.\n"
        "- Do not promote low-capacity alphabets, scalar/sign/low-period phase sources, free-torsor possibility alone, canonical representatives, quotient-section conventions, Fourier phase/power, or inventory hits to strict phase/origin sourcehood, strict damping/compression bridge, selector closure, role transfer, `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must construct an explicit strict free-`12`-torsor source with nonimported basepoint/polarity and coupling to the `9/5` variational density, pivot to a genuinely different typed object outside quotient-section/source-alphabet/orbit/Fourier/inventory constructions, or preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
