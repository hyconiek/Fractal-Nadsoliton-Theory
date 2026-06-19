#!/usr/bin/env python3
"""P2893/S1843: free-orbit section obstruction theorem for 9/5 carriers.

P2892 honestly preserved no-new-live-frontier after the P2888-P2891 lane.  This
step is not another carrier replay: it extracts the finite group-action theorem
behind the obstruction.  For every free Z12 translation orbit of P2888/P2889
9/5 carriers, a source-neutral quotient-level law cannot choose an embedded
representative because an invariant section would have to choose a point fixed
by the whole translation group.  The finite certificate recomputes all 50 free
orbits and verifies the section obstruction orbit-by-orbit.
"""
from __future__ import annotations

import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2889_s1839_translation_orbit_source_law_9_over_5_no_go_audit import N, P2888, orbit, regenerate_p2888_target_triples

P2889 = GEN / "p2889_s1839_translation_orbit_source_law_9_over_5_no_go_audit.json"
P2892 = GEN / "p2892_s1842_post_phase_origin_inventory_state_map_no_new_live_frontier_certificate.json"
OUT = GEN / "p2893_s1843_free_orbit_section_obstruction_theorem.json"
MD = GEN / "p2893_s1843_free_orbit_section_obstruction_theorem.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def orbit_representatives() -> list[tuple[Any, ...]]:
    unseen = set(regenerate_p2888_target_triples())
    orbits: list[tuple[Any, ...]] = []
    while unseen:
        representative = min(unseen)
        orb = orbit(representative)
        orbits.append(orb)
        unseen.difference_update(orb)
    return sorted(orbits, key=lambda members: min(members))


def stabilizer_size(point: Any) -> int:
    return sum(1 for shift in range(N) if point.rotate(shift) == point)


def quotient_section_audit() -> dict[str, Any]:
    orbits = orbit_representatives()
    rows = []
    invariant_section_count_total = 0
    equivariant_endomorphism_count_total = 0
    for index, members in enumerate(orbits):
        fixed_candidates = [point for point in members if stabilizer_size(point) == N]
        invariant_section_count = len(fixed_candidates)
        # A Z12-equivariant endomorphism of a free transitive orbit exists for each
        # output choice at one basepoint; these are shifts, not quotient-level
        # source-neutral representative sections.
        equivariant_endomorphism_count = len(members)
        invariant_section_count_total += invariant_section_count
        equivariant_endomorphism_count_total += equivariant_endomorphism_count
        if index < 8:
            rows.append(
                {
                    "orbit_index": index,
                    "orbit_size": len(members),
                    "point_stabilizer_size_histogram": {str(size): [stabilizer_size(point) for point in members].count(size) for size in sorted({stabilizer_size(point) for point in members})},
                    "quotient_to_embedded_invariant_section_count": invariant_section_count,
                    "equivariant_endomorphism_count_not_source_neutral": equivariant_endomorphism_count,
                    "sample_representative": min(members).primitive_record(),
                }
            )
    return {
        "target_triple_count": sum(len(members) for members in orbits),
        "translation_orbit_count": len(orbits),
        "orbit_size_histogram": {str(size): [len(members) for members in orbits].count(size) for size in sorted({len(members) for members in orbits})},
        "point_stabilizer_size_histogram_global": {str(size): [stabilizer_size(point) for members in orbits for point in members].count(size) for size in sorted({stabilizer_size(point) for members in orbits for point in members})},
        "quotient_to_embedded_invariant_section_count_total": invariant_section_count_total,
        "equivariant_endomorphism_count_total_not_source_neutral": equivariant_endomorphism_count_total,
        "all_orbits_free": all(len(members) == N for members in orbits),
        "all_quotient_to_embedded_invariant_sections_absent": invariant_section_count_total == 0,
        "sample_orbit_section_rows": rows,
    }


def build_payload(p2892: dict[str, Any]) -> dict[str, Any]:
    audit = quotient_section_audit()
    facts = {
        "p2892_rechecked": p2892.get("status") == "P2892_POST_PHASE_ORIGIN_INVENTORY_STATE_MAP_NO_NEW_LIVE_FRONTIER_CERTIFICATE",
        "p2888_target_count_reproduced": audit["target_triple_count"] == 600,
        "p2889_orbit_count_reproduced": audit["translation_orbit_count"] == 50,
        "all_orbits_free": audit["all_orbits_free"],
        "no_source_neutral_quotient_to_embedded_section": audit["all_quotient_to_embedded_invariant_sections_absent"],
    }
    return {
        "status": "P2893_FREE_ORBIT_SECTION_OBSTRUCTION_THEOREM_NO_CLOSURE",
        "input_hashes": {"P2888": sha(P2888), "P2889": sha(P2889), "P2892": sha(P2892)},
        "free_orbit_section_obstruction_theorem": {
            "input_status_rechecked": p2892.get("status"),
            "candidate_class": "source-neutral quotient-to-embedded representative section for free Z12 translation orbits of the P2888/P2889 9/5 carriers",
            "finite_section_audit": audit,
            "proof_certificate": {
                "group_action_fact": "For a free transitive Z12 orbit O, a quotient-level invariant section */Z12 -> O would have to choose a point fixed by every translation.  Free orbits have no such point.",
                "finite_result": "The recomputed 600 carriers form 50 orbits, every orbit has size 12 and point stabilizer size 1, and the total invariant quotient-to-embedded section count is 0.",
                "endomorphism_boundary": "There are 12 equivariant endomorphisms per free orbit after a basepoint is chosen, but these are shifted representatives and do not descend to a source-neutral quotient section.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_source_neutral_embedded_representative_law": False,
            "exports_translation_breaking_source": False,
            "exports_nonimported_9_over_5_variational_density": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "source_neutral_embedded_representative_law_exported": False,
                "translation_breaking_source_exported": False,
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
            "reason": "P2893 turns the P2889/P2892 orbit obstruction into an explicit finite section theorem.  All 50 target carrier orbits are free Z12 orbits, so no quotient-level source-neutral law has an invariant embedded representative section; the audited section count is exactly 0.  Equivariant endomorphisms exist only after a basepoint/imported representative is chosen, so they do not source the missing origin/phase.",
            "next_honest_step": "Do not try to repair the 9/5 carrier lane by choosing an equivariant orbit endomorphism, canonical representative, lexicographic representative, or other quotient-section convention.  A next proof-grade move must introduce one genuine translation-breaking strict source with a nonconventional phase/sign and a coupling theorem to the 9/5 variational density, or pivot to a genuinely different typed object outside quotient-section/support/orbit/Fourier/inventory constructions; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["free_orbit_section_obstruction_theorem"]["finite_section_audit"]
    lines = [
        "# P2893/S1843 free-orbit section obstruction theorem",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite section audit",
        f"- target triples reproduced: `{audit['target_triple_count']}`",
        f"- translation orbit count: `{audit['translation_orbit_count']}`",
        f"- orbit size histogram: `{audit['orbit_size_histogram']}`",
        f"- global point stabilizer histogram: `{audit['point_stabilizer_size_histogram_global']}`",
        f"- invariant quotient-to-embedded section count: `{audit['quotient_to_embedded_invariant_section_count_total']}`",
        f"- equivariant endomorphisms after basepoint import: `{audit['equivariant_endomorphism_count_total_not_source_neutral']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2892))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2893/S1843 free-orbit section obstruction theorem",
        "## P2893/S1843 free-orbit section obstruction theorem\n\n"
        "`P2893/S1843` upgrades the P2889/P2892 orbit obstruction to a finite section theorem.  The recomputed `600` `9/5` carriers form `50` free `Z12` translation orbits; every point stabilizer has size `1`, so a quotient-level source-neutral section to an embedded representative has `0` invariant choices.  Equivariant orbit endomorphisms exist only after importing a basepoint/representative and therefore do not export a strict phase/origin source, nonimported `9/5` variational density, localized action density, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2893/S1843 free-orbit section obstruction `L_total` guard",
        "## P2893/S1843 free-orbit section obstruction `L_total` guard\n\n"
        "`P2893/S1843` is an invariant-section no-go theorem for source-neutral representative choice, not a strict action term.  It supplies no localized unit-bearing density and no variational chain rule into nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current free-orbit section obstruction guardrail (P2893/S1843, 2026-06-19)",
        "## Current free-orbit section obstruction guardrail (P2893/S1843, 2026-06-19)\n\n"
        "- P2893 upgrades the P2889/P2892 orbit obstruction to a finite section theorem for source-neutral quotient-to-embedded representative selection on the `9/5` carrier orbits.\n"
        "- The recomputed `600` carriers form `50` free `Z12` orbits with point stabilizer size `1`; the invariant quotient-to-embedded section count is `0`.\n"
        "- Do not promote equivariant orbit endomorphisms, canonical/lexicographic representatives, quotient-section conventions, support/origin/density choices, Fourier phase/power, or inventory hits to strict phase/origin sourcehood, strict damping/compression bridge, selector closure, role transfer, `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must introduce a genuine translation-breaking strict source with nonconventional phase/sign and a coupling theorem to the `9/5` variational density, pivot to a genuinely different typed object outside quotient-section/support/orbit/Fourier/inventory constructions, or preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
