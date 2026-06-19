#!/usr/bin/env python3
"""P2889/S1839: translation-orbit source-law 9/5 no-go audit.

P2888 showed that non-C12 supports plus intrinsic distance-profile origins and
bounded integer densities can represent 9/5, but only as many unsourced carriers.
This packet tests the next stricter obligation: can the P2888 target carriers be
made into a source-neutral strict law without importing an absolute origin?

The audit regenerates exactly the P2888 target triples with support size 5 and
bounded density values {-2,-1,0,1,2} whose uniform unit average is 9/5.  It then
quotients them by the Z12 translation action on support, selected origin, and
localized density.  If a carrier is in a nontrivial translation orbit, any
translation-neutral rule can select at most the orbit; selecting one embedded
representative is an origin import.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2888_s1838_non_c12_origin_law_unit_coupling_9_over_5_no_go_audit import (
    DENSITY_VALUES,
    N,
    P2887,
    TARGET,
    intrinsic_unique_origins,
    sample_assignment_with_sum,
    support_points,
)

P2888 = GEN / "p2888_s1838_non_c12_origin_law_unit_coupling_9_over_5_no_go_audit.json"
OUT = GEN / "p2889_s1839_translation_orbit_source_law_9_over_5_no_go_audit.json"
MD = GEN / "p2889_s1839_translation_orbit_source_law_9_over_5_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


@dataclass(frozen=True, order=True)
class TargetTriple:
    support: tuple[int, ...]
    origin: int
    density_by_node: tuple[tuple[int, int], ...]

    def rotate(self, shift: int) -> "TargetTriple":
        density = tuple(sorted((((node + shift) % N), value) for node, value in self.density_by_node))
        support = tuple(node for node, _value in density)
        return TargetTriple(support=support, origin=(self.origin + shift) % N, density_by_node=density)

    def primitive_record(self) -> dict[str, Any]:
        return {
            "support": list(self.support),
            "origin": self.origin,
            "density_by_node": [[node, value] for node, value in self.density_by_node],
            "unit_average": str(TARGET),
        }


def density_assignments_sum(k: int, target_sum: int) -> list[tuple[int, ...]]:
    assignments: list[tuple[int, ...]] = []

    def rec(prefix: tuple[int, ...], remaining_slots: int, remaining_sum: int) -> None:
        if remaining_slots == 0:
            if remaining_sum == 0:
                assignments.append(prefix)
            return
        rest = remaining_slots - 1
        low = min(DENSITY_VALUES) * rest
        high = max(DENSITY_VALUES) * rest
        for value in DENSITY_VALUES:
            new_remaining = remaining_sum - value
            if low <= new_remaining <= high:
                rec(prefix + (value,), rest, new_remaining)

    rec((), k, target_sum)
    return assignments


def regenerate_p2888_target_triples() -> list[TargetTriple]:
    triples: list[TargetTriple] = []
    k = TARGET.denominator
    target_sum = TARGET.numerator
    assignments = density_assignments_sum(k, target_sum)
    for mask in range(1, 1 << N):
        points = support_points(mask)
        if len(points) != k:
            continue
        origins = intrinsic_unique_origins(points)
        if len(origins) != 1:
            continue
        for assignment in assignments:
            triples.append(
                TargetTriple(
                    support=points,
                    origin=origins[0],
                    density_by_node=tuple(zip(points, assignment)),
                )
            )
    return sorted(triples)


def orbit(triple: TargetTriple) -> tuple[TargetTriple, ...]:
    return tuple(sorted({triple.rotate(shift) for shift in range(N)}))


def orbit_audit(triples: list[TargetTriple]) -> dict[str, Any]:
    unseen = set(triples)
    orbits: list[tuple[TargetTriple, ...]] = []
    while unseen:
        current = min(unseen)
        orb = orbit(current)
        orbits.append(orb)
        unseen.difference_update(orb)
    sizes = [len(orb) for orb in orbits]
    stabilizer_sizes = [N // size for size in sizes]
    orbit_samples = []
    for orb in orbits[:8]:
        representative = min(orb)
        orbit_samples.append(
            {
                "orbit_size": len(orb),
                "translation_stabilizer_size": N // len(orb),
                "representative": representative.primitive_record(),
            }
        )
    return {
        "target_triple_count": len(triples),
        "translation_orbit_count": len(orbits),
        "translation_orbit_size_histogram": {str(size): sizes.count(size) for size in sorted(set(sizes))},
        "translation_stabilizer_size_histogram": {str(size): stabilizer_sizes.count(size) for size in sorted(set(stabilizer_sizes))},
        "all_target_triples_have_free_translation_orbit": all(size == N for size in sizes),
        "sample_orbits": orbit_samples,
    }


def build_payload(p2888: dict[str, Any]) -> dict[str, Any]:
    triples = regenerate_p2888_target_triples()
    audit = orbit_audit(triples)
    previous_target_count = p2888["non_c12_origin_law_unit_coupling_9_over_5_no_go_audit"]["unit_coupling_audit"]["target_9_over_5_record_count"]
    facts = {
        "p2888_rechecked": p2888.get("status") == "P2888_NON_C12_ORIGIN_LAW_UNIT_COUPLING_9_OVER_5_NO_GO_AUDIT_NO_CLOSURE",
        "p2888_target_count_reproduced": audit["target_triple_count"] == previous_target_count == 600,
        "all_target_triples_are_in_free_translation_orbits": audit["all_target_triples_have_free_translation_orbit"],
        "translation_neutral_rule_can_select_only_orbit_not_representative": audit["all_target_triples_have_free_translation_orbit"],
        "strict_source_law_missing": True,
    }
    return {
        "status": "P2889_TRANSLATION_ORBIT_SOURCE_LAW_9_OVER_5_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2887": sha(P2887), "P2888": sha(P2888)},
        "translation_orbit_source_law_9_over_5_no_go_audit": {
            "input_status_rechecked": p2888.get("status"),
            "candidate_class": "P2888 target support-origin-density triples quotiented by Z12 translation action",
            "orbit_audit": audit,
            "proof_certificate": {
                "reproduction_result": "The P2888 target carrier family is exactly reproduced: 600 embedded target triples.",
                "orbit_result": "The 600 triples form translation orbits of size 12; no target triple has a nontrivial translation stabilizer.",
                "sourcehood_obstruction": "A translation-neutral source law can select at most a 12-element orbit.  Selecting one embedded representative requires an absolute origin/phase convention or an additional strict source not present in P2888.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_translation_neutral_strict_source_law": False,
            "exports_unique_embedded_support_origin_density_triple": False,
            "exports_nonimported_9_over_5_variational_chain_rule": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "translation_neutral_strict_source_law_exported": False,
                "unique_embedded_support_origin_density_triple_exported": False,
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
            "reason": "P2889 reproduces the P2888 9/5 carriers and quotients them by Z12 translations.  All carriers lie in free 12-element translation orbits, so a translation-neutral law can pick at most an orbit, not a unique embedded support-origin-density representative.  A representative choice would import an absolute origin or a new strict symmetry-breaking source.",
            "next_honest_step": "Do not replay P2888 carrier enumeration, translation orbit representatives, distance-profile origin selection, bounded density coefficients, or C12-neutral unit measures as strict sourcehood.  A next proof-grade move must supply one new strict translation-breaking source/phase law coupled to the 9/5 carrier and variational density, or pivot to a genuinely different typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["translation_orbit_source_law_9_over_5_no_go_audit"]["orbit_audit"]
    lines = [
        "# P2889/S1839 translation-orbit source-law 9/5 no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite orbit audit",
        f"- target triples reproduced: `{audit['target_triple_count']}`",
        f"- translation orbit count: `{audit['translation_orbit_count']}`",
        f"- translation orbit size histogram: `{audit['translation_orbit_size_histogram']}`",
        f"- translation stabilizer size histogram: `{audit['translation_stabilizer_size_histogram']}`",
        f"- all target triples free under translation: `{audit['all_target_triples_have_free_translation_orbit']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2888))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2889/S1839 translation-orbit source-law 9/5 no-go audit",
        "## P2889/S1839 translation-orbit source-law 9/5 no-go audit\n\n"
        "`P2889/S1839` reproduces the `P2888` target carrier family and quotients the support/origin/density triples by the `Z12` translation action.  The `600` embedded `9/5` carriers form `50` free translation orbits of size `12`; hence a translation-neutral law can select at most an orbit, not a unique embedded representative, without importing an absolute origin or a new strict symmetry-breaking source.  No nonimported `9/5` variational chain rule, localized action density, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2889/S1839 translation-orbit source-law `L_total` guard",
        "## P2889/S1839 translation-orbit source-law `L_total` guard\n\n"
        "`P2889/S1839` adds no strict action term.  The `P2888` carriers survive only as free translation orbits; choosing one embedded representative requires an extra translation-breaking source and does not provide a unit-bearing localized action density or variational chain rule into nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current translation-orbit source-law 9/5 no-go guardrail (P2889/S1839, 2026-06-19)",
        "## Current translation-orbit source-law 9/5 no-go guardrail (P2889/S1839, 2026-06-19)\n\n"
        "- P2889 reproduces the P2888 `9/5` support-origin-density carriers and quotients them by the `Z12` translation action.\n"
        "- The `600` embedded carriers form `50` free translation orbits of size `12`; a translation-neutral law can select an orbit but not a unique embedded representative without importing an absolute origin or a new strict translation-breaking source.\n"
        "- Do not promote P2888 carrier enumeration, translation orbit representatives, distance-profile origin selection, bounded density coefficients, or `C12`-neutral unit measures to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must supply one new strict translation-breaking source/phase law coupled to the `9/5` carrier and variational density, pivot to a genuinely different typed object, or preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
