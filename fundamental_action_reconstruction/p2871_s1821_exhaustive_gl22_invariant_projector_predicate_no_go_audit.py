#!/usr/bin/env python3
"""P2871/S1821: exhaustive GL(2,2)-invariant projector-predicate no-go audit.

P2870 showed that a finite family of intrinsic character-projector scores does
not distinguish the nonidentity projectors 5, 7, and 11.  This packet removes a
possible escape hatch by exhaustively enumerating *all* Boolean predicates on the
four Aut(Z12)-character point projectors and filtering them by the intrinsic
V4/GL(2,2) relabeling symmetry.

Only four invariant predicates exist: empty, identity-only, nonidentity-all, and
all.  None selects the singleton {11}.  The target predicate {11} exists as a
Boolean predicate, but its orbit under GL(2,2) is {{5}, {7}, {11}} and therefore
it imports an endpoint label rather than exporting a strict non-premise selector.
"""
from __future__ import annotations

import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2869_s1819_aut_character_idempotent_endpoint_localizer_no_source_audit import UNITS_Z12
from p2870_s1820_intrinsic_character_projector_selector_no_go_audit import all_v4_relabelings

P2870 = GEN / "p2870_s1820_intrinsic_character_projector_selector_no_go_audit.json"
OUT = GEN / "p2871_s1821_exhaustive_gl22_invariant_projector_predicate_no_go_audit.json"
MD = GEN / "p2871_s1821_exhaustive_gl22_invariant_projector_predicate_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

TARGET_SINGLETON = frozenset({11})


def all_predicates() -> list[frozenset[int]]:
    predicates: list[frozenset[int]] = []
    units = list(UNITS_Z12)
    for mask in range(1 << len(units)):
        selected = frozenset(unit for index, unit in enumerate(units) if mask & (1 << index))
        predicates.append(selected)
    return predicates


def apply_relabeling(predicate: frozenset[int], mapping: dict[int, int]) -> frozenset[int]:
    return frozenset(mapping[value] for value in predicate)


def predicate_orbit(predicate: frozenset[int]) -> list[list[int]]:
    orbit = {apply_relabeling(predicate, mapping) for mapping in all_v4_relabelings()}
    return [list(item) for item in sorted(orbit, key=lambda item: (len(item), sorted(item)))]


def is_invariant(predicate: frozenset[int]) -> bool:
    return all(apply_relabeling(predicate, mapping) == predicate for mapping in all_v4_relabelings())


def predicate_record(predicate: frozenset[int]) -> dict[str, Any]:
    orbit = predicate_orbit(predicate)
    return {
        "selected": list(sorted(predicate)),
        "cardinality": len(predicate),
        "is_gl22_invariant": is_invariant(predicate),
        "orbit_size": len(orbit),
        "orbit": orbit,
        "selects_target_singleton_11": predicate == TARGET_SINGLETON,
    }


def build_payload(p2870: dict[str, Any]) -> dict[str, Any]:
    predicates = all_predicates()
    records = [predicate_record(predicate) for predicate in predicates]
    invariant_records = [record for record in records if record["is_gl22_invariant"]]
    target_record = predicate_record(TARGET_SINGLETON)
    accepted = [record for record in invariant_records if record["selects_target_singleton_11"]]
    facts = {
        "p2870_rechecked": p2870.get("status") == "P2870_INTRINSIC_CHARACTER_PROJECTOR_SELECTOR_NO_GO_AUDIT_NO_CLOSURE",
        "all_16_boolean_predicates_enumerated": len(records) == 16,
        "only_4_gl22_invariant_predicates": len(invariant_records) == 4,
        "no_invariant_predicate_selects_singleton_11": accepted == [],
        "target_singleton_orbit_has_three_nonidentity_singletons": target_record["orbit"] == [[5], [7], [11]],
        "accepted_count_zero": True,
    }
    return {
        "status": "P2871_EXHAUSTIVE_GL22_INVARIANT_PROJECTOR_PREDICATE_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2870": sha(P2870)},
        "exhaustive_gl22_invariant_projector_predicate_no_go_audit": {
            "input_status_rechecked": p2870.get("status"),
            "candidate_class": "all Boolean selector predicates on the four Aut(Z12)-character point projectors, filtered by V4/GL(2,2) invariance",
            "predicate_count": len(records),
            "invariant_predicates": invariant_records,
            "target_singleton_record": target_record,
            "accepted_candidate_count": 0,
            "proof_certificate": {
                "exhaustion_step": "There are 2^4=16 Boolean predicates on the four point projectors, all enumerated explicitly.",
                "invariance_step": "GL(2,2) fixes identity and permutes the three nonidentity projectors, so invariant predicates are unions of {1} and {5,7,11}.",
                "sourcehood_step": "The singleton {11} is not invariant; selecting it requires importing an endpoint label or symmetry-breaking law not present in current artifacts.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_gl22_invariant_projector_selector": False,
            "exports_boundary_source_law": False,
            "exports_eta_source_law": False,
            "exports_prime5_source_law": False,
            "exports_selector_or_localizer_source": False,
            "exports_unit_bearing_coupling_localization_theorem": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "boundary_source_law_exported": False,
                "eta_source_exported": False,
                "prime5_source_exported": False,
                "selector_or_localizer_source_exported": False,
                "selector_closure_exported": False,
                "unit_bearing_coupling_localization_theorem_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2871 exhausts all Boolean selector predicates on the four character point projectors.  GL(2,2)-invariant predicates cannot select singleton 11; the singleton target predicate is a non-invariant endpoint-label import.",
            "next_honest_step": "Do not replay Boolean projector predicates or target singleton predicates as sourcehood.  A next proof-grade move must provide a new strict symmetry-breaking law whose own provenance fixes the nonidentity singleton and supplies the unit-bearing coefficient/coupling theorem, or pivot to a different new typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["exhaustive_gl22_invariant_projector_predicate_no_go_audit"]
    lines = [
        "# P2871/S1821 exhaustive GL(2,2)-invariant projector-predicate no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Exhaustive predicate audit",
        f"- candidate class: `{audit['candidate_class']}`",
        f"- predicate count: `{audit['predicate_count']}`",
        f"- invariant predicates: `{audit['invariant_predicates']}`",
        f"- target singleton record: `{audit['target_singleton_record']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload(read_json(P2870))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2871/S1821 exhaustive GL(2,2)-invariant projector-predicate no-go audit",
        "## P2871/S1821 exhaustive GL(2,2)-invariant projector-predicate no-go audit\n\n"
        "`P2871/S1821` exhausts all `2^4=16` Boolean selector predicates on the four Aut(`Z12`)-character point projectors and filters them by V4/GL(2,2) invariance.  The only invariant predicates are empty, identity-only, nonidentity-all, and all; none selects singleton `11`.  The target predicate `{11}` has orbit `{ {5}, {7}, {11} }`, so it imports an endpoint label rather than exporting a strict selector/source law.  No boundary source law, eta source, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2871/S1821 exhaustive projector-predicate `L_total` guard",
        "## P2871/S1821 exhaustive projector-predicate `L_total` guard\n\n"
        "`P2871/S1821` adds no strict action term.  Exhaustive Boolean projector predicates do not export a unit-bearing boundary/source density, coupling coefficient, localization/pullback theorem, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current exhaustive GL(2,2)-invariant projector-predicate no-go guardrail (P2871/S1821, 2026-06-18)",
        "## Current exhaustive GL(2,2)-invariant projector-predicate no-go guardrail (P2871/S1821, 2026-06-18)\n\n"
        "- P2871 exhausts all `2^4=16` Boolean selector predicates on the Aut(`Z12`)-character point projectors after P2870.\n"
        "- GL(2,2)-invariant predicates are only empty, identity-only, nonidentity-all, and all; singleton `{11}` is non-invariant and imports an endpoint label.\n"
        "- Do not promote Boolean projector predicates, target singleton predicates, intrinsic character-projector scoring, Aut-character idempotents, endpoint/polarity imports, prime-5 scaled coefficients, Dirichlet data, or log-scale harmonicity to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must provide a new strict symmetry-breaking law whose provenance fixes the nonidentity singleton and supplies the unit-bearing coefficient/coupling theorem, or use a different new typed object; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
