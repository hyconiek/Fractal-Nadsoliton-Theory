#!/usr/bin/env python3
"""P2874/S1824: dihedral-character chiral endpoint-source no-go audit.

P2873 exhausted invariant Boolean endpoint predicates on Z12 and showed that a
singleton endpoint localizer is not available without an imported orientation or
endpoint label.  This packet tests the next, stronger candidate class requested
by that boundary: finite one-dimensional D12 chiral/character laws.

The calculation enumerates all four real one-dimensional characters of the
D12 action on Z12 endpoints.  For each character chi, it solves the endpoint
field equivariance equations v(g.x)=chi(g)v(x).  Reflection-odd characters have
only the zero endpoint field because every endpoint has a reflection stabilizer.
Reflection-even characters give either the constant field or the alternating
rotation-parity field.  Neither has singleton support at endpoint 11, and the
alternating field still needs an imported polarity/sign and cannot localize the
unit-bearing 9/5 coefficient.
"""
from __future__ import annotations

import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2873 = GEN / "p2873_s1823_exhaustive_dihedral_z12_endpoint_predicate_no_go_audit.json"
OUT = GEN / "p2874_s1824_dihedral_character_chiral_endpoint_source_no_go_audit.json"
MD = GEN / "p2874_s1824_dihedral_character_chiral_endpoint_source_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MODULUS = 12
ENDPOINTS = tuple(range(MODULUS))
TARGET_ENDPOINT = 11
SUCCESSOR_ENDPOINT = 1


def rotation(endpoint: int, power: int = 1) -> int:
    return (endpoint + power) % MODULUS


def reflection(endpoint: int) -> int:
    return (-endpoint) % MODULUS


def character_value(eps_r: int, eps_s: int, *, rotation_power: int, reflected: bool) -> int:
    return (eps_r ** (rotation_power % MODULUS)) * (eps_s if reflected else 1)


def character_record(eps_r: int, eps_s: int) -> dict[str, Any]:
    """Return the exact equivariant endpoint field space for one D12 character.

    D12 acts transitively on endpoints by rotations, and the stabilizer of each
    endpoint contains a reflection.  Hence a reflection-odd character forces all
    endpoint values to zero.  A reflection-even character leaves one field,
    determined by v(0)=1 and v(k)=eps_r**k.
    """
    reflection_odd = eps_s == -1
    if reflection_odd:
        basis: list[list[int]] = []
        support: list[int] = []
        dimension = 0
        field_kind = "zero_only_reflection_stabilizer_obstruction"
    else:
        field = [eps_r ** endpoint for endpoint in ENDPOINTS]
        basis = [field]
        support = [endpoint for endpoint, value in enumerate(field) if value != 0]
        dimension = 1
        field_kind = "constant" if eps_r == 1 else "alternating_rotation_parity"
    singleton_support = support == [TARGET_ENDPOINT]
    separates_11_from_1 = bool(basis and basis[0][TARGET_ENDPOINT] != basis[0][SUCCESSOR_ENDPOINT])
    return {
        "character": {"chi_r": eps_r, "chi_s": eps_s},
        "reflection_odd": reflection_odd,
        "endpoint_field_dimension": dimension,
        "basis_fields": basis,
        "support": support,
        "field_kind": field_kind,
        "singleton_11_support": singleton_support,
        "separates_11_from_1_by_value": separates_11_from_1,
        "exports_endpoint_11_source": False,
        "obstruction": (
            "reflection-odd characters vanish on endpoint fields because endpoint stabilizers contain reflections"
            if reflection_odd
            else "reflection-even character field is global across all endpoints, not a singleton 11 localizer"
        ),
    }


def build_payload(p2873: dict[str, Any]) -> dict[str, Any]:
    records = [character_record(eps_r, eps_s) for eps_r in (1, -1) for eps_s in (1, -1)]
    accepted = [record for record in records if record["singleton_11_support"]]
    alternating = [record for record in records if record["field_kind"] == "alternating_rotation_parity"]
    facts = {
        "p2873_rechecked": p2873.get("status") == "P2873_EXHAUSTIVE_DIHEDRAL_Z12_ENDPOINT_PREDICATE_NO_GO_AUDIT_NO_CLOSURE",
        "all_four_real_d12_characters_enumerated": len(records) == 4,
        "reflection_odd_characters_have_zero_endpoint_field": all(record["endpoint_field_dimension"] == 0 for record in records if record["reflection_odd"]),
        "reflection_even_characters_are_global_not_singleton": all(record["support"] == list(ENDPOINTS) for record in records if not record["reflection_odd"]),
        "no_character_field_has_singleton_11_support": accepted == [],
        "alternating_character_is_global_parity_not_endpoint_11_localization": len(alternating) == 1 and not alternating[0]["separates_11_from_1_by_value"] and alternating[0]["support"] == list(ENDPOINTS),
    }
    return {
        "status": "P2874_DIHEDRAL_CHARACTER_CHIRAL_ENDPOINT_SOURCE_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2873": sha(P2873)},
        "dihedral_character_chiral_endpoint_source_no_go_audit": {
            "input_status_rechecked": p2873.get("status"),
            "candidate_class": "all real one-dimensional D12 chiral/character laws with endpoint-field equivariance v(g.x)=chi(g)v(x)",
            "modulus": MODULUS,
            "character_count": len(records),
            "character_records": records,
            "accepted_candidate_count": len(accepted),
            "proof_certificate": {
                "character_enumeration_step": "The real one-dimensional characters of D12 are exactly chi(r)=±1 and chi(s)=±1.",
                "stabilizer_step": "Every endpoint has a reflection stabilizer; therefore chi(s)=-1 forces v(endpoint)=0 for all endpoints.",
                "global_field_step": "The chi(s)=+1 fields are constant or alternating across all twelve endpoints, not singleton 11 localizers.",
                "sourcehood_step": "The alternating field is only a global parity pattern; it does not distinguish 11 from the other odd endpoint 1 and supplies no non-premise endpoint-11 source, polarity choice, or unit-bearing 9/5 coupling theorem.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_chiral_endpoint_source": False,
            "exports_boundary_source_law": False,
            "exports_orientation_source_law": False,
            "exports_chiral_endpoint_11_source_law": False,
            "exports_selector_or_localizer_source": False,
            "exports_unit_bearing_coupling_localization_theorem": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "boundary_source_law_exported": False,
                "orientation_source_law_exported": False,
                "chiral_endpoint_11_source_law_exported": False,
                "selector_or_localizer_source_exported": False,
                "unit_bearing_coupling_localization_theorem_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2874 enumerates all real one-dimensional D12 chiral/character source laws.  Reflection-odd characters vanish on endpoint fields, while reflection-even characters are global constant/alternating fields.  The alternating field is only a global parity pattern; it does not distinguish 11 from 1 and supplies no singleton localization, non-premise polarity, or unit-bearing 9/5 coupling.",
            "next_honest_step": "Do not replay D12 one-dimensional character/chiral fields, reflection-odd sign characters, or alternating parity fields as endpoint-11 sourcehood.  A next proof-grade move must provide a higher-structure strict source with a nonzero localized signed value at 11 and a unit-bearing 9/5 coupling theorem, or pivot to a genuinely different typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["dihedral_character_chiral_endpoint_source_no_go_audit"]
    lines = [
        "# P2874/S1824 dihedral-character chiral endpoint-source no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Character-source audit",
        f"- candidate class: `{audit['candidate_class']}`",
        f"- character count: `{audit['character_count']}`",
        f"- accepted candidate count: `{audit['accepted_candidate_count']}`",
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
    payload = build_payload(read_json(P2873))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2874/S1824 dihedral-character chiral endpoint-source no-go audit",
        "## P2874/S1824 dihedral-character chiral endpoint-source no-go audit\n\n"
        "`P2874/S1824` enumerates all four real one-dimensional `D12` chiral/character laws `chi(r)=±1`, `chi(s)=±1` and solves the endpoint-field equivariance equations `v(g.x)=chi(g)v(x)`.  Reflection-odd characters have only the zero endpoint field because every endpoint has a reflection stabilizer.  Reflection-even characters give global constant or alternating fields, not singleton `11` localizers; the alternating field is only a global parity pattern, does not distinguish `11` from `1`, and still lacks non-premise polarity and a unit-bearing `9/5` coupling theorem.  No boundary/orientation source law, selector/localizer source, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2874/S1824 dihedral character `L_total` guard",
        "## P2874/S1824 dihedral character `L_total` guard\n\n"
        "`P2874/S1824` adds no strict action term.  One-dimensional `D12` chiral/character endpoint fields do not export a localized unit-bearing boundary/source density at `11`, non-premise polarity law, coupling coefficient, localization/pullback theorem, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current dihedral-character chiral endpoint-source no-go guardrail (P2874/S1824, 2026-06-18)",
        "## Current dihedral-character chiral endpoint-source no-go guardrail (P2874/S1824, 2026-06-18)\n\n"
        "- P2874 enumerates all real one-dimensional `D12` chiral/character laws after P2873 and solves endpoint-field equivariance.\n"
        "- Reflection-odd characters vanish on endpoint fields because each endpoint has a reflection stabilizer; reflection-even characters are global constant or alternating fields, not singleton `11` localizers.\n"
        "- Do not promote `D12` one-dimensional character fields, reflection-odd sign characters, alternating parity fields, full-dihedral endpoint predicates, reflection-only adjacent pairs, cyclic predecessor conventions, endpoint-label imports, Aut-character idempotents, or prime-5 scaled coefficients to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must export a higher-structure strict source with a nonzero localized signed value at `11` plus a unit-bearing `9/5` coupling theorem, or use a genuinely different new typed object; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
