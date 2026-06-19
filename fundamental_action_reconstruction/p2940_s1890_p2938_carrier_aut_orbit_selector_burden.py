#!/usr/bin/env python3
"""P2940/S1890: P2938 carrier Aut-orbit selector burden.

P2939 proved that the P2938 carrier is algebraically ready and Aut-breaking.
P2940 makes the provenance obstruction sharper: it computes the U(12) orbit
quotient on nodes 1..11 and checks whether the carrier is orbit-constant.

If a carrier is not orbit-constant, strict provenance cannot be Aut-invariant;
it must supply an explicit selector/symmetry-breaking source for the affected
orbits.  This does not replay selector closure.  It quantifies the exact burden
that any future provenance theorem for the P2938 carrier must discharge.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2939 = GEN / "p2939_s1889_p2938_carrier_aut_breaking_provenance_boundary.json"
OUT = GEN / "p2940_s1890_p2938_carrier_aut_orbit_selector_burden.json"
MD = GEN / "p2940_s1890_p2938_carrier_aut_orbit_selector_burden.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MODULUS = 12
NODES = list(range(1, 12))
UNITS = [1, 5, 7, 11]
PRIME_VECTOR = [1, 2, 2, 2, 2]
PRIMES = [2, 3, 5, 7, 11]


def factor_vector(n: int) -> list[int]:
    remaining = n
    vector = []
    for prime in PRIMES:
        exponent = 0
        while remaining % prime == 0:
            remaining //= prime
            exponent += 1
        vector.append(exponent)
    if remaining != 1:
        raise ValueError(f"node {n} has unaudited factor {remaining}")
    return vector


def carrier_value(n: int) -> int:
    return sum(exp * coord for exp, coord in zip(factor_vector(n), PRIME_VECTOR))


def orbit_of(seed: int) -> list[int]:
    return sorted({(u * seed) % MODULUS for u in UNITS})


def orbit_rows() -> list[dict[str, Any]]:
    seen: set[int] = set()
    rows = []
    for node in NODES:
        if node in seen:
            continue
        orbit = orbit_of(node)
        seen.update(orbit)
        values = {member: carrier_value(member) for member in orbit}
        distinct_values = sorted(set(values.values()))
        rows.append({
            "orbit_representative": node,
            "orbit": orbit,
            "carrier_values": values,
            "distinct_values": distinct_values,
            "orbit_constant": len(distinct_values) == 1,
            "selector_burden": len(distinct_values) > 1,
        })
    return rows


def acceptance_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    burden_rows = [row for row in rows if row["selector_burden"]]
    return [
        {
            "criterion": "Aut_orbit_quotient_computed",
            "satisfied": True,
            "evidence": "U(12) orbits on nodes 1..11 are explicitly enumerated",
        },
        {
            "criterion": "carrier_orbit_constant",
            "satisfied": not burden_rows,
            "evidence": f"{len(burden_rows)} orbits have nonconstant carrier values",
        },
        {
            "criterion": "selector_burden_quantified",
            "satisfied": True,
            "evidence": "all nonconstant orbits are listed with their value splits",
        },
        {
            "criterion": "strict_selector_or_symmetry_breaking_source_exported",
            "satisfied": False,
            "evidence": "no current artifact exports a strict source selecting representatives inside the nonconstant orbits",
        },
        {
            "criterion": "strict_nadsoliton_provenance_theorem_for_P2938_carrier",
            "satisfied": False,
            "evidence": "a provenance theorem would have to include the missing selector/symmetry-breaking source",
        },
    ]


def build_payload(p2939: dict[str, Any]) -> dict[str, Any]:
    rows = orbit_rows()
    burden_rows = [row for row in rows if row["selector_burden"]]
    criteria = acceptance_rows(rows)
    accepted = all(row["satisfied"] for row in criteria)
    return {
        "status": "P2940_P2938_CARRIER_AUT_ORBIT_SELECTOR_BURDEN_NO_STRICT_PROVENANCE",
        "input_hashes": {"P2939": hashlib.sha256(P2939.read_bytes()).hexdigest() if P2939.exists() else None},
        "constructed_theoretical_objects": {
            "carrier_under_test": "P2938_UnitCharacter_Enriched_Z12_PrimeCoordinate_Carrier",
            "prime_coordinate_vector_order_2_3_5_7_11": PRIME_VECTOR,
            "aut_orbit_rows": rows,
            "selector_burden_rows": burden_rows,
            "acceptance_rows": criteria,
        },
        "orbit_certificate": {
            "orbit_count": len(rows),
            "selector_burden_orbit_count": len(burden_rows),
            "orbit_constant_count": sum(1 for row in rows if row["orbit_constant"]),
            "all_orbits_constant": not burden_rows,
            "strict_selector_source_exported": False,
            "strict_provenance_theorem_exported": accepted,
        },
        "decision": {
            "positive_witnesses": {
                "Aut_orbit_quotient_computed": True,
                "selector_burden_quantified": True,
                "P2938_carrier_requires_symmetry_breaking_for_provenance": bool(burden_rows),
            },
            "negative_export_flags": {
                "strict_selector_source_exported": False,
                "strict_aut_breaking_prime_coordinate_source_law_exported": False,
                "strict_prime_log_value_source_exported": False,
                "strict_delta_eta_source_exported": False,
                "strict_beta_eta_coupling_theorem_exported": False,
                "strict_damping_beta_eta_source_packet_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_hamiltonian_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2940 computes the U(12) orbit quotient for the exact P2938 carrier and finds nonconstant carrier values on multiple orbits.  Therefore any strict provenance theorem for this carrier must include a genuine selector/symmetry-breaking source.  No such source is exported here, so P2938 remains a finite carrier only.",
            "next_honest_step": "A next admissible move must either export a strict selector/symmetry-breaking provenance theorem for the listed nonconstant orbits, or pivot to a genuinely new typed object outside the prime-coordinate carrier/provenance lane.  Otherwise preserve the P2929-P2940 no-new-live-frontier boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["orbit_certificate"]
    lines = [
        "# P2940/S1890 P2938 carrier Aut-orbit selector burden",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Orbit certificate",
        f"- U(12) orbit count: `{cert['orbit_count']}`",
        f"- selector-burden orbits: `{cert['selector_burden_orbit_count']}`",
        f"- orbit-constant rows: `{cert['orbit_constant_count']}`",
        f"- all orbits constant: `{cert['all_orbits_constant']}`",
        f"- strict selector source exported: `{cert['strict_selector_source_exported']}`",
        f"- strict provenance theorem exported: `{cert['strict_provenance_theorem_exported']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2939))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2940/S1890 P2938 carrier Aut-orbit selector burden", "## P2940/S1890 P2938 carrier Aut-orbit selector burden\n\n`P2940/S1890` computes the `U(12)` orbit quotient for the exact P2938 carrier `[1,2,2,2,2]`.  The carrier is not orbit-constant on several orbits, so any strict provenance theorem for this carrier must include a genuine selector/symmetry-breaking source for the listed nonconstant orbits.  No such source is exported; therefore P2938 remains a finite carrier only and no strict `L_p`, damping packet, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2940/S1890 Aut-orbit selector burden `L_total` guard", "## P2940/S1890 Aut-orbit selector burden `L_total` guard\n\n`P2940/S1890` shows that the P2938 carrier requires selector/symmetry-breaking provenance on nonconstant `U(12)` orbits.  Since no strict selector/provenance theorem or damping coupling is exported, the carrier cannot enter nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current P2938 Aut-orbit selector-burden guardrail (P2940/S1890, 2026-06-19)", "## Current P2938 Aut-orbit selector-burden guardrail (P2940/S1890, 2026-06-19)\n\n- P2940 computes the `U(12)` orbit quotient of the exact P2938 carrier `[1,2,2,2,2]` on nodes `1..11`.\n- The carrier is not orbit-constant on multiple orbits; strict provenance would therefore require a genuine selector/symmetry-breaking source for those listed nonconstant orbits.\n- P2940 does not export selector closure, strict provenance, delta/eta source, beta/eta coupling, strict `L_p`, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- A next admissible move must export that selector/provenance theorem or pivot to a genuinely new typed object outside the prime-coordinate carrier/provenance lane; otherwise preserve the P2929-P2940 no-new-live-frontier boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
