#!/usr/bin/env python3
"""P2882/S1832: Euler source-ratio 9/5 forcing no-go audit.

P2881 showed that a finite family of source-neutral local quadratic/unit laws does
not derive center coefficient 9/5.  This packet tests the next more analytic
possibility: perhaps an Euler equation from a one-coordinate quadratic action
can force 9/5 through a stiffness/source balance.

For a scalar quadratic action E(x)=1/2 A x^2 - J x with positive integer
stiffness A and integer source J, the exact Euler equation is A*x=J.  Hence
x=9/5 occurs exactly when J/A=9/5, i.e. after reduction the exported
source-to-stiffness ratio is already 9:5.  The finite scan below confirms this
on a bounded integer box and records the proof obligation: a future closure must
export a strict law for that source/stiffness ratio rather than merely choosing
it.
"""
from __future__ import annotations

import json
from fractions import Fraction
from math import gcd
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2881 = GEN / "p2881_s1831_variational_unit_law_9_over_5_derivation_no_go_audit.json"
OUT = GEN / "p2882_s1832_euler_source_ratio_9_over_5_forcing_no_go_audit.json"
MD = GEN / "p2882_s1832_euler_source_ratio_9_over_5_forcing_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

TARGET = Fraction(9, 5)
STIFFNESS_RANGE = range(1, 61)
SOURCE_RANGE = range(-108, 109)


def primitive_ratio(source: int, stiffness: int) -> tuple[int, int]:
    divisor = gcd(abs(source), stiffness)
    if divisor == 0:
        return (0, 1)
    return (source // divisor, stiffness // divisor)


def euler_records() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for stiffness in STIFFNESS_RANGE:
        for source in SOURCE_RANGE:
            solution = Fraction(source, stiffness)
            primitive = primitive_ratio(source, stiffness)
            records.append(
                {
                    "stiffness_A": stiffness,
                    "source_J": source,
                    "solution": str(solution),
                    "primitive_source_to_stiffness_ratio": list(primitive),
                    "forces_9_over_5": solution == TARGET,
                    "imports_9_to_5_ratio": primitive == (9, 5),
                }
            )
    return records


def build_payload(p2881: dict[str, Any]) -> dict[str, Any]:
    records = euler_records()
    target_records = [record for record in records if record["forces_9_over_5"]]
    primitive_target_records = [record for record in target_records if record["imports_9_to_5_ratio"]]
    without_import_records = [record for record in records if not record["imports_9_to_5_ratio"] and record["forces_9_over_5"]]
    target_stiffnesses = [record["stiffness_A"] for record in target_records]
    target_sources = [record["source_J"] for record in target_records]
    facts = {
        "p2881_rechecked": p2881.get("status") == "P2881_VARIATIONAL_UNIT_LAW_9_OVER_5_DERIVATION_NO_GO_AUDIT_NO_CLOSURE",
        "integer_euler_box_checked": len(records) == len(STIFFNESS_RANGE) * len(SOURCE_RANGE),
        "target_occurs_only_when_source_stiffness_ratio_imports_9_to_5": len(target_records) == len(primitive_target_records),
        "no_9_over_5_without_imported_9_to_5_source_ratio": len(without_import_records) == 0,
        "target_family_is_scaled_not_canonical": len(set((r["source_J"], r["stiffness_A"]) for r in target_records)) > 1,
    }
    return {
        "status": "P2882_EULER_SOURCE_RATIO_9_OVER_5_FORCING_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2881": sha(P2881)},
        "euler_source_ratio_9_over_5_forcing_no_go_audit": {
            "input_status_rechecked": p2881.get("status"),
            "candidate_class": "scalar quadratic Euler laws A*x=J with positive integer stiffness A and integer source J",
            "stiffness_range": [min(STIFFNESS_RANGE), max(STIFFNESS_RANGE)],
            "source_range": [min(SOURCE_RANGE), max(SOURCE_RANGE)],
            "candidate_count": len(records),
            "target_solution": "9/5",
            "target_record_count": len(target_records),
            "target_records_with_primitive_9_to_5_ratio_count": len(primitive_target_records),
            "target_records_without_primitive_9_to_5_ratio_count": len(without_import_records),
            "target_stiffnesses": target_stiffnesses,
            "target_sources": target_sources,
            "sample_target_records": target_records[:12],
            "proof_certificate": {
                "euler_identity": "For E(x)=1/2 A x^2 - J x, the Euler equation is A*x=J, so x=J/A.",
                "exact_target_condition": "x=9/5 if and only if 5J=9A; for integer A,J this reduces to primitive source/stiffness ratio 9:5.",
                "finite_result": "Every target record in the scanned box has primitive ratio 9:5; no record forces 9/5 without importing that ratio, and the target records form a scaled family rather than a canonical source.",
                "sourcehood_step": "A strict derivation must export the source/stiffness ratio 9:5 from strict dimensional/unit data; the Euler equation alone only transmits an already chosen ratio.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_independent_euler_source_law_for_9_over_5": False,
            "exports_source_stiffness_ratio_9_to_5": False,
            "exports_unit_bearing_9_over_5_coupling_theorem": False,
            "exports_unit_bearing_action_density": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "independent_euler_source_law_exported": False,
                "source_stiffness_ratio_9_to_5_exported": False,
                "unit_bearing_9_over_5_coupling_theorem_exported": False,
                "unit_bearing_action_density_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2882 audits the analytic Euler source-ratio route after P2881.  The exact equation A*x=J forces x=9/5 only when the source/stiffness ratio J:A is already 9:5.  The finite integer scan confirms that every target record imports this primitive ratio and that the target records are a scaled family, not a canonical strict source.",
            "next_honest_step": "Do not replay scalar Euler laws A*x=J, local quadratic minimization, denominator-5 coefficient boxes, endpoint pins, or imported 9:5 source/stiffness ratios as sourcehood.  A next proof-grade move must export a strict dimensional/unit law that fixes the primitive source/stiffness ratio 9:5 from independent data, or pivot to a genuinely different typed object outside the endpoint/coefficient/source-ratio family; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["euler_source_ratio_9_over_5_forcing_no_go_audit"]
    lines = [
        "# P2882/S1832 Euler source-ratio 9/5 forcing no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Euler source-ratio audit",
        f"- candidate class: `{audit['candidate_class']}`",
        f"- candidate count: `{audit['candidate_count']}`",
        f"- target record count: `{audit['target_record_count']}`",
        f"- target records with primitive 9:5 ratio: `{audit['target_records_with_primitive_9_to_5_ratio_count']}`",
        f"- target records without primitive 9:5 ratio: `{audit['target_records_without_primitive_9_to_5_ratio_count']}`",
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
    payload = build_payload(read_json(P2881))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2882/S1832 Euler source-ratio 9/5 forcing no-go audit",
        "## P2882/S1832 Euler source-ratio 9/5 forcing no-go audit\n\n"
        "`P2882/S1832` audits the analytic scalar Euler route after `P2881`: for a quadratic action `E(x)=1/2 A x^2 - J x`, the Euler equation is `A*x=J`.  Exact arithmetic and a finite integer scan show `x=9/5` occurs only when the primitive source/stiffness ratio is already `J:A=9:5`; the target records form a scaled family rather than a canonical strict source.  No independent source/stiffness law, unit-bearing `9/5` coupling theorem, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2882/S1832 Euler source-ratio 9/5 forcing `L_total` guard",
        "## P2882/S1832 Euler source-ratio 9/5 forcing `L_total` guard\n\n"
        "`P2882/S1832` adds no strict action term.  The scalar Euler identity `A*x=J` can transmit `9/5` only from an already exported source/stiffness ratio `9:5`; it does not export that ratio, a localized unit-bearing boundary/source density, a derived coupling coefficient, a nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current Euler source-ratio 9/5 forcing no-go guardrail (P2882/S1832, 2026-06-18)",
        "## Current Euler source-ratio 9/5 forcing no-go guardrail (P2882/S1832, 2026-06-18)\n\n"
        "- P2882 tests the post-P2881 analytic Euler route `A*x=J` for a scalar quadratic action and scans positive integer stiffness/source pairs exactly.\n"
        "- `x=9/5` occurs only when the primitive source/stiffness ratio is already `J:A=9:5`; the target records are scaled representatives, not a canonical strict source.\n"
        "- Do not promote scalar Euler laws, local quadratic minimization, denominator-5 coefficient boxes, endpoint pins, or imported `9:5` source/stiffness ratios to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must export a strict dimensional/unit law fixing primitive ratio `9:5` from independent data, or pivot outside the endpoint/coefficient/source-ratio family; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
