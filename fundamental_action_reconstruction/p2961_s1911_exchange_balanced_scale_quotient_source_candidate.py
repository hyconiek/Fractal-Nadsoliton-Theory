#!/usr/bin/env python3
"""P2961/S1911: exchange-balanced scale quotient source candidate.

P2960 admitted developmental localizer laws but did not prove one.  P2961 takes
one concrete step: construct the canonical scale quotient for positive integer
provenance weights (a,b), quotient by common positive scale gcd(a,b), then test
the involution exchanging the two provenance sectors K and C.

The finite computation is backed by an unbounded lemma: in primitive positive
integer pairs, the only fixed point of exchange (a,b)->(b,a) is (1,1), because
a=b and gcd(a,b)=1 imply a=b=1.  Thus the developmental self-balance/minimal
source candidate is sharpened from a convention into a precise quotient-fixed
source candidate.  It still is not a strict export until the theory proves that
the nadsoliton really supplies this exchange symmetry as a source law and until
a unit-bearing nonproxy coupling is constructed.
"""
from __future__ import annotations

import hashlib, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2959_s1909_p2938_u12_aggregate_localizer_acceptance_no_go import K, C, TARGET
from p2960_s1910_developmental_ontology_localizer_law_intake import OUT as P2960

OUT = GEN / "p2961_s1911_exchange_balanced_scale_quotient_source_candidate.json"
MD = GEN / "p2961_s1911_exchange_balanced_scale_quotient_source_candidate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def vector(a: int, b: int) -> list[int]:
    return [a * x + b * y for x, y in zip(K, C)]


def primitive_pair(a: int, b: int) -> tuple[int, int]:
    g = math.gcd(a, b)
    return (a // g, b // g)


def quotient_rows(max_weight: int = 12) -> list[dict[str, Any]]:
    rows = []
    for a in range(1, max_weight + 1):
        for b in range(1, max_weight + 1):
            p = primitive_pair(a, b)
            orbit = sorted([list(p), [p[1], p[0]]])
            rows.append({
                "a": a,
                "b": b,
                "gcd_scale": math.gcd(a, b),
                "primitive_representative": list(p),
                "exchange_image": [p[1], p[0]],
                "exchange_orbit": orbit,
                "is_primitive": math.gcd(a, b) == 1,
                "is_exchange_fixed_after_quotient": p == (p[1], p[0]),
                "vector": vector(*p),
                "equals_target_after_quotient": vector(*p) == TARGET,
            })
    return rows


def quotient_orbits(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    seen = {}
    for row in rows:
        key = tuple(tuple(x) for x in row["exchange_orbit"])
        seen.setdefault(key, row)
    out = []
    for key, row in sorted(seen.items()):
        rep = tuple(row["primitive_representative"])
        fixed = row["is_exchange_fixed_after_quotient"]
        out.append({
            "orbit": [list(x) for x in key],
            "representative": list(rep),
            "is_fixed_orbit": fixed,
            "selected_by_exchange_balance": fixed,
            "vector": vector(*rep),
            "sum": sum(vector(*rep)),
            "equals_target": vector(*rep) == TARGET,
        })
    return out


def proof_obligation_rows(orbits: list[dict[str, Any]]) -> list[dict[str, Any]]:
    fixed = [o for o in orbits if o["is_fixed_orbit"]]
    return [
        {"obligation": "canonical_gcd_scale_quotient_constructed", "satisfied": True, "evidence": "positive integer pairs are quotiented by gcd(a,b) to primitive representatives"},
        {"obligation": "exchange_action_computed_on_quotient", "satisfied": True, "evidence": "the involution (a,b)->(b,a) is evaluated on every primitive orbit in the bounded audit"},
        {"obligation": "unique_exchange_fixed_primitive_orbit", "satisfied": fixed == [o for o in fixed if o["representative"] == [1, 1]], "evidence": "finite rows agree with the unbounded lemma gcd(a,a)=a, so primitive fixed implies a=1"},
        {"obligation": "target_vector_emerges_from_fixed_orbit", "satisfied": len(fixed) == 1 and fixed[0]["equals_target"] and fixed[0]["sum"] == 9, "evidence": "the fixed primitive orbit gives K+C=[1,2,2,2,2] and sum 9"},
        {"obligation": "strict_nadsoliton_exchange_symmetry_source_exported", "satisfied": False, "evidence": "P2961 constructs the candidate source law but does not prove that the nadsoliton exports K/C exchange symmetry"},
        {"obligation": "unit_bearing_nonproxy_coupling_exported", "satisfied": False, "evidence": "no action-density coupling or L_total term is constructed"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["canonical_gcd_scale_quotient_constructed", "exchange_action_computed_on_quotient", "unique_exchange_fixed_primitive_orbit", "target_vector_emerges_from_fixed_orbit", "strict_nadsoliton_exchange_symmetry_source_exported", "unit_bearing_nonproxy_coupling_exported"]
    rows = []
    for mask in range(1 << len(names)):
        present = {name: bool(mask & (1 << i)) for i, name in enumerate(names)}
        rows.append({
            "mask": mask,
            "present": present,
            "accepts_as_developmental_source_theorem_candidate": all(present[n] for n in names[:4]),
            "accepts_as_strict_ratio_package_source": all(present.values()),
        })
    return rows


def build_payload(p2960: dict[str, Any]) -> dict[str, Any]:
    qrows = quotient_rows()
    orbits = quotient_orbits(qrows)
    obligations = proof_obligation_rows(orbits)
    matrix = acceptance_matrix()
    fixed = [o for o in orbits if o["is_fixed_orbit"]]
    return {
        "status": "P2961_EXCHANGE_BALANCED_SCALE_QUOTIENT_SOURCE_CANDIDATE_NO_STRICT_EXPORT",
        "input_hashes": {"P2960": hashlib.sha256(P2960.read_bytes()).hexdigest() if P2960.exists() else None},
        "constructed_theoretical_objects": {
            "candidate_object": "ExchangeBalanced_PrimitiveScaleQuotient_SourceCandidate",
            "unbounded_lemma": "For positive integers a,b, the exchange-fixed primitive quotient classes satisfy a=b and gcd(a,b)=1, hence (a,b)=(1,1).",
            "quotient_rows": qrows,
            "exchange_orbit_rows": orbits,
            "proof_obligation_rows": obligations,
            "finite_acceptance_matrix": matrix,
        },
        "quotient_certificate": {
            "bounded_pair_rows": len(qrows),
            "exchange_orbit_count": len(orbits),
            "fixed_orbits": fixed,
            "unique_fixed_orbit_is_target": len(fixed) == 1 and fixed[0]["representative"] == [1, 1] and fixed[0]["equals_target"],
            "target_sum_from_fixed_orbit": fixed[0]["sum"] if fixed else None,
            "developmental_source_candidate_exported": all(row["satisfied"] for row in obligations[:4]),
            "strict_nadsoliton_exchange_symmetry_source_exported": False,
            "unit_bearing_nonproxy_coupling_exported": False,
            "strict_ratio_package_source_exported": all(row["satisfied"] for row in obligations),
            "acceptance_matrix_rows": len(matrix),
            "developmental_acceptance_rows": sum(1 for row in matrix if row["accepts_as_developmental_source_theorem_candidate"]),
            "strict_acceptance_rows": sum(1 for row in matrix if row["accepts_as_strict_ratio_package_source"]),
        },
        "decision": {
            "positive_progress": "P2961 upgrades one P2960 candidate: the scale quotient and exchange-fixed primitive orbit are now explicit, and the unique fixed orbit produces the exact P2938 vector and sum 9 without using a target_sum=9 cut.",
            "breakthrough": "This is a real source-candidate theorem, not full strict closure.  The missing step is to prove that the nadsoliton itself exports the K/C exchange symmetry and then couple the selected quotient to a unit-bearing nonproxy action density.",
            "negative_export_flags": {k: False for k in ["strict_nadsoliton_exchange_symmetry_source_exported", "strict_p2938_provenance_exported", "strict_ratio_package_source_exported", "damping_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay exchange-fixed quotient enumeration, target_sum=9 cuts, primitive equal-summand convention, bounded localizer variants, K+C decompositions, beta-scale normalization, scalar Euler insertion, or count/role aliases.  The next proof-grade move must either prove a strict nadsoliton source theorem for the K/C exchange symmetry used by P2961, or construct the unit-bearing nonproxy action-density coupling receiving the P2961 quotient-selected vector; otherwise pivot outside the ratio-package lane while preserving the P2929-P2961 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["quotient_certificate"]
    lines = [
        "# P2961/S1911 exchange-balanced scale quotient source candidate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Quotient certificate",
        f"- bounded pair rows: `{cert['bounded_pair_rows']}`",
        f"- exchange orbit count: `{cert['exchange_orbit_count']}`",
        f"- fixed orbits: `{cert['fixed_orbits']}`",
        f"- unique fixed orbit is target: `{cert['unique_fixed_orbit_is_target']}`",
        f"- target sum from fixed orbit: `{cert['target_sum_from_fixed_orbit']}`",
        f"- developmental source candidate exported: `{cert['developmental_source_candidate_exported']}`",
        f"- strict nadsoliton exchange-symmetry source exported: `{cert['strict_nadsoliton_exchange_symmetry_source_exported']}`",
        f"- unit-bearing nonproxy coupling exported: `{cert['unit_bearing_nonproxy_coupling_exported']}`",
        f"- strict ratio-package source exported: `{cert['strict_ratio_package_source_exported']}`",
        f"- acceptance matrix rows/developmental/strict: `{cert['acceptance_matrix_rows']}/{cert['developmental_acceptance_rows']}/{cert['strict_acceptance_rows']}`",
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
    payload = build_payload(read_json(P2960))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2961/S1911 exchange-balanced scale quotient source candidate", "## P2961/S1911 exchange-balanced scale quotient source candidate\n\n`P2961/S1911` takes one P2960 developmental candidate and constructs its canonical scale quotient.  Positive integer provenance weights `(a,b)` are quotiented by `gcd(a,b)`, then the exchange involution `(a,b)->(b,a)` is evaluated on primitive quotient classes.  The unbounded lemma is immediate: an exchange-fixed primitive class has `a=b` and `gcd(a,b)=1`, hence `(a,b)=(1,1)`.  The unique fixed class produces `K+C=[1,2,2,2,2]` and sum `9` without a target_sum cut.  This exports a developmental source-candidate theorem, not strict provenance: the strict nadsoliton source for the K/C exchange symmetry and the unit-bearing nonproxy coupling remain missing.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2961/S1911 exchange-balanced quotient `L_total` guard", "## P2961/S1911 exchange-balanced quotient `L_total` guard\n\n`P2961/S1911` constructs a canonical primitive scale quotient and exchange-fixed source candidate for the P2938 vector, but no strict nadsoliton exchange-symmetry source theorem or unit-bearing nonproxy action-density coupling is exported.  Therefore the selected vector still cannot enter `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE as a sourced damping term.\n")
    append_once(AGENTS, "Current exchange-balanced scale quotient source-candidate guardrail (P2961/S1911, 2026-06-20)", "## Current exchange-balanced scale quotient source-candidate guardrail (P2961/S1911, 2026-06-20)\n\n- P2961 takes one P2960 developmental candidate and makes it more proof-grade: quotient positive integer provenance weights `(a,b)` by `gcd(a,b)`, then impose the K/C exchange involution on primitive quotient classes.\n- The unbounded lemma is that an exchange-fixed primitive class must be `(1,1)`, so the quotient-selected vector is `K+C=[1,2,2,2,2]` with sum `9` without importing a target_sum cut.\n- This is a developmental source-candidate theorem, not strict closure: current artifacts still do not prove that the nadsoliton exports the K/C exchange symmetry or a unit-bearing nonproxy action-density coupling.\n- Do not promote P2961 to strict P2938 provenance, ratio-package source, damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE.  The next admissible move must prove the strict nadsoliton K/C exchange-symmetry source theorem, construct the unit-bearing nonproxy coupling receiving the P2961 quotient-selected vector, or pivot outside the ratio-package lane while preserving the P2929-P2961 boundary.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
