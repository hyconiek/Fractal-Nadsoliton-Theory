#!/usr/bin/env python3
"""P2980/S1930: Z12 CRT idempotent-decomposition source candidate.

P2979 closed the nilradical lane as bounded no-go on current artifacts and
required a genuinely new strict typed object outside that lane.  This audit
supplies one new finite ring-theoretic object: the Chinese-remainder idempotent
decomposition of Z/12Z into orthogonal projectors for the mod-3 and mod-4
factors.

The object is computationally exact: idempotents are enumerated directly, the
nontrivial pair 4 and 9 is orthogonal, and 4+9=1 mod 12.  This is new finite
carrier structure, not a strict source theorem: no current artifact exports
strict nadsoliton provenance, a nonpremise choice of factor semantics, coupling
to damping/selector/bridge/action atoms, or unit-bearing nonproxy L_total
installation.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from math import gcd
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2979_s1929_nilradical_action_variational_installation_obstruction import OUT as P2979

OUT = GEN / "p2980_s1930_z12_crt_idempotent_decomposition_source_candidate.json"
MD = GEN / "p2980_s1930_z12_crt_idempotent_decomposition_source_candidate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MODULUS = 12
FACTORS = [3, 4]


def residues(x: int) -> dict[str, int]:
    return {f"mod_{m}": x % m for m in FACTORS}


def idempotent_certificate() -> dict[str, Any]:
    idempotents = [x for x in range(MODULUS) if (x * x - x) % MODULUS == 0]
    nontrivial = [x for x in idempotents if x not in (0, 1)]
    orthogonal_pairs = []
    for a, b in product(nontrivial, repeat=2):
        if a < b and (a * b) % MODULUS == 0 and (a + b) % MODULUS == 1:
            orthogonal_pairs.append({"a": a, "b": b, "product_mod_12": 0, "sum_mod_12": 1, "residues": {str(a): residues(a), str(b): residues(b)}})
    multiplication_table = {f"{a}*{b}": (a * b) % MODULUS for a in idempotents for b in idempotents}
    return {
        "modulus": MODULUS,
        "factors": FACTORS,
        "coprime_factor_pair": gcd(*FACTORS) == 1,
        "idempotents": idempotents,
        "nontrivial_idempotents": nontrivial,
        "idempotent_count": len(idempotents),
        "orthogonal_completion_pairs": orthogonal_pairs,
        "orthogonal_completion_pair_count": len(orthogonal_pairs),
        "multiplication_table_on_idempotents": multiplication_table,
        "crt_projector_interpretation": {"4": "1 mod 3, 0 mod 4", "9": "0 mod 3, 1 mod 4"},
    }


def candidate_rows(cert: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "candidate": "Z12_CRT_idempotent_decomposition_object",
            "new_outside_nilradical_lane": True,
            "finite_idempotent_witness": True,
            "orthogonal_projector_pair": cert["orthogonal_completion_pair_count"] == 1,
            "strict_nadsoliton_provenance_exported": False,
            "nonpremise_factor_semantics_exported": False,
            "couples_to_named_source_atom": False,
            "action_variational_installation": False,
            "accepted_current_strict_source_object": False,
            "witness": "idempotents [0,1,4,9] with orthogonal completion 4+9=1 and 4*9=0",
        },
        {
            "candidate": "mod3_mod4_projector_semantics",
            "new_outside_nilradical_lane": True,
            "finite_idempotent_witness": True,
            "orthogonal_projector_pair": True,
            "strict_nadsoliton_provenance_exported": False,
            "nonpremise_factor_semantics_exported": False,
            "couples_to_named_source_atom": False,
            "action_variational_installation": False,
            "accepted_current_strict_source_object": False,
            "witness": "4 and 9 split the CRT factors, but factor labels are not physical source semantics",
        },
        {
            "candidate": "completed_strict_CRT_idempotent_source_schema",
            "new_outside_nilradical_lane": True,
            "finite_idempotent_witness": True,
            "orthogonal_projector_pair": True,
            "strict_nadsoliton_provenance_exported": True,
            "nonpremise_factor_semantics_exported": True,
            "couples_to_named_source_atom": True,
            "action_variational_installation": True,
            "accepted_current_strict_source_object": False,
            "witness": "schema row only; current artifacts do not export provenance/coupling/installation",
        },
    ]


def obligation_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    current = [r for r in rows if r["candidate"] != "completed_strict_CRT_idempotent_source_schema"]
    return [
        {"obligation": "new_outside_nilradical_lane", "satisfied": any(r["new_outside_nilradical_lane"] for r in current), "evidence": "CRT idempotent projectors are orthogonal ring projectors, not nilpotents or nilradical receivers"},
        {"obligation": "finite_idempotent_witness", "satisfied": any(r["finite_idempotent_witness"] for r in current), "evidence": "all x in Z12 are enumerated for x^2=x"},
        {"obligation": "orthogonal_projector_pair", "satisfied": any(r["orthogonal_projector_pair"] for r in current), "evidence": "4*9=0 and 4+9=1 mod 12"},
        {"obligation": "strict_nadsoliton_provenance_exported", "satisfied": any(r["strict_nadsoliton_provenance_exported"] for r in current), "evidence": "no theorem maps nadsoliton ontology to this CRT split"},
        {"obligation": "nonpremise_factor_semantics_exported", "satisfied": any(r["nonpremise_factor_semantics_exported"] for r in current), "evidence": "mod3/mod4 labels remain algebraic factor labels, not strict physical semantics"},
        {"obligation": "couples_to_named_source_atom", "satisfied": any(r["couples_to_named_source_atom"] for r in current), "evidence": "no selector, damping, bridge, or action source atom is coupled"},
        {"obligation": "action_variational_installation", "satisfied": any(r["action_variational_installation"] for r in current), "evidence": "no unit-bearing density or nonproxy variational chain is installed"},
        {"obligation": "accepted_current_strict_source_object", "satisfied": any(r["accepted_current_strict_source_object"] for r in current), "evidence": "finite projector readiness is not strict sourcehood"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["new_object", "finite_witness", "orthogonal_pair", "strict_provenance", "factor_semantics", "named_coupling", "action_installation"]
    return [{"present": dict(zip(names, bits)), "accepts_CRT_idempotent_strict_source": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p2979_path: Any) -> dict[str, Any]:
    cert = idempotent_certificate()
    rows = candidate_rows(cert)
    obligations = obligation_rows(rows)
    matrix = acceptance_matrix()
    return {
        "status": "P2980_Z12_CRT_IDEMPOTENT_DECOMPOSITION_SOURCE_CANDIDATE_DEVELOPMENTAL_NO_STRICT_EXPORT",
        "input_hashes": {"P2979": hashlib.sha256(p2979_path.read_bytes()).hexdigest() if p2979_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "Z12CRTIdempotentDecomposition_SourceCandidate",
            "idempotent_certificate": cert,
            "candidate_rows": rows,
            "proof_obligation_rows": obligations,
            "finite_acceptance_matrix": matrix,
        },
        "idempotent_summary": {
            "idempotents": cert["idempotents"],
            "nontrivial_idempotents": cert["nontrivial_idempotents"],
            "orthogonal_completion_pair_count": cert["orthogonal_completion_pair_count"],
            "orthogonal_completion_pairs": cert["orthogonal_completion_pairs"],
            "accepted_current_strict_source_objects": [r["candidate"] for r in rows if r["accepted_current_strict_source_object"]],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_CRT_idempotent_strict_source"]),
        },
        "decision": {
            "positive_progress": "P2980 supplies one genuinely new finite typed object outside the closed nilradical lane: the Z12 CRT idempotent decomposition with projectors 4 and 9.",
            "breakthrough": "Developmental only: the orthogonal projector split is exact, but current artifacts export no strict nadsoliton provenance, nonpremise factor semantics, named source-atom coupling, or action/variational installation.",
            "negative_export_flags": {k: False for k in ["CRT_idempotent_strict_source_exported", "strict_nadsoliton_provenance_exported", "factor_semantics_exported", "named_coupling_exported", "action_density_exported", "nonproxy_ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not promote CRT idempotent projectors, mod3/mod4 factor labels, or orthogonal decompositions to selector closure, damping source, bridge completion, role transfer, nonproxy L_total, or ToE.  A next proof-grade move may attack exactly one missing theorem for this object: strict nadsoliton provenance, nonpremise factor semantics, named source-atom coupling, or action/variational installation; otherwise preserve the P2929-P2980 developmental/no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["idempotent_summary"]
    lines = [
        "# P2980/S1930 Z12 CRT idempotent-decomposition source candidate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Idempotent certificate",
        f"- idempotents: `{cert['idempotents']}`",
        f"- nontrivial idempotents: `{cert['nontrivial_idempotents']}`",
        f"- orthogonal completion pair count: `{cert['orthogonal_completion_pair_count']}`",
        f"- accepted current strict source objects: `{cert['accepted_current_strict_source_objects']}`",
        f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`",
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
    read_json(P2979)
    payload = build_payload(P2979)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2980/S1930 Z12 CRT idempotent-decomposition source candidate", "## P2980/S1930 Z12 CRT idempotent-decomposition source candidate\n\n`P2980/S1930` follows the P2979 nilradical no-go by supplying one genuinely new finite typed object outside that lane: the `Z/12Z` CRT idempotent decomposition.  Exhaustive enumeration gives idempotents `[0,1,4,9]`; the nontrivial projectors `4` and `9` satisfy `4*9=0` and `4+9=1 mod 12`, corresponding algebraically to the mod-3 and mod-4 CRT factors.  This is finite carrier progress only: no strict nadsoliton provenance, nonpremise factor semantics, named source-atom coupling, action-density installation, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2980/S1930 CRT idempotent-decomposition `L_total` guard", "## P2980/S1930 CRT idempotent-decomposition `L_total` guard\n\n`P2980/S1930` adds no sourced term to `L_total`.  The CRT idempotents `4` and `9` are finite orthogonal projectors only; without strict provenance, nonpremise factor semantics, named coupling, unit-bearing density, and a nonproxy variational chain, they cannot enter EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current Z12 CRT idempotent-decomposition source-candidate guardrail (P2980/S1930, 2026-06-20)", "## Current Z12 CRT idempotent-decomposition source-candidate guardrail (P2980/S1930, 2026-06-20)\n\n- P2980 supplies one genuinely new finite typed object outside the closed nilradical lane: the `Z/12Z` CRT idempotent decomposition.\n- Exhaustive enumeration gives idempotents `[0,1,4,9]`; the nontrivial projectors `4` and `9` satisfy `4*9=0` and `4+9=1 mod 12`, with a `128`-profile acceptance matrix and only the full profile accepting.\n- This is developmental algebraic structure, not a strict source theorem: no strict nadsoliton provenance, nonpremise factor semantics, named source-atom coupling, action-density installation, nonproxy `L_total`, bridge closure, role transfer, or ToE is exported.\n- Do not promote CRT projectors, mod-3/mod-4 factor labels, or orthogonal decompositions to selector closure, damping source, bridge completion, role transfer, nonproxy `L_total`, or ToE.  A next admissible move may attack exactly one missing theorem for this object, or preserve the P2929-P2980 developmental/no-strict-export boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
