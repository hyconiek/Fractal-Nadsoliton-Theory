#!/usr/bin/env python3
"""P2982/S1932: CRT idempotent nonpremise factor-semantics obstruction.

After P2981, this audit attacks exactly one remaining CRT route: nonpremise
factor semantics for the Z12 CRT idempotent projector split. It does not replay
strict provenance, nilradical objects, named source-atom coupling, action
installation, selector closure, bridge completion, role transfer, or
L_total/ToE promotion.

The finite result is positive on algebra and negative on semantics: the
projectors 4 and 9 implement exact CRT projections to mod-3 and mod-4 factors,
and the factors are cardinality-distinguished. However, cardinality and residue
labels do not by themselves export a nonpremise strict physical semantics, a
named source atom, or a unit-bearing action installation.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2981_s1931_crt_idempotent_strict_nadsoliton_provenance_obstruction import OUT as P2981

OUT = GEN / "p2982_s1932_crt_factor_semantics_obstruction.json"
MD = GEN / "p2982_s1932_crt_factor_semantics_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
MODULUS = 12


def crt_factor_semantics_witness() -> dict[str, Any]:
    idempotents = [x for x in range(MODULUS) if (x * x - x) % MODULUS == 0]
    projectors = {4: "mod_3_identity_mod_4_zero", 9: "mod_3_zero_mod_4_identity"}
    rows = []
    for projector, label in projectors.items():
        image = sorted({(projector * x) % MODULUS for x in range(MODULUS)})
        kernel = [x for x in range(MODULUS) if (projector * x) % MODULUS == 0]
        rows.append({
            "projector": projector,
            "algebraic_label": label,
            "residue_signature": {"mod_3": projector % 3, "mod_4": projector % 4},
            "image": image,
            "image_size": len(image),
            "kernel": kernel,
            "kernel_size": len(kernel),
            "idempotent": (projector * projector - projector) % MODULUS == 0,
        })
    return {
        "modulus": MODULUS,
        "idempotents": idempotents,
        "nontrivial_projectors": sorted(projectors),
        "orthogonal_completion_pair": {"a": 4, "b": 9, "product_mod_12": (4 * 9) % MODULUS, "sum_mod_12": (4 + 9) % MODULUS},
        "projector_semantics_rows": rows,
        "factor_cardinalities": {"mod_3_factor": 3, "mod_4_factor": 4},
        "factor_cardinality_distinguishes_labels": True,
        "algebraic_factor_semantics_present": True,
        "nonpremise_physical_semantics_exported": False,
        "strict_nadsoliton_source_map_exported": False,
        "named_source_atom_coupling_exported": False,
    }


def candidate_rows(witness: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "candidate": "CRT_algebraic_factor_labels",
            "exact_projector_split": witness["orthogonal_completion_pair"]["product_mod_12"] == 0 and witness["orthogonal_completion_pair"]["sum_mod_12"] == 1,
            "algebraic_factor_semantics_present": True,
            "nonpremise_physical_semantics_exported": False,
            "strict_nadsoliton_source_map_exported": False,
            "couples_to_named_source_atom": False,
            "accepted_factor_semantics_theorem": False,
            "witness": "4 projects to the mod-3 side and 9 projects to the mod-4 side, but only as ring algebra",
        },
        {
            "candidate": "cardinality_distinguished_factor_meaning",
            "exact_projector_split": True,
            "algebraic_factor_semantics_present": True,
            "nonpremise_physical_semantics_exported": False,
            "strict_nadsoliton_source_map_exported": False,
            "couples_to_named_source_atom": False,
            "accepted_factor_semantics_theorem": False,
            "witness": "factor sizes 3 and 4 distinguish labels but do not name strict physical roles",
        },
        {
            "candidate": "completed_nonpremise_factor_semantics_schema",
            "exact_projector_split": True,
            "algebraic_factor_semantics_present": True,
            "nonpremise_physical_semantics_exported": True,
            "strict_nadsoliton_source_map_exported": True,
            "couples_to_named_source_atom": True,
            "accepted_factor_semantics_theorem": False,
            "witness": "schema row only; current artifacts do not export strict physical semantics or named coupling",
        },
    ]


def obligation_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    current = [r for r in rows if r["candidate"] != "completed_nonpremise_factor_semantics_schema"]
    return [
        {"obligation": "exact_projector_split", "satisfied": any(r["exact_projector_split"] for r in current), "evidence": "4*9=0 and 4+9=1 mod 12"},
        {"obligation": "algebraic_factor_semantics_present", "satisfied": any(r["algebraic_factor_semantics_present"] for r in current), "evidence": "residue signatures identify mod-3/mod-4 CRT factors"},
        {"obligation": "nonpremise_physical_semantics_exported", "satisfied": any(r["nonpremise_physical_semantics_exported"] for r in current), "evidence": "no theorem turns mod-3/mod-4 labels into strict physical roles"},
        {"obligation": "strict_nadsoliton_source_map_exported", "satisfied": any(r["strict_nadsoliton_source_map_exported"] for r in current), "evidence": "P2981 already blocks strict source provenance"},
        {"obligation": "couples_to_named_source_atom", "satisfied": any(r["couples_to_named_source_atom"] for r in current), "evidence": "no selector, damping, bridge, or action atom is coupled"},
        {"obligation": "accepted_factor_semantics_theorem", "satisfied": any(r["accepted_factor_semantics_theorem"] for r in current), "evidence": "algebraic semantics is readiness, not nonpremise physical semantics"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["exact_projector_split", "algebraic_factor_semantics", "nonpremise_physical_semantics", "strict_source_map", "named_coupling", "nonpromotion_boundary"]
    return [{"present": dict(zip(names, bits)), "accepts_CRT_nonpremise_factor_semantics": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p2981_path: Any) -> dict[str, Any]:
    witness = crt_factor_semantics_witness()
    rows = candidate_rows(witness)
    obligations = obligation_rows(rows)
    matrix = acceptance_matrix()
    return {
        "status": "P2982_CRT_IDEMPOTENT_NONPREMISE_FACTOR_SEMANTICS_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P2981": hashlib.sha256(p2981_path.read_bytes()).hexdigest() if p2981_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "CRTIdempotentNonpremiseFactorSemantics_ObstructionMatrix",
            "crt_factor_semantics_witness": witness,
            "candidate_rows": rows,
            "proof_obligation_rows": obligations,
            "finite_acceptance_matrix": matrix,
        },
        "factor_semantics_certificate": {
            "idempotents": witness["idempotents"],
            "projector_semantics_rows": witness["projector_semantics_rows"],
            "orthogonal_completion_pair": witness["orthogonal_completion_pair"],
            "algebraic_factor_semantics_present": witness["algebraic_factor_semantics_present"],
            "accepted_current_factor_semantics_theorems": [r["candidate"] for r in rows if r["accepted_factor_semantics_theorem"]],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_CRT_nonpremise_factor_semantics"]),
        },
        "decision": {
            "positive_progress": "P2982 attacks exactly one P2981 remaining CRT route: nonpremise factor semantics for the CRT idempotent projector split.",
            "breakthrough": "Bounded no-go: the mod-3/mod-4 algebraic factor labels are exact, but current artifacts do not export a nonpremise strict physical semantics, strict source map, or named source-atom coupling for them.",
            "negative_export_flags": {k: False for k in ["nonpremise_factor_semantics_exported", "CRT_strict_provenance_exported", "named_coupling_exported", "action_density_exported", "nonproxy_ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay CRT projector uniqueness, cardinality labels, residue signatures, nilradical lanes, selector replay, bridge maps, role transfer, or L_total placeholders as semantics.  A next proof-grade move may attack exactly one remaining CRT route: named source-atom coupling or action/variational installation; otherwise preserve the P2929-P2982 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["factor_semantics_certificate"]
    lines = [
        "# P2982/S1932 CRT idempotent nonpremise factor-semantics obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Factor semantics certificate",
        f"- idempotents: `{cert['idempotents']}`",
        f"- orthogonal completion pair: `{cert['orthogonal_completion_pair']}`",
        f"- projector semantics rows: `{cert['projector_semantics_rows']}`",
        f"- algebraic factor semantics present: `{cert['algebraic_factor_semantics_present']}`",
        f"- accepted current factor semantics theorems: `{cert['accepted_current_factor_semantics_theorems']}`",
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
    read_json(P2981)
    payload = build_payload(P2981)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2982/S1932 CRT idempotent nonpremise factor-semantics obstruction", "## P2982/S1932 CRT idempotent nonpremise factor-semantics obstruction\n\n`P2982/S1932` attacks exactly one P2981 remaining CRT route: nonpremise factor semantics for the CRT idempotent projector split.  The finite side is positive: projectors `4` and `9` implement exact mod-3/mod-4 CRT projections, satisfy `4*9=0` and `4+9=1 mod 12`, and give algebraic factor labels.  The source/semantics side remains blocked: no current artifact exports nonpremise strict physical semantics, strict nadsoliton source map, named source-atom coupling, action-density installation, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2982/S1932 CRT factor-semantics `L_total` guard", "## P2982/S1932 CRT factor-semantics `L_total` guard\n\n`P2982/S1932` adds no CRT factor-semantics term to `L_total`.  Residue signatures and cardinality-distinguished mod-3/mod-4 labels are algebraic semantics only; without nonpremise physical semantics, strict source map, named coupling, unit-bearing density, and nonproxy variational chain, EOM, Hamiltonian, bridge closure, role transfer, and ToE remain unexported.\n")
    append_once(AGENTS, "Current CRT idempotent nonpremise factor-semantics obstruction guardrail (P2982/S1932, 2026-06-20)", "## Current CRT idempotent nonpremise factor-semantics obstruction guardrail (P2982/S1932, 2026-06-20)\n\n- P2982 attacks exactly one P2981 remaining CRT route: nonpremise factor semantics for the `Z/12Z` CRT idempotent projector split.\n- Finite positives remain real: projectors `4` and `9` implement exact mod-3/mod-4 CRT projections, satisfy `4*9=0` and `4+9=1 mod 12`, and the acceptance matrix has `64` profiles with only the full profile accepting.\n- The current route is bounded no-go: residue signatures and cardinality labels are algebraic semantics only; no nonpremise strict physical semantics, strict source map, or named source-atom coupling is exported.\n- Do not promote CRT projector uniqueness, cardinality labels, residue signatures, nilradical replay, selector replay, bridge maps, role transfer, nonproxy `L_total`, or ToE.  A next admissible move may attack one remaining CRT route (named source-atom coupling or action/variational installation) or preserve the P2929-P2982 no-strict-export boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
