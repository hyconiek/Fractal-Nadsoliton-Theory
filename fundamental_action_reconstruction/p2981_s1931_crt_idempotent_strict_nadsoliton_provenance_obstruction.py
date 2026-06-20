#!/usr/bin/env python3
"""P2981/S1931: CRT idempotent strict-nadsoliton provenance obstruction.

P2980 supplied a new finite typed object outside the closed nilradical lane: the
Z12 CRT idempotent decomposition.  This audit attacks exactly one missing
P2980 theorem route: strict nadsoliton provenance for that CRT projector split.
It does not replay nilradical objects, factor-semantics promotion, source-atom
coupling, action installation, selector closure, bridge completion, role
transfer, or L_total/ToE promotion.

The finite result is positive on algebra and negative on sourcehood: the
nontrivial idempotents 4 and 9 are the unique orthogonal completion pair of 1
inside Z/12Z, and the mod-3/mod-4 factors are non-isomorphic by size.  However,
current artifacts do not export a theorem that the nadsoliton sources this CRT
split, nor a nonpremise physical meaning for the factor labels.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2980_s1930_z12_crt_idempotent_decomposition_source_candidate import OUT as P2980

OUT = GEN / "p2981_s1931_crt_idempotent_strict_nadsoliton_provenance_obstruction.json"
MD = GEN / "p2981_s1931_crt_idempotent_strict_nadsoliton_provenance_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
MODULUS = 12


def crt_provenance_witness() -> dict[str, Any]:
    idempotents = [x for x in range(MODULUS) if (x * x - x) % MODULUS == 0]
    pairs = []
    for a, b in product(idempotents, repeat=2):
        if a < b and a not in (0, 1) and b not in (0, 1):
            pairs.append({
                "a": a,
                "b": b,
                "orthogonal": (a * b) % MODULUS == 0,
                "sum_is_one": (a + b) % MODULUS == 1,
                "residues": {str(a): {"mod_3": a % 3, "mod_4": a % 4}, str(b): {"mod_3": b % 3, "mod_4": b % 4}},
                "completes_one": (a * b) % MODULUS == 0 and (a + b) % MODULUS == 1,
            })
    completion_pairs = [p for p in pairs if p["completes_one"]]
    return {
        "modulus": MODULUS,
        "idempotents": idempotents,
        "nontrivial_idempotents": [x for x in idempotents if x not in (0, 1)],
        "nontrivial_pair_rows": pairs,
        "orthogonal_completion_pairs": completion_pairs,
        "unique_nontrivial_completion_pair": len(completion_pairs) == 1,
        "factor_cardinalities": {"mod_3": 3, "mod_4": 4},
        "factor_cardinality_distinguishes_labels": 3 != 4,
        "strict_nadsoliton_source_map_exported": False,
        "nonpremise_factor_semantics_exported": False,
    }


def candidate_rows(witness: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "candidate": "finite_CRT_projector_uniqueness",
            "finite_projector_uniqueness": witness["unique_nontrivial_completion_pair"],
            "factor_label_distinction": witness["factor_cardinality_distinguishes_labels"],
            "strict_nadsoliton_source_map_exported": False,
            "nonpremise_factor_semantics_exported": False,
            "couples_to_named_source_atom": False,
            "accepted_strict_provenance_theorem": False,
            "witness": "4 and 9 are the unique nontrivial orthogonal idempotents summing to 1",
        },
        {
            "candidate": "cardinality_distinguished_mod3_mod4_labels",
            "finite_projector_uniqueness": True,
            "factor_label_distinction": True,
            "strict_nadsoliton_source_map_exported": False,
            "nonpremise_factor_semantics_exported": False,
            "couples_to_named_source_atom": False,
            "accepted_strict_provenance_theorem": False,
            "witness": "mod3 and mod4 factors are algebraically distinguishable by size, but not physically sourced",
        },
        {
            "candidate": "completed_strict_CRT_projector_provenance_schema",
            "finite_projector_uniqueness": True,
            "factor_label_distinction": True,
            "strict_nadsoliton_source_map_exported": True,
            "nonpremise_factor_semantics_exported": True,
            "couples_to_named_source_atom": True,
            "accepted_strict_provenance_theorem": False,
            "witness": "schema row only; current artifacts do not export strict provenance or semantics",
        },
    ]


def obligation_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    current = [r for r in rows if r["candidate"] != "completed_strict_CRT_projector_provenance_schema"]
    return [
        {"obligation": "finite_projector_uniqueness", "satisfied": any(r["finite_projector_uniqueness"] for r in current), "evidence": "the only nontrivial orthogonal completion pair is (4,9)"},
        {"obligation": "factor_label_distinction", "satisfied": any(r["factor_label_distinction"] for r in current), "evidence": "mod3 and mod4 factors have different cardinalities"},
        {"obligation": "strict_nadsoliton_source_map_exported", "satisfied": any(r["strict_nadsoliton_source_map_exported"] for r in current), "evidence": "no theorem maps nadsoliton ontology to this CRT projector split"},
        {"obligation": "nonpremise_factor_semantics_exported", "satisfied": any(r["nonpremise_factor_semantics_exported"] for r in current), "evidence": "algebraic mod3/mod4 labels are not a strict physical semantics theorem"},
        {"obligation": "couples_to_named_source_atom", "satisfied": any(r["couples_to_named_source_atom"] for r in current), "evidence": "no selector, damping, bridge, or action atom is coupled"},
        {"obligation": "accepted_strict_provenance_theorem", "satisfied": any(r["accepted_strict_provenance_theorem"] for r in current), "evidence": "finite uniqueness is readiness, not sourcehood"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_uniqueness", "factor_distinction", "strict_source_map", "factor_semantics", "named_coupling", "nonpromotion_boundary"]
    return [{"present": dict(zip(names, bits)), "accepts_CRT_strict_provenance": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p2980_path: Any) -> dict[str, Any]:
    witness = crt_provenance_witness()
    rows = candidate_rows(witness)
    obligations = obligation_rows(rows)
    matrix = acceptance_matrix()
    return {
        "status": "P2981_CRT_IDEMPOTENT_STRICT_NADSOLITON_PROVENANCE_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P2980": hashlib.sha256(p2980_path.read_bytes()).hexdigest() if p2980_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "CRTIdempotentStrictNadsolitonProvenance_ObstructionMatrix",
            "crt_provenance_witness": witness,
            "candidate_rows": rows,
            "proof_obligation_rows": obligations,
            "finite_acceptance_matrix": matrix,
        },
        "provenance_certificate": {
            "idempotents": witness["idempotents"],
            "orthogonal_completion_pairs": witness["orthogonal_completion_pairs"],
            "unique_nontrivial_completion_pair": witness["unique_nontrivial_completion_pair"],
            "factor_cardinality_distinguishes_labels": witness["factor_cardinality_distinguishes_labels"],
            "accepted_current_strict_provenance_theorems": [r["candidate"] for r in rows if r["accepted_strict_provenance_theorem"]],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_CRT_strict_provenance"]),
        },
        "decision": {
            "positive_progress": "P2981 attacks exactly one P2980 missing theorem: strict nadsoliton provenance for the CRT idempotent projector split.",
            "breakthrough": "Bounded no-go: finite uniqueness and factor cardinality distinction are real, but no current artifact exports a strict nadsoliton source map, nonpremise factor semantics, or named coupling for the CRT split.",
            "negative_export_flags": {k: False for k in ["CRT_strict_provenance_exported", "factor_semantics_exported", "named_coupling_exported", "action_density_exported", "nonproxy_ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay CRT projector uniqueness, cardinality-distinguished factor labels, nilradical lanes, selector replay, bridge maps, role transfer, or L_total placeholders as provenance.  A next proof-grade move may attack exactly one remaining CRT route: nonpremise factor semantics, named source-atom coupling, or action/variational installation; otherwise preserve the P2929-P2981 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["provenance_certificate"]
    lines = [
        "# P2981/S1931 CRT idempotent strict-nadsoliton provenance obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Provenance certificate",
        f"- idempotents: `{cert['idempotents']}`",
        f"- orthogonal completion pairs: `{cert['orthogonal_completion_pairs']}`",
        f"- unique nontrivial completion pair: `{cert['unique_nontrivial_completion_pair']}`",
        f"- factor cardinality distinguishes labels: `{cert['factor_cardinality_distinguishes_labels']}`",
        f"- accepted current strict provenance theorems: `{cert['accepted_current_strict_provenance_theorems']}`",
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
    read_json(P2980)
    payload = build_payload(P2980)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2981/S1931 CRT idempotent strict-nadsoliton provenance obstruction", "## P2981/S1931 CRT idempotent strict-nadsoliton provenance obstruction\n\n`P2981/S1931` attacks exactly one P2980 missing theorem: strict nadsoliton provenance for the CRT idempotent projector split.  The finite side is positive: `[0,1,4,9]` are the idempotents, `(4,9)` is the unique nontrivial orthogonal completion pair of `1`, and the mod-3/mod-4 factors are cardinality-distinguished.  The source side remains blocked: no current artifact exports a strict nadsoliton source map, nonpremise factor semantics, named source-atom coupling, action-density installation, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2981/S1931 CRT provenance `L_total` guard", "## P2981/S1931 CRT provenance `L_total` guard\n\n`P2981/S1931` adds no CRT projector term to `L_total`.  Finite idempotent uniqueness and factor cardinality distinction are provenance-readiness only; without a strict source map, nonpremise factor semantics, named coupling, unit-bearing density, and nonproxy variational chain, EOM, Hamiltonian, bridge closure, role transfer, and ToE remain unexported.\n")
    append_once(AGENTS, "Current CRT idempotent strict-nadsoliton provenance obstruction guardrail (P2981/S1931, 2026-06-20)", "## Current CRT idempotent strict-nadsoliton provenance obstruction guardrail (P2981/S1931, 2026-06-20)\n\n- P2981 attacks exactly one P2980 missing theorem: strict nadsoliton provenance for the `Z/12Z` CRT idempotent projector split.\n- Finite positives remain real: `[0,1,4,9]` are the idempotents, `(4,9)` is the unique nontrivial orthogonal completion pair of `1`, and mod-3/mod-4 factors are cardinality-distinguished; the acceptance matrix has `64` profiles with only the full profile accepting.\n- The current route is bounded no-go: no strict nadsoliton source map, nonpremise factor semantics, or named source-atom coupling is exported.\n- Do not promote CRT projector uniqueness, factor cardinality labels, nilradical replay, selector replay, bridge maps, role transfer, nonproxy `L_total`, or ToE.  A next admissible move may attack one remaining CRT route (nonpremise factor semantics, named source-atom coupling, or action/variational installation) or preserve the P2929-P2981 no-strict-export boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
