#!/usr/bin/env python3
"""P2708/S1658: Z12 boundary 1-cocycle selector-source obstruction.

A new typed-object test after P2707: use the oriented boundary 1-cocycle line
of the Z12 cycle as a possible strict symmetry-breaking provider, then compute
whether it canonically selects +1 rather than -1 without an external sign.
"""
from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2708_s1658_z12_boundary_1_cocycle_selector_source_obstruction.json"
MD = GEN / "p2708_s1658_z12_boundary_1_cocycle_selector_source_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2707": GEN / "p2707_s1657_post_p2706_no_new_live_frontier_reconciliation.json",
    "P2706": GEN / "p2706_s1656_damping_to_selector_interface_obstruction_witness_table.json",
    "P2700": GEN / "p2700_s1650_exhaustive_aut_invariant_selector_functional_enumeration_no_go.json",
    "P2699": GEN / "p2699_s1649_z12_fractal_information_aut_invariant_selector_source_no_go.json",
}

N = 12
NEGATIVE_EXPORT_FLAGS = [
    "qw2191_discharged",
    "non_premise_selector_provider_exported",
    "canonical_orientation_selected",
    "pair12_strict_core_upgrade_exported",
    "ltotal_promoted",
    "toe_closure_exported",
]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def incidence_matrix_cycle(n: int = N) -> list[list[int]]:
    # Rows are vertices, columns are oriented edges e_i: i -> i+1.
    mat = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        mat[i][i] = -1
        mat[(i + 1) % n][i] = 1
    return mat


def rank_fraction(matrix: list[list[int]]) -> int:
    a = [[Fraction(x) for x in row] for row in matrix]
    if not a:
        return 0
    rows, cols = len(a), len(a[0])
    rank = 0
    for col in range(cols):
        pivot = None
        for row in range(rank, rows):
            if a[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        pivot_value = a[rank][col]
        a[rank] = [x / pivot_value for x in a[rank]]
        for row in range(rows):
            if row != rank and a[row][col] != 0:
                factor = a[row][col]
                a[row] = [x - factor * y for x, y in zip(a[row], a[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


def cohomology_certificate() -> dict[str, Any]:
    d0 = incidence_matrix_cycle(N)
    rank_d0 = rank_fraction(d0)
    c1_dim = N
    h1_dim = c1_dim - rank_d0
    circulation = [1 for _ in range(N)]
    boundary_pairing = sum(circulation)
    return {
        "chain_model": "Z12 cycle graph with vertices v_i and oriented boundary edges e_i: i -> i+1",
        "vertex_count": N,
        "edge_count": N,
        "rank_d0": rank_d0,
        "h1_dimension": h1_dim,
        "primitive_circulation_vector": circulation,
        "primitive_boundary_pairing": boundary_pairing,
        "cohomology_reading": "There is a one-dimensional oriented circulation line, so a premise sign can orient the circle; the line itself still contains the pair ±omega.",
    }


def automorphism_action_certificate() -> dict[str, Any]:
    omega = [1 for _ in range(N)]
    rotation_action = omega[:]  # i -> i+1 preserves edge orientation.
    inversion_action = [-x for x in reversed(omega)]  # i -> -i reverses the oriented boundary cycle.
    rotation_fixed = rotation_action == omega
    inversion_fixed = inversion_action == omega
    invariant_nonzero_possible = inversion_fixed
    sign_pair = {"plus_orientation_pairing": sum(omega), "minus_orientation_pairing": -sum(omega)}
    return {
        "rotation_fixed": rotation_fixed,
        "inversion_sends_omega_to_minus_omega": inversion_action == [-x for x in omega],
        "aut_z12_invariant_nonzero_orientation_exists": invariant_nonzero_possible,
        "invariant_subspace_dimension_under_inversion": 0 if not invariant_nonzero_possible else 1,
        "sign_pair": sign_pair,
        "selector_reading": "The boundary 1-cocycle line gives a two-point orientation pair.  Aut(Z12) inversion exchanges the two signs, so choosing +omega over -omega is still an extra sign/branch premise unless a new strict source breaks inversion.",
    }


def premise_vs_strict_table(cohomology: dict[str, Any], aut: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "candidate_object": "oriented_boundary_1_cocycle_line_H1_Z12",
            "positive_content": "nonzero H1 line exists and can carry a premise orientation sign",
            "computed_evidence": {"h1_dimension": cohomology["h1_dimension"], "primitive_boundary_pairing": cohomology["primitive_boundary_pairing"]},
            "strict_provider_obligation_met": False,
            "blocker": "H1 supplies ±omega, not a canonical choice of +omega.",
        },
        {
            "candidate_object": "Aut_Z12_invariant_orientation",
            "positive_content": "rotation preserves the circulation vector",
            "computed_evidence": {"rotation_fixed": aut["rotation_fixed"], "inversion_sends_omega_to_minus_omega": aut["inversion_sends_omega_to_minus_omega"], "invariant_dimension": aut["invariant_subspace_dimension_under_inversion"]},
            "strict_provider_obligation_met": False,
            "blocker": "Full Aut(Z12) includes inversion, which kills every nonzero invariant orientation.",
        },
        {
            "candidate_object": "premise_signed_orientation",
            "positive_content": "if an external premise selects +omega, it orients +1 over -1",
            "computed_evidence": aut["sign_pair"],
            "strict_provider_obligation_met": False,
            "blocker": "This is premise-based sign selection, not a non-premise strict-sourced provider.",
        },
    ]


def source_boundary() -> dict[str, Any]:
    p2707 = read_json(INPUTS["P2707"])
    p2706 = read_json(INPUTS["P2706"])
    p2700 = read_json(INPUTS["P2700"])
    p2699 = read_json(INPUTS["P2699"])
    return {
        "p2707_no_new_live_frontier": p2707.get("decision", {}).get("no_new_live_frontier_certificate") is True,
        "p2706_damping_no_selector": p2706.get("decision", {}).get("damping_transport_exports_directed_selector") is False,
        "p2700_bounded_no_go": p2700.get("decision", {}).get("bounded_no_go_now") is True,
        "p2699_bounded_no_go": p2699.get("decision", {}).get("bounded_no_go_now") is True,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2708/S1658 Z12 boundary 1-cocycle selector-source obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Computation",
        f"- rank_d0={payload['cohomology_certificate']['rank_d0']}",
        f"- h1_dimension={payload['cohomology_certificate']['h1_dimension']}",
        f"- invariant_subspace_dimension_under_inversion={payload['automorphism_action_certificate']['invariant_subspace_dimension_under_inversion']}",
        "",
        "## Candidate table",
    ]
    for row in payload["premise_vs_strict_table"]:
        lines.append(f"- `{row['candidate_object']}`: strict_provider_obligation_met={row['strict_provider_obligation_met']}. {row['blocker']}")
    lines.extend(["", "## Decision", payload["decision"]["current_unlock_reading"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    cohomology = cohomology_certificate()
    aut = automorphism_action_certificate()
    table = premise_vs_strict_table(cohomology, aut)
    boundary = source_boundary()
    no_unlock = cohomology["h1_dimension"] == 1 and aut["invariant_subspace_dimension_under_inversion"] == 0 and all(not row["strict_provider_obligation_met"] for row in table) and all(boundary.values())
    payload = {
        "status": "P2708_Z12_BOUNDARY_1_COCYCLE_OBSTRUCTION_NO_STRICT_SELECTOR" if no_unlock else "P2708_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "cohomology_certificate": cohomology,
        "automorphism_action_certificate": aut,
        "premise_vs_strict_table": table,
        "source_boundary": boundary,
        "decision": {
            "boundary_1_cocycle_exports_nonpremise_selector": False,
            "current_unlock_reading": "The Z12 boundary 1-cocycle is a real new typed object with H1 dimension 1, but it exports only an orientation line ±omega.  Inversion in Aut(Z12) sends omega to -omega, so the Aut-invariant nonzero orientation subspace is zero-dimensional.  Therefore this object can support a premise orientation but cannot discharge QW-2191 or export a non-premise strict selector on current artifacts.",
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "next_honest_step": "A further admissible move must provide a strict source for the missing sign of the boundary 1-cocycle, or a different genuinely new typed object/provider outside the closed lanes.  Without that, preserve the P2697-P2708 no-new-live-frontier/no-strict-selector certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2708/S1658 Z12 boundary 1-cocycle selector-source obstruction", "## P2708/S1658 Z12 boundary 1-cocycle selector-source obstruction\n\n`P2708/S1658` tests a genuinely new typed candidate after P2707: the oriented boundary 1-cocycle line of the Z12 cycle.  The finite chain computation gives `rank(d0)=11` and `dim H1=1`, so a premise sign can orient the circle.  However Aut(Z12) inversion sends the primitive circulation `omega` to `-omega`; the nonzero Aut-invariant orientation subspace is therefore empty.  This candidate does not discharge `QW-2191`, export a non-premise selector provider, upgrade pair12 strict-core, promote `L_total`, or imply ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2708/S1658 boundary cocycle Ltotal guard", "## P2708/S1658 boundary cocycle Ltotal guard\n\n`P2708/S1658` is a finite Z12 cohomology/orientation obstruction.  It identifies a premise-orientable boundary 1-cocycle line but no strict source for its sign; it is not a variational source and does not promote `L_total`, selector closure, pair12 strict-core upgrade, role transfer, bridge closure, or ToE closure.\n")
    append_once(AGENTS, "Current Z12 boundary 1-cocycle selector-source obstruction guardrail (P2708/S1658, 2026-06-14)", "## Current Z12 boundary 1-cocycle selector-source obstruction guardrail (P2708/S1658, 2026-06-14)\n\n- P2708 tests a new typed candidate: the Z12 oriented boundary 1-cocycle line.  The finite computation gives `rank(d0)=11`, `dim H1=1`, and a real premise-orientable sign pair `±omega`.\n- Aut(Z12) inversion sends `omega` to `-omega`, leaving no nonzero Aut-invariant orientation; therefore the candidate does not export a non-premise strict selector or discharge `QW-2191`.\n- Do not promote this boundary-cocycle candidate to pair12 strict-core, `L_total`, role transfer, bridge closure, or ToE closure without a new strict source for the missing sign.\n")
    return payload


if __name__ == "__main__":
    main()
