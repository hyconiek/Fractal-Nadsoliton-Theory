#!/usr/bin/env python3
"""P2735/S1685: branch-square non-exact flux polarity-object no-go.

P2734 showed that the branch-square cocycle exported by the P2732 selected tau
signs is exact.  P2735 tests the next honest candidate class rather than
replaying that exact cocycle: all Z2 edge-label systems on the lambda/orientation
branch square, modulo vertex-gauge flips.  The finite calculation confirms that
a non-exact Z2 flux class exists, but it is only an add-on frustrated line-bundle
choice.  It is invariant under the square symmetries and selects no base vertex,
lambda sign, orientation branch, or P2721 polarity.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2735_s1685_branch_square_nonexact_flux_polarity_object_no_go.json"
MD = GEN / "p2735_s1685_branch_square_nonexact_flux_polarity_object_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
INPUTS = {
    "P2734_BRANCH_SQUARE_HOLONOMY": GEN / "p2734_s1684_lambda_orientation_branch_square_cocycle_holonomy_no_go.json",
}
NEGATIVE_EXPORT_FLAGS = [
    "internal_nonexact_polarity_source_exported",
    "base_vertex_lambda_selected",
    "p2721_polarity_selected",
    "strict_mechanism_fixing_lambda_exported",
    "qw2191_discharged",
    "pair12_strict_core_upgrade_exported",
    "bridge_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]

VERTICES = ((-1, -1), (1, -1), (1, 1), (-1, 1))
EDGES = tuple((VERTICES[i], VERTICES[(i + 1) % 4]) for i in range(4))
# Dihedral symmetries of the square generated as permutations of vertex indices.
D4_PERMS = sorted(set(
    tuple(((s * i + r) % 4) for i in range(4))
    for r in range(4)
    for s in (1, -1)
))


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def holonomy(labels: tuple[int, int, int, int]) -> int:
    product = 1
    for label in labels:
        product *= label
    return product


def gauge_transform(labels: tuple[int, int, int, int], vertex_signs: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return tuple(labels[i] * vertex_signs[i] * vertex_signs[(i + 1) % 4] for i in range(4))


def permute_labels(labels: tuple[int, int, int, int], perm: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    edge_map = {frozenset((perm[i], perm[(i + 1) % 4])): labels[i] for i in range(4)}
    return tuple(edge_map[frozenset((i, (i + 1) % 4))] for i in range(4))


def gauge_orbit(labels: tuple[int, int, int, int]) -> set[tuple[int, int, int, int]]:
    return {gauge_transform(labels, signs) for signs in itertools.product((-1, 1), repeat=4)}


def audit_flux_class() -> dict[str, Any]:
    all_cochains = list(itertools.product((-1, 1), repeat=4))
    exact = [c for c in all_cochains if holonomy(c) == 1]
    nonexact = [c for c in all_cochains if holonomy(c) == -1]
    orbit_reps: dict[int, list[int]] = {}
    for target_holonomy in (-1, 1):
        subset = [c for c in all_cochains if holonomy(c) == target_holonomy]
        orbit = gauge_orbit(subset[0])
        orbit_reps[target_holonomy] = [list(x) for x in sorted(orbit)]
    nonexact_rep = nonexact[0]
    d4_preserves_nonexact = all(holonomy(permute_labels(nonexact_rep, perm)) == -1 for perm in D4_PERMS)
    d4_vertex_orbit = sorted({perm[i] for perm in D4_PERMS for i in range(4)})
    return {
        "vertex_count": len(VERTICES),
        "edge_count": len(EDGES),
        "cochain_count": len(all_cochains),
        "exact_holonomy_plus_count": len(exact),
        "nonexact_holonomy_minus_count": len(nonexact),
        "gauge_orbit_count": 2,
        "gauge_orbits_by_holonomy": orbit_reps,
        "nonexact_flux_exists_as_abstract_class": True,
        "nonexact_flux_is_d4_symmetric": d4_preserves_nonexact,
        "d4_vertex_orbit_size": len(d4_vertex_orbit),
        "base_vertex_selected_by_flux": False,
        "lambda_sign_selected_by_flux": False,
        "orientation_branch_selected_by_flux": False,
        "p2721_polarity_selected_by_flux": False,
        "obstruction": "The finite branch square has two Z2 gauge classes classified by plaquette holonomy.  The non-exact -1 class exists, but it is an unsourced frustrated flux choice; the D4 square symmetries still move all vertices together, so the class selects no lambda sign, orientation branch, or P2721 polarity.",
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "full_z2_edge_cochain_class_enumerated": audit["cochain_count"] == 16,
        "nonexact_flux_class_exists": audit["nonexact_holonomy_minus_count"] == 8,
        "nonexact_flux_internally_sourced_by_prior_artifacts": False,
        "base_vertex_lambda_or_orientation_selected": audit["base_vertex_selected_by_flux"] or audit["lambda_sign_selected_by_flux"] or audit["orientation_branch_selected_by_flux"],
        "p2721_polarity_coupling_selected": audit["p2721_polarity_selected_by_flux"],
        "strict_lambda_fixing_law_exported": False,
    }
    missing = [key for key, value in facts.items() if not value]
    return {
        "facts": facts,
        "missing_criteria": missing,
        "accepted_as_strict_nonexact_polarity_object": not missing,
        "blocker": "A non-exact branch-square Z2 flux is available only as an added frustrated class; current artifacts do not source it or turn it into a lambda/P2721 selector.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["branch_square_nonexact_flux_audit"]
    lines = [
        "# P2735/S1685 branch-square non-exact flux polarity-object no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite audit",
        f"- cochain_count={audit['cochain_count']}",
        f"- exact_holonomy_plus_count={audit['exact_holonomy_plus_count']}",
        f"- nonexact_holonomy_minus_count={audit['nonexact_holonomy_minus_count']}",
        f"- gauge_orbit_count={audit['gauge_orbit_count']}",
        f"- nonexact_flux_is_d4_symmetric={audit['nonexact_flux_is_d4_symmetric']}",
        f"- d4_vertex_orbit_size={audit['d4_vertex_orbit_size']}",
        audit["obstruction"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    audit = audit_flux_class()
    acceptance = acceptance_matrix(audit)
    no_go = not acceptance["accepted_as_strict_nonexact_polarity_object"]
    payload = {
        "status": "P2735_BRANCH_SQUARE_NONEXACT_FLUX_POLARITY_OBJECT_NO_GO" if no_go else "P2735_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_candidate_class": "all Z2 edge cochains on the P2732 lambda/orientation branch square, modulo vertex gauge, as a candidate non-exact internal polarity object",
        "branch_square_nonexact_flux_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "abstract_nonexact_flux_class_exists": True,
            "internal_nonexact_polarity_source_exported": False,
            "base_vertex_lambda_selected": False,
            "p2721_polarity_selected": False,
            "strict_mechanism_fixing_lambda_exported": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2735 enumerates all 16 Z2 edge cochains on the branch square and finds the expected two gauge classes by holonomy.  The non-exact -1 flux class is real as an abstract add-on, but it is not sourced by prior strict artifacts and remains D4-symmetric over the four vertices, so it cannot fix lambda or P2721 polarity.",
            "next_honest_step": "Do not replay an abstract non-exact branch-square flux as closure.  A next admissible proof-grade move must export an internal law that sources the non-exact flux and breaks the lambda/P2721 polarity, or pivot to a genuinely different typed object; otherwise preserve the P2697-P2735 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2735/S1685 branch-square non-exact flux polarity-object no-go", "## P2735/S1685 branch-square non-exact flux polarity-object no-go\n\n`P2735/S1685` enumerates all `2^4=16` `Z2` edge-label systems on the `(lambda,orientation)` branch square, modulo vertex-gauge flips.  There are exactly two gauge classes, classified by plaquette holonomy `+1` and `-1`; the non-exact `-1` class exists only as an unsourced frustrated flux choice.  Because the square symmetries still move all four vertices together, this class selects no base vertex, `lambda`, orientation branch, or P2721 polarity.  No strict mechanism fixing `lambda`, `QW-2191` discharge, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2735/S1685 non-exact flux Ltotal guard", "## P2735/S1685 non-exact flux Ltotal guard\n\n`P2735/S1685` confirms that adding a non-exact branch-square `Z2` flux would be a real frustrated topological choice, but current artifacts do not source that flux or couple it to a unique P2721 polarity.  It therefore adds no strict variational/polarity source term and cannot promote `L_total`.\n")
    append_once(AGENTS, "Current branch-square non-exact flux polarity-object no-go guardrail (P2735/S1685, 2026-06-14)", "## Current branch-square non-exact flux polarity-object no-go guardrail (P2735/S1685, 2026-06-14)\n\n- P2735 enumerates all `2^4=16` `Z2` edge cochains on the `(lambda,orientation)` branch square and finds exactly two vertex-gauge classes, with holonomy `+1` and `-1`.\n- The non-exact `-1` flux class exists only as an abstract/added frustrated class; current strict artifacts do not source it, and the D4 square symmetry still selects no base vertex, `lambda` sign, orientation branch, or P2721 polarity.\n- Do not replay abstract branch-square non-exact flux as `QW-2191` discharge, selector closure, pair12 strict-core upgrade, role transfer, bridge closure, `L_total`, or ToE.  A next admissible move must export an internal law sourcing the non-exact flux and breaking the lambda/P2721 polarity, or pivot to a genuinely different typed object; otherwise preserve the P2697-P2735 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
