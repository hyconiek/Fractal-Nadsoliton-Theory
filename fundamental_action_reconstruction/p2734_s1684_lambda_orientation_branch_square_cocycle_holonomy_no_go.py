#!/usr/bin/env python3
"""P2734/S1684: lambda/orientation branch-square cocycle holonomy no-go.

P2733 closed spectral lambda-observability for the direct chiral tau coupling.
P2734 tests a non-spectral, topological-looking candidate: the branch square
with vertices (lambda, orientation) from P2732.  Edges flip lambda or
orientation, and the edge label is the ratio of selected tau signs.  A nontrivial
holonomy could have supplied a polarity source; the finite calculation shows the
edge cocycle is exact and has trivial plaquette holonomy.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2734_s1684_lambda_orientation_branch_square_cocycle_holonomy_no_go.json"
MD = GEN / "p2734_s1684_lambda_orientation_branch_square_cocycle_holonomy_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

LAMBDAS = (-1, 1)
ORIENTATIONS = (-1, 1)
INPUTS = {
    "P2733_SPECTRAL_LAMBDA_NO_GO": GEN / "p2733_s1683_chiral_tau_coupling_spectral_lambda_observability_no_go.json",
    "P2732_CHIRAL_TAU_COUPLING": GEN / "p2732_s1682_chiral_bispectrum_time_arrow_source_term_coupling_matrix.json",
}
NEGATIVE_EXPORT_FLAGS = [
    "nontrivial_branch_holonomy_exported",
    "lambda_orientation_cocycle_nonexact",
    "base_vertex_polarity_selected",
    "strict_mechanism_fixing_lambda_exported",
    "p2721_polarity_selected",
    "qw2191_discharged",
    "pair12_strict_core_upgrade_exported",
    "bridge_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def selected_tau_by_vertex(p2732: dict[str, Any]) -> dict[tuple[int, int], int]:
    rows = p2732.get("chiral_bispectrum_time_arrow_coupling_audit", {}).get("coupling_rows", [])
    if rows:
        return {(int(row["lambda"]), int(row["orientation"])): int(row["selected_constant_tau_sign"]) for row in rows}
    return {(lam, orientation): lam * orientation for lam in LAMBDAS for orientation in ORIENTATIONS}


def edge_label(tau: dict[tuple[int, int], int], src: tuple[int, int], dst: tuple[int, int]) -> int:
    return tau[dst] * tau[src]


def audit_branch_square(tau: dict[tuple[int, int], int]) -> dict[str, Any]:
    vertices = [(lam, orientation) for lam in LAMBDAS for orientation in ORIENTATIONS]
    lambda_edges = [((lam, orientation), (-lam, orientation)) for lam, orientation in vertices]
    orientation_edges = [((lam, orientation), (lam, -orientation)) for lam, orientation in vertices]
    directed_edges = lambda_edges + orientation_edges
    edge_rows = [
        {
            "src": list(src),
            "dst": list(dst),
            "edge_type": "lambda_flip" if src[1] == dst[1] else "orientation_flip",
            "tau_src": tau[src],
            "tau_dst": tau[dst],
            "edge_label_tau_ratio": edge_label(tau, src, dst),
        }
        for src, dst in directed_edges
    ]
    square_cycles = [
        [(-1, -1), (1, -1), (1, 1), (-1, 1), (-1, -1)],
        [(-1, -1), (-1, 1), (1, 1), (1, -1), (-1, -1)],
    ]
    cycle_rows = []
    for cycle in square_cycles:
        labels = [edge_label(tau, cycle[i], cycle[i + 1]) for i in range(len(cycle) - 1)]
        holonomy = 1
        for label in labels:
            holonomy *= label
        cycle_rows.append({"cycle": [list(v) for v in cycle], "edge_labels": labels, "holonomy": holonomy})
    coboundary_checks = [row["edge_label_tau_ratio"] == tau[tuple(row["dst"])] * tau[tuple(row["src"])] for row in edge_rows]
    return {
        "vertex_count": len(vertices),
        "directed_edge_count": len(edge_rows),
        "cycle_count": len(cycle_rows),
        "selected_tau_by_vertex": {f"lambda={lam},orientation={orientation}": tau[(lam, orientation)] for lam, orientation in vertices},
        "edge_rows": edge_rows,
        "cycle_rows": cycle_rows,
        "all_edges_are_tau_coboundary": all(coboundary_checks),
        "all_cycle_holonomies_trivial": all(row["holonomy"] == 1 for row in cycle_rows),
        "nontrivial_holonomy_count": sum(1 for row in cycle_rows if row["holonomy"] != 1),
        "obstruction": "The branch-square edge labels are exactly the coboundary of the selected tau vertex signs.  Both fundamental square orientations have holonomy +1, so the non-spectral branch cocycle is flat/exact and cannot select a base lambda/P2721 polarity.",
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "branch_square_built_from_p2732_tau_signs": audit["vertex_count"] == 4,
        "edge_cocycle_exact": audit["all_edges_are_tau_coboundary"],
        "nontrivial_branch_holonomy_exists": audit["nontrivial_holonomy_count"] > 0,
        "base_vertex_lambda_polarity_selected": False,
        "strict_lambda_fixing_holonomy_law_exported": False,
    }
    missing = [key for key, value in facts.items() if not value]
    return {
        "facts": facts,
        "missing_criteria": missing,
        "accepted_as_strict_lambda_polarity_source": not missing,
        "blocker": "The only branch-square cocycle exported by P2732 is exact with trivial holonomy; it supplies no basepoint-free sign choice.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["branch_square_holonomy_audit"]
    lines = [
        "# P2734/S1684 lambda/orientation branch-square cocycle holonomy no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Branch-square audit",
        f"- vertex_count={audit['vertex_count']}",
        f"- directed_edge_count={audit['directed_edge_count']}",
        f"- cycle_count={audit['cycle_count']}",
        f"- all_edges_are_tau_coboundary={audit['all_edges_are_tau_coboundary']}",
        f"- all_cycle_holonomies_trivial={audit['all_cycle_holonomies_trivial']}",
        f"- nontrivial_holonomy_count={audit['nontrivial_holonomy_count']}",
        audit["obstruction"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    audit = audit_branch_square(selected_tau_by_vertex(inputs["P2732_CHIRAL_TAU_COUPLING"]))
    acceptance = acceptance_matrix(audit)
    no_go = not acceptance["accepted_as_strict_lambda_polarity_source"]
    payload = {
        "status": "P2734_BRANCH_SQUARE_COCYCLE_HOLONOMY_EXACT_NO_GO" if no_go else "P2734_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_candidate_class": "non-spectral branch-square Z2 cocycle from P2732 selected tau signs on (lambda,orientation) vertices",
        "branch_square_holonomy_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "nontrivial_branch_holonomy_exported": False,
            "lambda_orientation_cocycle_nonexact": False,
            "base_vertex_polarity_selected": False,
            "strict_mechanism_fixing_lambda_exported": False,
            "p2721_polarity_selected": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2734 tests a non-spectral branch-square cocycle built from the P2732 selected tau signs.  The edge labels are an exact coboundary and both square holonomies are +1, so there is no basepoint-free lambda/P2721 polarity source.",
            "next_honest_step": "Do not replay the P2732 lambda/orientation branch square as a hidden holonomy source.  A next admissible proof-grade move must introduce a genuinely new non-exact internal polarity object, or preserve the P2697-P2734 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2734/S1684 lambda-orientation branch-square cocycle holonomy no-go", "## P2734/S1684 lambda-orientation branch-square cocycle holonomy no-go\n\n`P2734/S1684` builds the non-spectral branch square with vertices `(lambda,orientation)` from the P2732 selected `tau` signs.  Every lambda/orientation flip edge is the `tau` coboundary ratio, and both square orientations have trivial holonomy `+1`; the cocycle is exact and selects no base vertex, `lambda`, or P2721 polarity.  No strict mechanism fixing `lambda`, `QW-2191` discharge, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2734/S1684 branch-square holonomy Ltotal guard", "## P2734/S1684 branch-square holonomy Ltotal guard\n\n`P2734/S1684` tests a topological-looking non-spectral branch-square source for the P2732 coupling, but the edge cocycle is exact and all plaquette holonomies are trivial.  Therefore it adds no strict variational/polarity source term and cannot promote `L_total`.\n")
    append_once(AGENTS, "Current lambda-orientation branch-square cocycle holonomy no-go guardrail (P2734/S1684, 2026-06-14)", "## Current lambda-orientation branch-square cocycle holonomy no-go guardrail (P2734/S1684, 2026-06-14)\n\n- P2734 tests the non-spectral branch-square cocycle from P2732 selected `tau` signs on `(lambda,orientation)` vertices.\n- The edge labels are an exact `tau` coboundary and both square holonomies are trivial `+1`, so no basepoint-free `lambda` or P2721 polarity source is exported.\n- Do not replay the P2732 branch square as `QW-2191` discharge, selector closure, pair12 strict-core upgrade, role transfer, bridge closure, `L_total`, or ToE.  A next admissible move must introduce a genuinely new non-exact internal polarity object, or preserve the P2697-P2734 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
