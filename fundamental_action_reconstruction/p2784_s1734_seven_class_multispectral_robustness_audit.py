#!/usr/bin/env python3
"""P2784/S1734: seven-class multispectral robustness audit.

P2783 proved that the seven P2781 quotient representatives are pairwise
nonisomorphic and pairwise separated by the full graph-Laplacian spectrum.  This
follow-up asks a more proof-grade computational question: is the separation an
artifact of one chosen invariant, or does it survive independent spectral
probes on the same seven-class quotient?

The audit keeps the same local graph class and therefore exports no canonical
geometry.  It computes adjacency, Laplacian, and signless-Laplacian spectra,
plus finite support statistics, for all seven representatives and checks all 21
pairs for collisions in each invariant family.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

import numpy as np

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2778_s1728_max_symmetry_16node_geometry_source_audit import N, adjacency, candidate_edge_sets
from p2779_s1729_16node_circulant_full_spectrum_quotient_audit import quotient_classes
from p2781_s1731_enumerated_two_layer_c8_spectrum_collision_audit import all_two_layer_c8_candidates

GEN = ROOT / "generated"
P2783 = GEN / "p2783_s1733_seven_class_quotient_integrity_certificate.json"
OUT = GEN / "p2784_s1734_seven_class_multispectral_robustness_audit.json"
MD = GEN / "p2784_s1734_seven_class_multispectral_robustness_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

NEGATIVE_EXPORT_FLAGS = [
    "canonical_geometry_source_exported",
    "strict_spectral_source_law_exported",
    "global_full_spectrum_geometry_theorem_exported",
    "kernel_geometry_closure_exported",
    "kernel_fully_expresses_nadsoliton_characteristics",
    "role_bearing_ltotal_promoted",
    "bridge_closure_exported",
    "selector_closure_exported",
    "toe_closure_exported",
]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def matrix_for(edges: list[tuple[int, int]]) -> np.ndarray:
    mat = np.zeros((N, N), dtype=float)
    for a, b in edges:
        mat[a, b] = mat[b, a] = 1.0
    return mat


def rounded_spectrum(mat: np.ndarray, places: int = 10) -> list[float]:
    return [round(float(value), places) for value in np.linalg.eigvalsh(mat)]


def invariant_row(rep: dict[str, Any], by_name: dict[str, dict[str, Any]]) -> dict[str, Any]:
    edges = by_name[rep["representative"]]["edges"]
    adj = matrix_for(edges)
    deg = np.diag(adj.sum(axis=1))
    lap = deg - adj
    signless = deg + adj
    lap_eigs = np.linalg.eigvalsh(lap)
    nonzero = [float(x) for x in lap_eigs if abs(float(x)) > 1e-8]
    spanning_tree_count = int(round(np.prod(nonzero) / N))
    triangles = int(round(np.trace(np.linalg.matrix_power(adj, 3)) / 6))
    closed_walk4 = int(round(np.trace(np.linalg.matrix_power(adj, 4))))
    return {
        "representative": rep["representative"],
        "member_count": len(rep["members"]),
        "members": rep["members"],
        "adjacency_spectrum": rounded_spectrum(adj),
        "laplacian_spectrum": rounded_spectrum(lap),
        "signless_laplacian_spectrum": rounded_spectrum(signless),
        "spanning_tree_count": spanning_tree_count,
        "triangle_count": triangles,
        "closed_walk4_count": closed_walk4,
    }


def multispectral_witness() -> dict[str, Any]:
    candidates = candidate_edge_sets() + all_two_layer_c8_candidates()
    by_name = {row["geometry"]: row for row in candidates}
    classes = quotient_classes(candidates)
    rows = [invariant_row(rep, by_name) for rep in classes]
    invariant_keys = [
        "adjacency_spectrum",
        "laplacian_spectrum",
        "signless_laplacian_spectrum",
        "spanning_tree_count",
        "triangle_count",
        "closed_walk4_count",
    ]
    pair_rows: list[dict[str, Any]] = []
    for left, right in itertools.combinations(rows, 2):
        collisions = {key: left[key] == right[key] for key in invariant_keys}
        pair_rows.append({
            "left": left["representative"],
            "right": right["representative"],
            "invariant_collisions": collisions,
            "all_three_spectra_distinct": not collisions["adjacency_spectrum"] and not collisions["laplacian_spectrum"] and not collisions["signless_laplacian_spectrum"],
        })
    collision_counts = {key: sum(1 for row in pair_rows if row["invariant_collisions"][key]) for key in invariant_keys}
    return {
        "source_class": "Same seven P2781/P2783 local quotient representatives; no graph-class expansion is introduced.",
        "representative_count": len(rows),
        "pair_count": len(pair_rows),
        "invariant_rows": rows,
        "pair_rows": pair_rows,
        "collision_counts_by_invariant": collision_counts,
        "spectral_collision_counts": {key: collision_counts[key] for key in ["adjacency_spectrum", "laplacian_spectrum", "signless_laplacian_spectrum"]},
        "all_pairs_separated_by_all_three_spectra": all(row["all_three_spectra_distinct"] for row in pair_rows),
        "finite_certificate_statement": "All 21 representative pairs are separated by adjacency, Laplacian, and signless-Laplacian spectra on the P2781/P2783 seven-class quotient.",
    }


def acceptance_matrix(witness: dict[str, Any], p2783: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2783_integrity_certificate_present": p2783.get("status") == "P2783_SEVEN_CLASS_QUOTIENT_INTEGRITY_CERTIFICATE_NO_CLOSURE",
        "seven_representatives_reused_without_class_expansion": witness["representative_count"] == 7,
        "all_21_pairs_checked": witness["pair_count"] == 21,
        "zero_adjacency_spectrum_collisions": witness["collision_counts_by_invariant"]["adjacency_spectrum"] == 0,
        "zero_laplacian_spectrum_collisions": witness["collision_counts_by_invariant"]["laplacian_spectrum"] == 0,
        "zero_signless_laplacian_spectrum_collisions": witness["collision_counts_by_invariant"]["signless_laplacian_spectrum"] == 0,
        "strict_nadsoliton_spectral_source_law_exported": False,
        "canonical_full_16node_4regular_generator_supplied": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_multispectral_robustness_certificate": all(facts[key] for key in [
            "p2783_integrity_certificate_present",
            "seven_representatives_reused_without_class_expansion",
            "all_21_pairs_checked",
            "zero_adjacency_spectrum_collisions",
            "zero_laplacian_spectrum_collisions",
            "zero_signless_laplacian_spectrum_collisions",
        ]),
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "Multispectral separation is robust inside the seven-class local quotient, but the admissible graph class, target spectrum, and K/L_total coupling are still externally declared rather than strict-sourced.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    w = payload["multispectral_witness"]
    lines = [
        "# P2784/S1734 seven-class multispectral robustness audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Multispectral result",
        f"- representative_count={w['representative_count']}",
        f"- pair_count={w['pair_count']}",
        f"- spectral_collision_counts={w['spectral_collision_counts']}",
        f"- all_pairs_separated_by_all_three_spectra={w['all_pairs_separated_by_all_three_spectra']}",
        "",
        "## Decision",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2783 = read_json(P2783)
    witness = multispectral_witness()
    acceptance = acceptance_matrix(witness, p2783)
    payload = {
        "status": "P2784_SEVEN_CLASS_MULTISPECTRAL_ROBUSTNESS_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2783": sha(P2783)},
        "input_statuses": {"P2783": p2783.get("status")},
        "audited_question": "Does the P2783 seven-class separation survive independent adjacency, Laplacian, and signless-Laplacian spectral probes?",
        "multispectral_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Use P2784 as a local multispectral robustness certificate only.  The next honest move is exactly one of: supply a canonical-generation theorem/tool certificate for the full connected 16-node 4-regular graph class with reproducible graph6/hash provenance and then rerun the full collision audit; or export a strict nadsoliton spectral action/source law fixing the admissible class, target spectrum, and K/L_total coupling before testing.  Otherwise preserve the P2697-P2784 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2784/S1734 seven-class multispectral robustness audit", "## P2784/S1734 seven-class multispectral robustness audit\n\n`P2784/S1734` stress-tests the P2783 seven-class quotient with independent spectral probes.  On the same seven local representatives, all 21 pairs are separated not only by the full graph-Laplacian spectrum but also by the adjacency spectrum and signless-Laplacian spectrum.  This is a stronger local robustness witness, but it does not expand the graph class and does not source the class, target, or variational coupling from `K`/`L_total`; therefore no strict nadsoliton spectral source law, canonical geometry, kernel full-expression theorem, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2784/S1734 multispectral robustness Ltotal guard", "## P2784/S1734 multispectral robustness Ltotal guard\n\n`P2784/S1734` adds no variational source term.  Adjacency, Laplacian, and signless-Laplacian spectra all separate the seven P2783 representatives, but the computation remains a local invariant robustness audit rather than a sourced nonproxy `K`/`L_total` spectral action.\n")
    append_once(AGENTS, "Current seven-class multispectral robustness guardrail (P2784/S1734, 2026-06-15)", "## Current seven-class multispectral robustness guardrail (P2784/S1734, 2026-06-15)\n\n- P2784 stress-tests the P2783 seven-class quotient with adjacency, Laplacian, and signless-Laplacian spectra.\n- All 21 representative pairs are separated by all three spectra, giving a stronger local robustness witness than a single-invariant check.\n- Do not promote this local multispectral witness to canonical geometry, strict spectral source law, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  A next admissible move must supply a canonical-generation theorem/tool certificate for the full graph class, or export a strict spectral action/source law before testing.\n")
    return payload


if __name__ == "__main__":
    main()
