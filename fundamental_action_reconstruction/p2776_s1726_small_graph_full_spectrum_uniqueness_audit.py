#!/usr/bin/env python3
"""P2776/S1726: small-graph full-spectrum uniqueness audit.

P2775 showed that the full graph-Laplacian spectrum distinguishes the exact
P2774 pair, but only as a pair-local discriminator.  This bounded follow-up
runs a finite graph-class audit: all connected simple unlabeled geometries on
4 and 5 vertices, obtained by exhaustive labeled enumeration and canonical
permutation quotienting.  The full Laplacian spectrum is injective on this tiny
class, but this is still not a canonical nadsoliton geometry theorem: the class
is too small, not 16-point/strict-sourced, and has no variational coupling to
K or L_total.
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

GEN = ROOT / "generated"
P2775 = GEN / "p2775_s1725_full_laplacian_spectrum_pair_discriminator.json"
OUT = GEN / "p2776_s1726_small_graph_full_spectrum_uniqueness_audit.json"
MD = GEN / "p2776_s1726_small_graph_full_spectrum_uniqueness_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

NEGATIVE_EXPORT_FLAGS = [
    "canonical_geometry_source_exported",
    "full_spectrum_global_geometry_theorem_exported",
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


def edge_list(n: int) -> list[tuple[int, int]]:
    return list(itertools.combinations(range(n), 2))


def adjacency_from_mask(n: int, mask: int) -> list[list[int]]:
    adj = [[0] * n for _ in range(n)]
    for bit, (a, b) in enumerate(edge_list(n)):
        if (mask >> bit) & 1:
            adj[a][b] = adj[b][a] = 1
    return adj


def is_connected(adj: list[list[int]]) -> bool:
    seen = {0}
    stack = [0]
    while stack:
        node = stack.pop()
        for nxt, flag in enumerate(adj[node]):
            if flag and nxt not in seen:
                seen.add(nxt)
                stack.append(nxt)
    return len(seen) == len(adj)


def canonical_code(adj: list[list[int]]) -> str:
    """Brute-force unlabeled graph code, safe for n<=5 in this audit."""
    n = len(adj)
    best: str | None = None
    for perm in itertools.permutations(range(n)):
        code = "".join("1" if adj[perm[i]][perm[j]] else "0" for i in range(n) for j in range(i + 1, n))
        if best is None or code < best:
            best = code
    if best is None:
        raise ValueError("empty permutation set")
    return best


def laplacian_spectrum(adj: list[list[int]], places: int = 10) -> list[float]:
    degrees = [sum(row) for row in adj]
    mat = np.diag(degrees) - np.array(adj, dtype=float)
    return [round(float(value), places) for value in np.linalg.eigvalsh(mat)]


def degree_sequence(adj: list[list[int]]) -> list[int]:
    return sorted(sum(row) for row in adj)


def audit_order(n: int) -> dict[str, Any]:
    total_labeled = 1 << len(edge_list(n))
    connected_labeled = 0
    class_by_code: dict[str, dict[str, Any]] = {}
    for mask in range(total_labeled):
        adj = adjacency_from_mask(n, mask)
        if not is_connected(adj):
            continue
        connected_labeled += 1
        code = canonical_code(adj)
        if code not in class_by_code:
            spectrum = laplacian_spectrum(adj)
            class_by_code[code] = {
                "canonical_code": code,
                "edge_count": sum(sum(row) for row in adj) // 2,
                "degree_sequence": degree_sequence(adj),
                "laplacian_spectrum_rounded": spectrum,
                "spectrum_key": json.dumps(spectrum),
            }
    spectrum_to_codes: dict[str, list[str]] = {}
    for code, row in class_by_code.items():
        spectrum_to_codes.setdefault(row["spectrum_key"], []).append(code)
    collisions = [codes for codes in spectrum_to_codes.values() if len(codes) > 1]
    sample_rows = sorted(class_by_code.values(), key=lambda row: (row["edge_count"], row["degree_sequence"], row["canonical_code"]))[:8]
    return {
        "vertex_count": n,
        "total_labeled_graphs": total_labeled,
        "connected_labeled_graphs": connected_labeled,
        "connected_unlabeled_classes": len(class_by_code),
        "distinct_laplacian_spectra": len(spectrum_to_codes),
        "laplacian_cospectral_nonisomorphic_collision_count": len(collisions),
        "full_spectrum_injective_on_unlabeled_connected_classes": len(collisions) == 0,
        "sample_unlabeled_rows": sample_rows,
    }


def uniqueness_witness() -> dict[str, Any]:
    rows = [audit_order(4), audit_order(5)]
    total_classes = sum(row["connected_unlabeled_classes"] for row in rows)
    total_spectra = sum(row["distinct_laplacian_spectra"] for row in rows)
    total_collisions = sum(row["laplacian_cospectral_nonisomorphic_collision_count"] for row in rows)
    return {
        "audited_class": "all connected simple unlabeled graphs on 4 and 5 vertices, quotiented by brute-force canonical permutation code",
        "orders": rows,
        "total_connected_unlabeled_classes": total_classes,
        "total_distinct_laplacian_spectra_by_order": total_spectra,
        "total_cospectral_nonisomorphic_collision_count": total_collisions,
        "full_spectrum_injective_on_this_finite_class": total_collisions == 0,
        "finite_positive_statement": "No Laplacian-cospectral nonisomorphic collisions occur among connected simple graphs on 4 or 5 vertices in this exhaustive audit.",
    }


def acceptance_matrix(witness: dict[str, Any], p2775: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2775_pair_local_discriminator_present": p2775.get("status") == "P2775_FULL_LAPLACIAN_SPECTRUM_PAIR_DISCRIMINATOR_NO_CLOSURE",
        "finite_graph_class_exhaustively_quotiented": True,
        "full_spectrum_injective_on_this_finite_class": witness["full_spectrum_injective_on_this_finite_class"],
        "strict_nadsoliton_spectral_source_law_exported": False,
        "sixteen_point_or_strict_graph_class_audited": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_small_class_spectral_uniqueness_witness": facts["full_spectrum_injective_on_this_finite_class"],
        "accepted_as_global_full_spectrum_geometry_theorem": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_kernel_full_expression_theorem": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "The full Laplacian spectrum is injective only on the audited tiny connected graph class n=4,5.  Current artifacts still lack a strict spectral source law, a 16-point/strict graph-class uniqueness audit, and a variational coupling to K/L_total.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    witness = payload["small_graph_full_spectrum_uniqueness_witness"]
    lines = [
        "# P2776/S1726 small-graph full-spectrum uniqueness audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite class result",
        f"- total_connected_unlabeled_classes={witness['total_connected_unlabeled_classes']}",
        f"- total_cospectral_nonisomorphic_collision_count={witness['total_cospectral_nonisomorphic_collision_count']}",
        f"- full_spectrum_injective_on_this_finite_class={witness['full_spectrum_injective_on_this_finite_class']}",
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
    p2775 = read_json(P2775)
    witness = uniqueness_witness()
    acceptance = acceptance_matrix(witness, p2775)
    payload = {
        "status": "P2776_SMALL_GRAPH_FULL_SPECTRUM_UNIQUENESS_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2775": sha(P2775)},
        "input_statuses": {"P2775": p2775.get("status")},
        "audited_question": "Is the full Laplacian spectrum injective beyond the P2775 pair on a finite connected graph class?",
        "small_graph_full_spectrum_uniqueness_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Treat P2776 as a small-class positive uniqueness witness only.  The next honest move is exactly one of: extend the cospectral-degeneracy audit to a declared 16-point/regular/strict graph class with an optimized canonical quotient, or supply a strict nadsoliton spectral action/source law whose target spectrum and admissible graph class are specified before testing.  Without that, preserve the P2697-P2776 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2776/S1726 small-graph full-spectrum uniqueness audit", "## P2776/S1726 small-graph full-spectrum uniqueness audit\n\n`P2776/S1726` extends P2775 from a pair-local discriminator to an exhaustive tiny-class audit: all connected simple unlabeled graphs on 4 and 5 vertices are generated by labeled enumeration and brute-force canonical permutation quotienting.  The full Laplacian spectrum is injective on this finite class, with zero cospectral nonisomorphic collisions.  This is not canonical nadsoliton geometry: no strict spectral source law, 16-point/strict graph-class uniqueness theorem, kernel full-expression theorem, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2776/S1726 small spectral uniqueness Ltotal guard", "## P2776/S1726 small spectral uniqueness Ltotal guard\n\n`P2776/S1726` adds no variational source term.  It confirms full-Laplacian-spectrum uniqueness only on the tiny connected graph class with 4 and 5 vertices, without a sourced spectral action or nonproxy coupling to `K`/`L_total`; therefore it cannot promote role-bearing `L_total` or canonical nadsoliton geometry.\n")
    append_once(AGENTS, "Current small-graph full-spectrum uniqueness audit guardrail (P2776/S1726, 2026-06-15)", "## Current small-graph full-spectrum uniqueness audit guardrail (P2776/S1726, 2026-06-15)\n\n- P2776 extends P2775 from the exact P2774 pair to an exhaustive finite quotient of all connected simple unlabeled graphs on 4 and 5 vertices.\n- The full Laplacian spectrum is injective on this tiny class, with zero cospectral nonisomorphic collisions, but this does not supply a strict nadsoliton spectral source law, 16-point/strict graph-class uniqueness, or `K`/`L_total` variational coupling.\n- Do not promote this small-class spectral uniqueness witness to canonical nadsoliton geometry, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure.  A next admissible move must either extend the cospectral audit to a declared 16-point/regular/strict graph class or provide a sourced spectral action law before testing.\n")
    return payload


if __name__ == "__main__":
    main()
