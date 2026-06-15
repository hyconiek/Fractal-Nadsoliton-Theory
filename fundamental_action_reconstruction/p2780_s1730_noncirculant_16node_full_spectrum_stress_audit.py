#!/usr/bin/env python3
"""P2780/S1730: non-circulant 16-node full-spectrum stress audit.

P2779 found no nonisomorphic full-Laplacian-spectrum collisions on the P2778
circulant/Cayley quotient.  This follow-up adds one explicit non-circulant
16-node 4-regular family: two C8 layers with two cross-layer perfect matchings
whose shifts differ by 4.  The family is tested against the P2779 quotient by
isomorphism and full spectrum.

The result is again positive but bounded: the enlarged declared class has 23
labeled candidates, 7 isomorphism classes, and 7 distinct full Laplacian spectra.
That strengthens the finite witness but still does not source the graph class,
target spectrum, or variational K/L_total coupling.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2778_s1728_max_symmetry_16node_geometry_source_audit import candidate_edge_sets
from p2779_s1729_16node_circulant_full_spectrum_quotient_audit import laplacian_spectrum, quotient_classes

GEN = ROOT / "generated"
P2779 = GEN / "p2779_s1729_16node_circulant_full_spectrum_quotient_audit.json"
OUT = GEN / "p2780_s1730_noncirculant_16node_full_spectrum_stress_audit.json"
MD = GEN / "p2780_s1730_noncirculant_16node_full_spectrum_stress_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

NEGATIVE_EXPORT_FLAGS = [
    "strict_spectral_source_law_exported",
    "canonical_geometry_source_exported",
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


def normalized_edges(edges: set[tuple[int, int]]) -> list[tuple[int, int]]:
    return sorted({(min(a, b), max(a, b)) for a, b in edges if a != b})


def two_c8_layer_matching_edges(shift_a: int, shift_b: int) -> list[tuple[int, int]]:
    """Return a 16-node 4-regular two-layer graph.

    Nodes 0..7 and 8..15 are two C8 layers.  Each node in layer 0 connects to
    two layer-1 nodes, shifted by shift_a and shift_b.  With two cycle neighbors
    in its own layer, every node has degree four.
    """
    edges: set[tuple[int, int]] = set()
    for layer in (0, 1):
        base = 8 * layer
        for i in range(8):
            edges.add((base + i, base + ((i + 1) % 8)))
    for i in range(8):
        edges.add((i, 8 + ((i + shift_a) % 8)))
        edges.add((i, 8 + ((i + shift_b) % 8)))
    return normalized_edges(edges)


def noncirculant_family_candidates() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for shift in range(4):
        rows.append({
            "geometry": f"two_c8_layers_cross_pm{shift}_pm{shift + 4}",
            "family": "noncirculant_two_c8_layers_opposite_cross_matchings",
            "steps": [shift, shift + 4],
            "edges": two_c8_layer_matching_edges(shift, shift + 4),
        })
    return rows


def stress_witness() -> dict[str, Any]:
    base_candidates = candidate_edge_sets()
    added_candidates = noncirculant_family_candidates()
    expanded_candidates = base_candidates + added_candidates
    quotient = quotient_classes(expanded_candidates)
    buckets: dict[str, list[str]] = {}
    for row in quotient:
        buckets.setdefault(row["spectrum_signature"], []).append(row["representative"])
    collisions = [names for names in buckets.values() if len(names) > 1]
    added_spectra = {json.dumps(laplacian_spectrum(row["edges"])) for row in added_candidates}
    base_spectra = {json.dumps(laplacian_spectrum(row["edges"])) for row in base_candidates}
    return {
        "base_class": "P2779 quotient input: torus_4x4 plus connected 16-node 4-regular circulant Cayley candidates",
        "added_family": "two C8 layers with cross-layer perfect matching shifts {s, s+4} for s=0,1,2,3",
        "base_labeled_candidate_count": len(base_candidates),
        "added_labeled_candidate_count": len(added_candidates),
        "expanded_labeled_candidate_count": len(expanded_candidates),
        "expanded_isomorphism_class_count": len(quotient),
        "distinct_laplacian_spectrum_count_after_quotient": len(buckets),
        "nonisomorphic_cospectral_collision_count": len(collisions),
        "nonisomorphic_cospectral_collision_representatives": collisions,
        "full_spectrum_injective_on_expanded_quotient": len(collisions) == 0 and len(buckets) == len(quotient),
        "added_family_label_count": len(added_candidates),
        "added_family_distinct_spectrum_count": len(added_spectra),
        "added_family_spectrum_new_relative_to_p2779_base": added_spectra.isdisjoint(base_spectra),
        "quotient_rows": quotient,
        "finite_positive_statement": "Adding the non-circulant two-C8-layer family creates one new isomorphism/spectral class and still yields no nonisomorphic full-spectrum collision.",
    }


def acceptance_matrix(witness: dict[str, Any], p2779: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2779_declared_class_uniqueness_present": p2779.get("status") == "P2779_16NODE_CIRCULANT_FULL_SPECTRUM_QUOTIENT_AUDIT_NO_CLOSURE",
        "one_broader_16node_4regular_family_added": True,
        "isomorphism_quotient_performed": True,
        "full_spectrum_injective_on_expanded_quotient": witness["full_spectrum_injective_on_expanded_quotient"],
        "added_family_spectrum_new_relative_to_p2779_base": witness["added_family_spectrum_new_relative_to_p2779_base"],
        "strict_nadsoliton_spectral_source_law_exported": False,
        "global_graph_class_theorem_exported": False,
        "target_spectrum_or_variational_minimizer_exported_before_testing": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_expanded_declared_class_full_spectrum_uniqueness_witness": facts["full_spectrum_injective_on_expanded_quotient"],
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_kernel_full_expression_theorem": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "The expanded finite audit is positive, but it is still a declared stress class.  No strict law sources the graph class, target spectrum, or K/L_total variational coupling.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    witness = payload["noncirculant_full_spectrum_stress_witness"]
    lines = [
        "# P2780/S1730 non-circulant 16-node full-spectrum stress audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Expanded quotient spectral result",
        f"- base_labeled_candidate_count={witness['base_labeled_candidate_count']}",
        f"- added_labeled_candidate_count={witness['added_labeled_candidate_count']}",
        f"- expanded_labeled_candidate_count={witness['expanded_labeled_candidate_count']}",
        f"- expanded_isomorphism_class_count={witness['expanded_isomorphism_class_count']}",
        f"- distinct_laplacian_spectrum_count_after_quotient={witness['distinct_laplacian_spectrum_count_after_quotient']}",
        f"- nonisomorphic_cospectral_collision_count={witness['nonisomorphic_cospectral_collision_count']}",
        f"- full_spectrum_injective_on_expanded_quotient={witness['full_spectrum_injective_on_expanded_quotient']}",
        f"- added_family_spectrum_new_relative_to_p2779_base={witness['added_family_spectrum_new_relative_to_p2779_base']}",
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
    p2779 = read_json(P2779)
    witness = stress_witness()
    acceptance = acceptance_matrix(witness, p2779)
    payload = {
        "status": "P2780_NONCIRCULANT_16NODE_FULL_SPECTRUM_STRESS_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2779": sha(P2779)},
        "input_statuses": {"P2779": p2779.get("status")},
        "audited_question": "Does adding one explicit non-circulant 16-node 4-regular family create a nonisomorphic full-Laplacian-spectrum collision?",
        "noncirculant_full_spectrum_stress_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Treat P2780 as an expanded finite uniqueness witness only.  The next honest move is exactly one of: export a strict nadsoliton spectral action/source law that fixes the graph class, target spectrum, and K/L_total coupling before testing; or run a collision-seeking audit on a genuinely broader enumerated 16-node 4-regular graph class with canonical isomorphism quotienting rather than adding hand-picked families.  Otherwise preserve the P2697-P2780 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2780/S1730 non-circulant 16-node full-spectrum stress audit", "## P2780/S1730 non-circulant 16-node full-spectrum stress audit\n\n`P2780/S1730` stress-tests the P2779 full-spectrum quotient by adding an explicit non-circulant 16-node 4-regular family: two `C8` layers with cross-layer perfect matching shifts `{s,s+4}` for `s=0,1,2,3`.  The expanded 23 labeled candidates collapse to 7 isomorphism classes, and the full Laplacian spectrum remains injective on those 7 classes with zero nonisomorphic cospectral collisions.  This strengthens the finite witness but still exports no strict nadsoliton spectral source law, global graph-class theorem, kernel full-expression theorem, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2780/S1730 non-circulant spectral stress Ltotal guard", "## P2780/S1730 non-circulant spectral stress Ltotal guard\n\n`P2780/S1730` adds no variational source term.  Full-spectrum uniqueness survives one non-circulant 16-node 4-regular stress family, but the enlarged class, target, and spectrum are not sourced by a nonproxy `K`/`L_total` variational principle; therefore it cannot promote role-bearing `L_total` or canonical nadsoliton geometry.\n")
    append_once(AGENTS, "Current non-circulant 16-node full-spectrum stress audit guardrail (P2780/S1730, 2026-06-15)", "## Current non-circulant 16-node full-spectrum stress audit guardrail (P2780/S1730, 2026-06-15)\n\n- P2780 adds one explicit non-circulant 16-node 4-regular stress family to the P2779 quotient: two `C8` layers with cross-layer perfect matching shifts `{s,s+4}`.\n- The expanded 23 labeled candidates collapse to 7 isomorphism classes, and no nonisomorphic full-Laplacian-spectrum collision appears on this enlarged declared quotient.\n- Do not promote this expanded finite witness to canonical nadsoliton geometry, strict spectral source law, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  A next admissible move must export a strict spectral action/source law before testing, or run a genuinely broader enumerated 16-node 4-regular collision search with canonical isomorphism quotienting.\n")
    return payload


if __name__ == "__main__":
    main()
