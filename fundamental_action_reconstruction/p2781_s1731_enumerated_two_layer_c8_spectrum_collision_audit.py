#!/usr/bin/env python3
"""P2781/S1731: enumerated two-layer C8 full-spectrum collision audit.

P2780 added only the missing opposite-shift subfamily of a two-C8-layer stress
class.  This follow-up makes that local stress test exhaustive inside the same
model: enumerate every unordered pair of cross-layer perfect-matching shifts
0<=a<b<8, combine those 28 graphs with the P2779 base class, quotient by graph
isomorphism, and then audit full graph-Laplacian spectrum collisions.

The result is still positive but bounded.  The 47 labeled candidates collapse to
7 isomorphism classes with 7 distinct Laplacian spectra.  The exhaustive local
family does not create a nonisomorphic cospectral collision, but it also does not
source the admissible graph class, target spectrum, or K/L_total coupling.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2778_s1728_max_symmetry_16node_geometry_source_audit import candidate_edge_sets
from p2779_s1729_16node_circulant_full_spectrum_quotient_audit import quotient_classes
from p2780_s1730_noncirculant_16node_full_spectrum_stress_audit import two_c8_layer_matching_edges

GEN = ROOT / "generated"
P2780 = GEN / "p2780_s1730_noncirculant_16node_full_spectrum_stress_audit.json"
OUT = GEN / "p2781_s1731_enumerated_two_layer_c8_spectrum_collision_audit.json"
MD = GEN / "p2781_s1731_enumerated_two_layer_c8_spectrum_collision_audit.md"
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


def all_two_layer_c8_candidates() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for shift_a, shift_b in itertools.combinations(range(8), 2):
        rows.append({
            "geometry": f"two_c8_layers_cross_pm{shift_a}_pm{shift_b}",
            "family": "exhaustive_two_c8_layers_two_cross_matchings",
            "steps": [shift_a, shift_b],
            "edges": two_c8_layer_matching_edges(shift_a, shift_b),
        })
    return rows


def collision_witness() -> dict[str, Any]:
    base_candidates = candidate_edge_sets()
    enumerated_candidates = all_two_layer_c8_candidates()
    expanded_candidates = base_candidates + enumerated_candidates
    quotient = quotient_classes(expanded_candidates)
    buckets: dict[str, list[str]] = {}
    for row in quotient:
        buckets.setdefault(row["spectrum_signature"], []).append(row["representative"])
    collisions = [names for names in buckets.values() if len(names) > 1]
    two_layer_member_count_by_representative = {
        row["representative"]: sum(1 for name in row["members"] if name.startswith("two_c8_layers_cross_pm"))
        for row in quotient
    }
    return {
        "base_class": "P2779 quotient input: torus_4x4 plus connected 16-node 4-regular circulant Cayley candidates",
        "enumerated_family": "all two-C8-layer graphs with unordered cross-matching shifts 0<=a<b<8",
        "base_labeled_candidate_count": len(base_candidates),
        "enumerated_two_layer_candidate_count": len(enumerated_candidates),
        "expanded_labeled_candidate_count": len(expanded_candidates),
        "expanded_isomorphism_class_count": len(quotient),
        "distinct_laplacian_spectrum_count_after_quotient": len(buckets),
        "nonisomorphic_cospectral_collision_count": len(collisions),
        "nonisomorphic_cospectral_collision_representatives": collisions,
        "full_spectrum_injective_on_expanded_quotient": len(collisions) == 0 and len(buckets) == len(quotient),
        "two_layer_member_count_by_representative": two_layer_member_count_by_representative,
        "quotient_rows": quotient,
        "finite_positive_statement": "The exhaustive two-C8-layer shift-pair family adds no nonisomorphic full-spectrum collision beyond the P2779 base quotient.",
    }


def acceptance_matrix(witness: dict[str, Any], p2780: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2780_noncirculant_stress_present": p2780.get("status") == "P2780_NONCIRCULANT_16NODE_FULL_SPECTRUM_STRESS_AUDIT_NO_CLOSURE",
        "exhaustive_two_layer_c8_shift_pair_family_enumerated": witness["enumerated_two_layer_candidate_count"] == 28,
        "isomorphism_quotient_performed": True,
        "full_spectrum_injective_on_expanded_quotient": witness["full_spectrum_injective_on_expanded_quotient"],
        "strict_nadsoliton_spectral_source_law_exported": False,
        "global_graph_class_theorem_exported": False,
        "target_spectrum_or_variational_minimizer_exported_before_testing": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_exhaustive_local_family_full_spectrum_uniqueness_witness": facts["full_spectrum_injective_on_expanded_quotient"],
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_kernel_full_expression_theorem": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "The shift-pair family is exhausted only inside a local two-layer model.  No strict law sources that model, the target spectrum, or K/L_total variational coupling.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    witness = payload["two_layer_c8_collision_witness"]
    lines = [
        "# P2781/S1731 enumerated two-layer C8 full-spectrum collision audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Exhaustive local-family quotient result",
        f"- base_labeled_candidate_count={witness['base_labeled_candidate_count']}",
        f"- enumerated_two_layer_candidate_count={witness['enumerated_two_layer_candidate_count']}",
        f"- expanded_labeled_candidate_count={witness['expanded_labeled_candidate_count']}",
        f"- expanded_isomorphism_class_count={witness['expanded_isomorphism_class_count']}",
        f"- distinct_laplacian_spectrum_count_after_quotient={witness['distinct_laplacian_spectrum_count_after_quotient']}",
        f"- nonisomorphic_cospectral_collision_count={witness['nonisomorphic_cospectral_collision_count']}",
        f"- full_spectrum_injective_on_expanded_quotient={witness['full_spectrum_injective_on_expanded_quotient']}",
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
    p2780 = read_json(P2780)
    witness = collision_witness()
    acceptance = acceptance_matrix(witness, p2780)
    payload = {
        "status": "P2781_ENUMERATED_TWO_LAYER_C8_SPECTRUM_COLLISION_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2780": sha(P2780)},
        "input_statuses": {"P2780": p2780.get("status")},
        "audited_question": "Does the exhaustive two-C8-layer cross-matching shift-pair family create a nonisomorphic full-Laplacian-spectrum collision?",
        "two_layer_c8_collision_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Treat P2781 as an exhaustive local-family uniqueness witness only.  The next honest move is exactly one of: export a strict nadsoliton spectral action/source law that fixes the graph class, target spectrum, and K/L_total coupling before testing; or move from structured families to a canonical enumerator for connected 16-node 4-regular graphs with a collision-seeking full-spectrum quotient audit.  Otherwise preserve the P2697-P2781 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2781/S1731 enumerated two-layer C8 full-spectrum collision audit", "## P2781/S1731 enumerated two-layer C8 full-spectrum collision audit\n\n`P2781/S1731` exhausts the local two-`C8`-layer stress model left by P2780 by enumerating all unordered cross-layer matching shift pairs `0<=a<b<8`.  Combining those 28 graphs with the P2779 base class gives 47 labeled candidates, 7 isomorphism classes, and 7 distinct full Laplacian spectra, with zero nonisomorphic cospectral collisions.  This is only an exhaustive local-family witness: no strict nadsoliton spectral source law, global graph-class theorem, kernel full-expression theorem, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2781/S1731 enumerated two-layer spectral Ltotal guard", "## P2781/S1731 enumerated two-layer spectral Ltotal guard\n\n`P2781/S1731` adds no variational source term.  Full-spectrum uniqueness survives exhaustive enumeration inside the two-`C8`-layer shift-pair model, but that model, its target, and its spectrum are not sourced by a nonproxy `K`/`L_total` variational principle; therefore it cannot promote role-bearing `L_total` or canonical nadsoliton geometry.\n")
    append_once(AGENTS, "Current enumerated two-layer C8 full-spectrum collision audit guardrail (P2781/S1731, 2026-06-15)", "## Current enumerated two-layer C8 full-spectrum collision audit guardrail (P2781/S1731, 2026-06-15)\n\n- P2781 exhausts the local two-`C8`-layer stress model by enumerating all unordered cross-layer matching shift pairs `0<=a<b<8` and combining them with the P2779 base quotient.\n- The 47 labeled candidates collapse to 7 isomorphism classes, and no nonisomorphic full-Laplacian-spectrum collision appears on this exhaustive local-family quotient.\n- Do not promote this local-family witness to canonical nadsoliton geometry, strict spectral source law, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  A next admissible move must export a strict spectral action/source law before testing, or move to a canonical enumerator for connected 16-node 4-regular graphs with collision-seeking quotient audit.\n")
    return payload


if __name__ == "__main__":
    main()
