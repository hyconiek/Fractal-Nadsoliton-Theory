#!/usr/bin/env python3
"""P2779/S1729: 16-node circulant full-spectrum quotient audit.

P2778 showed that maximal symmetry is not a safe geometry source on the declared
16-node 4-regular class.  This bounded follow-up returns to the spectral lane
left open by P2776/P2778 and asks a narrower computational question: after
quotienting the same finite 16-node candidate class by graph isomorphism, does
the full graph-Laplacian spectrum collide between nonisomorphic geometries?

The result is positive only for this declared class.  The quotient has six
isomorphism classes and six distinct Laplacian spectra, so no nonisomorphic
cospectral collision is found here.  This still does not export a strict
nadsoliton spectral source law, canonical geometry, or K/L_total coupling.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2778_s1728_max_symmetry_16node_geometry_source_audit import N, adjacency, candidate_edge_sets

GEN = ROOT / "generated"
P2778 = GEN / "p2778_s1728_max_symmetry_16node_geometry_source_audit.json"
OUT = GEN / "p2779_s1729_16node_circulant_full_spectrum_quotient_audit.json"
MD = GEN / "p2779_s1729_16node_circulant_full_spectrum_quotient_audit.md"
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


def laplacian_spectrum(edges: list[tuple[int, int]], places: int = 10) -> list[float]:
    mat = np.zeros((N, N), dtype=float)
    for a, b in edges:
        mat[a, b] = mat[b, a] = -1.0
        mat[a, a] += 1.0
        mat[b, b] += 1.0
    return [round(float(value), places) for value in np.linalg.eigvalsh(mat)]


def isomorphic(edges_a: list[tuple[int, int]], edges_b: list[tuple[int, int]]) -> bool:
    adj_a = [set(row) for row in adjacency(edges_a)]
    adj_b = [set(row) for row in adjacency(edges_b)]
    order = sorted(range(N), key=lambda node: (-len(adj_a[node]), node))
    mapping: dict[int, int] = {}
    used: set[int] = set()

    def compatible(source: int, target: int) -> bool:
        for mapped_source, mapped_target in mapping.items():
            if (mapped_source in adj_a[source]) != (mapped_target in adj_b[target]):
                return False
        return True

    def rec(index: int) -> bool:
        if index == N:
            return True
        source = order[index]
        for target in range(N):
            if target in used:
                continue
            if not compatible(source, target):
                continue
            mapping[source] = target
            used.add(target)
            if rec(index + 1):
                return True
            used.remove(target)
            del mapping[source]
        return False

    return rec(0)


def quotient_classes(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    classes: list[dict[str, Any]] = []
    for candidate in candidates:
        for row in classes:
            if isomorphic(candidate["edges"], row["representative_edges"]):
                row["members"].append(candidate["geometry"])
                row["member_steps"].append(candidate["steps"])
                break
        else:
            classes.append({
                "representative": candidate["geometry"],
                "representative_edges": candidate["edges"],
                "members": [candidate["geometry"]],
                "member_steps": [candidate["steps"]],
            })
    out: list[dict[str, Any]] = []
    for row in classes:
        spectrum = laplacian_spectrum(row["representative_edges"])
        out.append({
            "representative": row["representative"],
            "members": row["members"],
            "member_steps": row["member_steps"],
            "isomorphism_class_size_in_declared_labels": len(row["members"]),
            "laplacian_spectrum_rounded": spectrum,
            "spectrum_signature": json.dumps(spectrum),
        })
    return sorted(out, key=lambda row: row["representative"])


def spectral_quotient_witness() -> dict[str, Any]:
    candidates = candidate_edge_sets()
    classes = quotient_classes(candidates)
    buckets: dict[str, list[str]] = {}
    for row in classes:
        buckets.setdefault(row["spectrum_signature"], []).append(row["representative"])
    collisions = [names for names in buckets.values() if len(names) > 1]
    return {
        "candidate_class": "P2778 declared class: torus_4x4 plus connected 16-node 4-regular circulant Cayley graphs C16({±a,±b})",
        "labeled_candidate_count": len(candidates),
        "isomorphism_class_count": len(classes),
        "distinct_laplacian_spectrum_count_after_quotient": len(buckets),
        "nonisomorphic_cospectral_collision_count": len(collisions),
        "nonisomorphic_cospectral_collision_representatives": collisions,
        "full_spectrum_injective_on_declared_quotient": len(collisions) == 0 and len(buckets) == len(classes),
        "torus_4x4_representative_present": any(row["representative"] == "torus_4x4" for row in classes),
        "quotient_rows": classes,
        "finite_positive_statement": "After isomorphism quotienting the P2778 declared class, the full Laplacian spectrum has no nonisomorphic collision.",
    }


def acceptance_matrix(witness: dict[str, Any], p2778: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2778_max_symmetry_no_go_present": p2778.get("status") == "P2778_MAX_SYMMETRY_16NODE_GEOMETRY_SOURCE_AUDIT_NO_CLOSURE",
        "declared_16node_class_audited": True,
        "isomorphism_quotient_performed": True,
        "full_spectrum_injective_on_declared_quotient": witness["full_spectrum_injective_on_declared_quotient"],
        "strict_nadsoliton_spectral_source_law_exported": False,
        "graph_class_globality_beyond_declared_circulants_exported": False,
        "target_spectrum_or_variational_minimizer_exported_before_testing": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_declared_class_full_spectrum_uniqueness_witness": facts["full_spectrum_injective_on_declared_quotient"],
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_kernel_full_expression_theorem": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "The full spectrum is injective on the declared 16-node quotient, but the class and target are still externally declared; no strict spectral action/source law, global graph-class theorem, or K/L_total variational coupling is exported.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    witness = payload["full_spectrum_16node_quotient_witness"]
    lines = [
        "# P2779/S1729 16-node circulant full-spectrum quotient audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Quotient spectral result",
        f"- labeled_candidate_count={witness['labeled_candidate_count']}",
        f"- isomorphism_class_count={witness['isomorphism_class_count']}",
        f"- distinct_laplacian_spectrum_count_after_quotient={witness['distinct_laplacian_spectrum_count_after_quotient']}",
        f"- nonisomorphic_cospectral_collision_count={witness['nonisomorphic_cospectral_collision_count']}",
        f"- full_spectrum_injective_on_declared_quotient={witness['full_spectrum_injective_on_declared_quotient']}",
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
    p2778 = read_json(P2778)
    witness = spectral_quotient_witness()
    acceptance = acceptance_matrix(witness, p2778)
    payload = {
        "status": "P2779_16NODE_CIRCULANT_FULL_SPECTRUM_QUOTIENT_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2778": sha(P2778)},
        "input_statuses": {"P2778": p2778.get("status")},
        "audited_question": "Does the full Laplacian spectrum collide between nonisomorphic geometries on the P2778 declared 16-node quotient class?",
        "full_spectrum_16node_quotient_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Treat P2779 only as a declared-class positive spectral uniqueness witness.  The next honest move is exactly one of: export a strict nadsoliton spectral action/source law that fixes the admissible graph class, target spectrum, and variational coupling before testing; or stress-test full-spectrum uniqueness outside the circulant/Cayley class by adding one broader 16-node 4-regular non-circulant candidate family with an isomorphism quotient and cospectral audit.  Otherwise preserve the P2697-P2779 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2779/S1729 16-node circulant full-spectrum quotient audit", "## P2779/S1729 16-node circulant full-spectrum quotient audit\n\n`P2779/S1729` performs an isomorphism quotient and full graph-Laplacian-spectrum audit on the P2778 declared 16-node class: `torus_4x4` plus connected 4-regular circulant Cayley candidates.  The 19 labeled candidates collapse to 6 isomorphism classes, and the full Laplacian spectrum is injective on those 6 classes with zero nonisomorphic cospectral collisions.  This is only a declared-class spectral uniqueness witness: no strict nadsoliton spectral source law, global graph-class theorem, kernel full-expression theorem, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2779/S1729 spectral quotient Ltotal guard", "## P2779/S1729 spectral quotient Ltotal guard\n\n`P2779/S1729` adds no variational source term.  Full-spectrum uniqueness holds on the isomorphism quotient of the declared 16-node circulant class, but the class, target, and spectrum are not sourced by a nonproxy `K`/`L_total` variational principle; therefore it cannot promote role-bearing `L_total` or canonical nadsoliton geometry.\n")
    append_once(AGENTS, "Current 16-node circulant full-spectrum quotient audit guardrail (P2779/S1729, 2026-06-15)", "## Current 16-node circulant full-spectrum quotient audit guardrail (P2779/S1729, 2026-06-15)\n\n- P2779 quotients the P2778 declared 16-node class by graph isomorphism and audits the full graph-Laplacian spectrum.\n- The 19 labeled candidates collapse to 6 isomorphism classes, and no nonisomorphic cospectral collision appears on this declared quotient.\n- Do not promote this declared-class uniqueness witness to canonical nadsoliton geometry, strict spectral source law, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  A next admissible move must export a strict spectral action/source law before testing, or broaden the 16-node audit outside the circulant/Cayley class with a fresh quotient and collision test.\n")
    return payload


if __name__ == "__main__":
    main()
