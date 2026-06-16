#!/usr/bin/env python3
"""P2791/S1741: eight-class orbit lower-bound certificate.

P2790 supplied one explicit eighth connected 16-node 4-regular witness outside
the seven local representatives.  This follow-up converts that fact into a
finite group lower-bound certificate: verify all eight representatives are
pairwise non-isomorphic, compute exact automorphism-group sizes, and sum the
disjoint S_16 orbit sizes.

This proves a concrete lower bound on the labeled connected 16-node 4-regular
class covered by the eight certified orbits.  It is still not a full generator
and does not export a strict spectral source law.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2786_s1736_graph6_provenance_toolchain_gate import graph6_decode
from p2787_s1737_small_canonical_generator_pipeline_audit import isomorphic
from p2788_s1738_complement_duality_exact_spectral_certificate import normalize_edges
from p2789_s1739_orbit_stabilizer_exact_quotient_certificate import automorphism_count
from p2790_s1740_eighth_16node_witness_no_exhaustion_certificate import EIGHTH_WITNESS_EDGES

GEN = ROOT / "generated"
P2786 = GEN / "p2786_s1736_graph6_provenance_toolchain_gate.json"
P2790 = GEN / "p2790_s1740_eighth_16node_witness_no_exhaustion_certificate.json"
OUT = GEN / "p2791_s1741_eight_class_orbit_lower_bound_certificate.json"
MD = GEN / "p2791_s1741_eight_class_orbit_lower_bound_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 16
NEGATIVE_EXPORT_FLAGS = [
    "canonical_16node_generator_certified",
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


def certified_rows(p2786: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for row in p2786.get("graph6_provenance_witness", {}).get("provenance_rows", []):
        rows.append({
            "label": row["representative"],
            "source": "P2786_local_representative",
            "edges": normalize_edges(graph6_decode(row["graph6"])),
        })
    rows.append({
        "label": "p2790_eighth_witness",
        "source": "P2790_explicit_no_exhaustion_witness",
        "edges": EIGHTH_WITNESS_EDGES,
    })
    return rows


def orbit_lower_bound_witness(p2786: dict[str, Any]) -> dict[str, Any]:
    rows = certified_rows(p2786)
    representative_rows = []
    for row in rows:
        aut_size = automorphism_count(row["edges"], N)
        representative_rows.append({
            "label": row["label"],
            "source": row["source"],
            "edge_count": len(row["edges"]),
            "automorphism_group_size": aut_size,
            "orbit_size_by_orbit_stabilizer": math.factorial(N) // aut_size,
        })
    pair_rows = []
    for left_index, left in enumerate(rows):
        for right_index in range(left_index + 1, len(rows)):
            right = rows[right_index]
            pair_rows.append({
                "left": left["label"],
                "right": right["label"],
                "isomorphic": isomorphic(left["edges"], right["edges"], N),
            })
    lower_bound = sum(row["orbit_size_by_orbit_stabilizer"] for row in representative_rows)
    return {
        "source_class": "Seven P2786/P2785 local representatives plus the explicit P2790 eighth witness.",
        "representative_count": len(rows),
        "pair_count": len(pair_rows),
        "representative_rows": representative_rows,
        "pairwise_isomorphism_rows": pair_rows,
        "all_28_pairs_nonisomorphic": all(not row["isomorphic"] for row in pair_rows),
        "certified_disjoint_labeled_orbit_lower_bound": lower_bound,
        "largest_single_orbit_label": max(representative_rows, key=lambda row: row["orbit_size_by_orbit_stabilizer"])["label"],
        "finite_certificate_statement": "Eight certified pairwise non-isomorphic connected 16-node 4-regular representatives give a disjoint labeled-orbit lower bound of 13,463,256,807,000.",
    }


def acceptance_matrix(witness: dict[str, Any], p2786: dict[str, Any], p2790: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2786_graph6_gate_present": p2786.get("status") == "P2786_GRAPH6_PROVENANCE_TOOLCHAIN_GATE_NO_CLOSURE",
        "p2790_eighth_witness_present": p2790.get("status") == "P2790_EIGHTH_16NODE_WITNESS_NO_EXHAUSTION_CERTIFICATE_NO_CLOSURE",
        "eight_representatives_checked": witness["representative_count"] == 8,
        "all_28_pairs_checked_nonisomorphic": witness["pair_count"] == 28 and witness["all_28_pairs_nonisomorphic"],
        "positive_disjoint_labeled_orbit_lower_bound": witness["certified_disjoint_labeled_orbit_lower_bound"] == 13463256807000,
        "canonical_16node_generator_certified": False,
        "strict_nadsoliton_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_eight_class_orbit_lower_bound_certificate": all(facts[key] for key in [
            "p2786_graph6_gate_present",
            "p2790_eighth_witness_present",
            "eight_representatives_checked",
            "all_28_pairs_checked_nonisomorphic",
            "positive_disjoint_labeled_orbit_lower_bound",
        ]),
        "accepted_as_full_16node_canonical_generator_certificate": False,
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "The eight disjoint orbits provide an exact lower bound and strengthen the no-exhaustion result, but a lower bound is not a full connected 16-node 4-regular generator and exports no strict K/L_total spectral source law.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    w = payload["orbit_lower_bound_witness"]
    lines = [
        "# P2791/S1741 eight-class orbit lower-bound certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Exact eight-class lower-bound result",
        f"- representative_count={w['representative_count']}",
        f"- pair_count={w['pair_count']}",
        f"- all_28_pairs_nonisomorphic={w['all_28_pairs_nonisomorphic']}",
        f"- certified_disjoint_labeled_orbit_lower_bound={w['certified_disjoint_labeled_orbit_lower_bound']}",
        f"- largest_single_orbit_label={w['largest_single_orbit_label']}",
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
    p2786 = read_json(P2786)
    p2790 = read_json(P2790)
    witness = orbit_lower_bound_witness(p2786)
    acceptance = acceptance_matrix(witness, p2786, p2790)
    payload = {
        "status": "P2791_EIGHT_CLASS_ORBIT_LOWER_BOUND_CERTIFICATE_NO_CLOSURE",
        "input_hashes": {"P2786": sha(P2786), "P2790": sha(P2790)},
        "input_statuses": {"P2786": p2786.get("status"), "P2790": p2790.get("status")},
        "audited_question": "What exact disjoint labeled-orbit lower bound follows after adding the P2790 eighth witness to the seven local 16-node representatives?",
        "orbit_lower_bound_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Use P2791 only as an exact lower-bound/no-exhaustion certificate.  The next honest move is exactly one of: supply/import an actual certified full connected 16-node 4-regular generator artifact/toolchain with graph6/hash provenance and run full exact quotient/charpoly/complement/orbit auditing; or export a strict nadsoliton spectral action/source law fixing the admissible class, target spectrum, and K/L_total coupling before testing.  Otherwise preserve the P2697-P2791 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2791/S1741 eight-class orbit lower-bound certificate", "## P2791/S1741 eight-class orbit lower-bound certificate\n\n`P2791/S1741` converts the P2790 eighth witness into an exact disjoint-orbit lower-bound certificate.  The seven P2786/P2785 local representatives plus the P2790 witness are pairwise non-isomorphic across all 28 pairs, and their exact `S_16` orbit sizes sum to 13,463,256,807,000 labeled connected 16-node 4-regular graphs covered by these eight certified orbits.  This is still a lower bound/no-exhaustion certificate, not a full connected 16-node 4-regular generator, not a strict spectral source law, and not a `K`/`L_total` variational coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2791/S1741 eight-class lower-bound Ltotal guard", "## P2791/S1741 eight-class lower-bound Ltotal guard\n\n`P2791/S1741` adds no variational source term.  Exact disjoint-orbit lower-bound arithmetic strengthens finite graph bookkeeping, but a labeled-orbit lower bound is not a sourced nonproxy `K`/`L_total` spectral action, not a canonical geometry theorem, and not a full 16-node generator.\n")
    append_once(AGENTS, "Current eight-class orbit lower-bound guardrail (P2791/S1741, 2026-06-16)", "## Current eight-class orbit lower-bound guardrail (P2791/S1741, 2026-06-16)\n\n- P2791 verifies that the seven P2786/P2785 local representatives plus the P2790 eighth witness are pairwise non-isomorphic across all 28 pairs and computes an exact disjoint labeled-orbit lower bound of 13,463,256,807,000 under `S_16`.\n- This strengthens no-exhaustion bookkeeping, but it is not the required full connected 16-node 4-regular generator/toolchain and does not source geometry from `K`/`L_total`.\n- Do not promote the eight-class lower bound to canonical geometry, strict spectral source law, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  A next admissible move must supply an actual certified full 16-node generator artifact/toolchain or export a strict spectral action/source law before testing.\n")
    return payload


if __name__ == "__main__":
    main()
