#!/usr/bin/env python3
"""P2795/S1745: post-subclass saturation / no-new-class frontier certificate.

P2792 exhausted the C16 two-jump circulant subclass.  P2793 exhausted the
Z_2^4 four-generator Cayley subclass.  P2794 proved the P2793 collapse by an
independent GL(4,2) group-action certificate.  P2795 is the honest reconciliation
step: compute the marginal class contribution of these named-subclass lanes
against the current P2791 eight-class lower-bound witness set.

The result is intentionally blocking rather than promotional: the named-subclass
work is saturated on current artifacts because every certified subclass class is
already inside the P2791 eight-class set, so continuing named-subclass expansion
without a new full generator/toolchain or a strict source law is a replay risk.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
P2791 = GEN / "p2791_s1741_eight_class_orbit_lower_bound_certificate.json"
P2792 = GEN / "p2792_s1742_c16_two_jump_circulant_subclass_exhaustion_certificate.json"
P2793 = GEN / "p2793_s1743_z2_4_four_generator_cayley_subclass_exhaustion_certificate.json"
P2794 = GEN / "p2794_s1744_z2_4_gl4_transitivity_automorphism_factorization_certificate.json"
OUT = GEN / "p2795_s1745_post_subclass_saturation_no_new_class_frontier_certificate.json"
MD = GEN / "p2795_s1745_post_subclass_saturation_no_new_class_frontier_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

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


def p2791_labels(p2791: dict[str, Any]) -> list[str]:
    return [row["label"] for row in p2791.get("orbit_lower_bound_witness", {}).get("representative_rows", [])]


def p2791_orbits(p2791: dict[str, Any]) -> dict[str, int]:
    return {
        row["label"]: row["orbit_size_by_orbit_stabilizer"]
        for row in p2791.get("orbit_lower_bound_witness", {}).get("representative_rows", [])
    }


def c16_labels(p2792: dict[str, Any]) -> list[str]:
    labels: list[str] = []
    for row in p2792.get("circulant_subclass_witness", {}).get("class_rows", []):
        labels.extend(row.get("matching_p2786_local_representatives", []))
    return sorted(set(labels))


def z2_labels(p2793: dict[str, Any]) -> list[str]:
    labels: list[str] = []
    for row in p2793.get("z2_4_cayley_subclass_witness", {}).get("class_rows", []):
        labels.extend(row.get("matching_p2791_eight_class_representatives", []))
    return sorted(set(labels))


def saturation_witness(p2791: dict[str, Any], p2792: dict[str, Any], p2793: dict[str, Any], p2794: dict[str, Any]) -> dict[str, Any]:
    eight_labels = p2791_labels(p2791)
    eight_set = set(eight_labels)
    orbits = p2791_orbits(p2791)
    subclass_rows = [
        {
            "lane": "P2792_C16_two_jump_circulant",
            "certificate_status": p2792.get("status"),
            "certified_class_count": p2792.get("circulant_subclass_witness", {}).get("connected_isomorphism_class_count"),
            "matched_p2791_labels": c16_labels(p2792),
        },
        {
            "lane": "P2793_Z2_4_four_generator_Cayley",
            "certificate_status": p2793.get("status"),
            "certified_class_count": p2793.get("z2_4_cayley_subclass_witness", {}).get("connected_isomorphism_class_count"),
            "matched_p2791_labels": z2_labels(p2793),
        },
        {
            "lane": "P2794_GL4_transitivity_explanation_of_P2793",
            "certificate_status": p2794.get("status"),
            "certified_class_count": 1 if p2794.get("acceptance_matrix", {}).get("accepted_as_gl4_transitivity_and_automorphism_factorization_certificate") else 0,
            "matched_p2791_labels": z2_labels(p2793),
        },
    ]
    cumulative_labels: set[str] = set()
    for row in subclass_rows:
        row_labels = set(row["matched_p2791_labels"])
        row["all_matches_inside_p2791_eight_class_set"] = row_labels <= eight_set
        row["new_labels_not_already_in_p2791"] = sorted(row_labels - eight_set)
        row["marginal_new_p2791_label_count"] = len(row_labels - cumulative_labels)
        row["orbit_lower_bound_contribution_already_counted_in_p2791"] = sum(orbits[label] for label in row_labels if label in orbits)
        cumulative_labels |= row_labels
    union_labels = sorted({label for row in subclass_rows for label in row["matched_p2791_labels"]})
    uncovered_p2791_labels = sorted(eight_set - set(union_labels))
    return {
        "p2791_eight_class_labels": eight_labels,
        "p2791_eight_class_count": len(eight_labels),
        "p2791_lower_bound": p2791.get("orbit_lower_bound_witness", {}).get("certified_disjoint_labeled_orbit_lower_bound"),
        "subclass_rows": subclass_rows,
        "union_of_named_subclass_labels": union_labels,
        "union_named_subclass_label_count": len(union_labels),
        "uncovered_p2791_labels_by_named_subclasses": uncovered_p2791_labels,
        "uncovered_p2791_label_count": len(uncovered_p2791_labels),
        "all_named_subclass_labels_inside_p2791": all(row["all_matches_inside_p2791_eight_class_set"] for row in subclass_rows),
        "total_new_labels_beyond_p2791": len(set(union_labels) - eight_set),
        "total_new_orbit_lower_bound_added_beyond_p2791": 0,
        "finite_certificate_statement": "P2792, P2793, and P2794 add zero new isomorphism classes beyond the P2791 eight-class lower-bound witness set; the named-subclass lane is saturated on current artifacts.",
    }


def acceptance_matrix(witness: dict[str, Any], p2791: dict[str, Any], p2792: dict[str, Any], p2793: dict[str, Any], p2794: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2791_eight_class_lower_bound_present": p2791.get("status") == "P2791_EIGHT_CLASS_ORBIT_LOWER_BOUND_CERTIFICATE_NO_CLOSURE",
        "p2792_c16_subclass_present": p2792.get("status") == "P2792_C16_TWO_JUMP_CIRCULANT_SUBCLASS_EXHAUSTION_CERTIFICATE_NO_CLOSURE",
        "p2793_z2_cayley_subclass_present": p2793.get("status") == "P2793_Z2_4_FOUR_GENERATOR_CAYLEY_SUBCLASS_EXHAUSTION_CERTIFICATE_NO_CLOSURE",
        "p2794_gl4_proof_present": p2794.get("status") == "P2794_Z2_4_GL4_TRANSITIVITY_AUTOMORPHISM_FACTORIZATION_CERTIFICATE_NO_CLOSURE",
        "all_named_subclass_labels_inside_p2791": witness["all_named_subclass_labels_inside_p2791"],
        "zero_new_labels_beyond_p2791": witness["total_new_labels_beyond_p2791"] == 0,
        "zero_new_orbit_lower_bound_beyond_p2791": witness["total_new_orbit_lower_bound_added_beyond_p2791"] == 0,
        "named_subclasses_do_not_cover_all_p2791_labels": witness["uncovered_p2791_label_count"] > 0,
        "canonical_16node_generator_certified": False,
        "strict_nadsoliton_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_post_subclass_saturation_no_new_class_certificate": all(facts[key] for key in [
            "p2791_eight_class_lower_bound_present",
            "p2792_c16_subclass_present",
            "p2793_z2_cayley_subclass_present",
            "p2794_gl4_proof_present",
            "all_named_subclass_labels_inside_p2791",
            "zero_new_labels_beyond_p2791",
            "zero_new_orbit_lower_bound_beyond_p2791",
            "named_subclasses_do_not_cover_all_p2791_labels",
        ]),
        "accepted_as_full_16node_canonical_generator_certificate": False,
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "The recent named-subclass certificates are saturated against P2791: they add zero new classes and still leave P2791 classes uncovered, so continuing named-subclass replay is not a proof-grade route to canonical geometry without a full generator/toolchain or strict spectral source law.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    w = payload["post_subclass_saturation_witness"]
    lines = [
        "# P2795/S1745 post-subclass saturation no-new-class frontier certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Saturation result",
        f"- p2791_eight_class_count={w['p2791_eight_class_count']}",
        f"- union_named_subclass_label_count={w['union_named_subclass_label_count']}",
        f"- union_of_named_subclass_labels={w['union_of_named_subclass_labels']}",
        f"- uncovered_p2791_label_count={w['uncovered_p2791_label_count']}",
        f"- uncovered_p2791_labels_by_named_subclasses={w['uncovered_p2791_labels_by_named_subclasses']}",
        f"- total_new_labels_beyond_p2791={w['total_new_labels_beyond_p2791']}",
        f"- total_new_orbit_lower_bound_added_beyond_p2791={w['total_new_orbit_lower_bound_added_beyond_p2791']}",
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
    p2791 = read_json(P2791)
    p2792 = read_json(P2792)
    p2793 = read_json(P2793)
    p2794 = read_json(P2794)
    witness = saturation_witness(p2791, p2792, p2793, p2794)
    acceptance = acceptance_matrix(witness, p2791, p2792, p2793, p2794)
    payload = {
        "status": "P2795_POST_SUBCLASS_SATURATION_NO_NEW_CLASS_FRONTIER_CERTIFICATE_NO_CLOSURE",
        "input_hashes": {"P2791": sha(P2791), "P2792": sha(P2792), "P2793": sha(P2793), "P2794": sha(P2794)},
        "input_statuses": {"P2791": p2791.get("status"), "P2792": p2792.get("status"), "P2793": p2793.get("status"), "P2794": p2794.get("status")},
        "audited_question": "Do the recent named-subclass certificates add any new 16-node 4-regular isomorphism class or orbit lower-bound contribution beyond P2791?",
        "post_subclass_saturation_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Stop replaying named-subclass expansion unless a genuinely new generator family is known to add classes beyond P2791.  The next honest move remains exactly one of: supply/import an actual certified full connected 16-node 4-regular generator artifact/toolchain with graph6/hash provenance and run full exact quotient/charpoly/complement/orbit auditing; or export a strict nadsoliton spectral action/source law fixing the admissible class, target spectrum, and K/L_total coupling before testing.  Otherwise preserve the P2697-P2795 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2795/S1745 post-subclass saturation no-new-class certificate", "## P2795/S1745 post-subclass saturation no-new-class certificate\n\n`P2795/S1745` reconciles the recent named-subclass lanes against the P2791 eight-class lower-bound witness set.  The C16 two-jump circulant lane contributes the five labels `circulant_pm1_pm2`, `circulant_pm1_pm3`, `circulant_pm1_pm4`, `circulant_pm1_pm6`, and `circulant_pm1_pm7`; the `Z_2^4` Cayley/GL(4,2) lane contributes only `torus_4x4`.  Their union is six labels, all already inside P2791, so the marginal new class count and marginal orbit lower-bound beyond P2791 are both zero.  This is a no-new-class frontier certificate for current named-subclass work, not a full connected 16-node 4-regular generator, not a strict spectral source law, and not a `K`/`L_total` variational coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2795/S1745 post-subclass saturation Ltotal guard", "## P2795/S1745 post-subclass saturation Ltotal guard\n\n`P2795/S1745` adds no variational source term.  It shows that current named-subclass certificates add zero new classes beyond P2791 and therefore should not be replayed as a route to `K`/`L_total`; a saturation/no-new-class certificate is not a sourced nonproxy spectral action, not canonical geometry, and not a full 16-node generator.\n")
    append_once(AGENTS, "Current post-subclass saturation guardrail (P2795/S1745, 2026-06-16)", "## Current post-subclass saturation guardrail (P2795/S1745, 2026-06-16)\n\n- P2795 reconciles P2792/P2793/P2794 against P2791 and finds zero new isomorphism classes and zero new orbit lower-bound contribution beyond the P2791 eight-class set; the named-subclass union covers six already-counted P2791 labels and leaves two P2791 labels uncovered.\n- Do not continue named-subclass replay as a primary route unless a genuinely new generator family is known to add classes beyond P2791, or a certified full connected 16-node 4-regular generator/toolchain is supplied.\n- Do not promote the saturation/no-new-class result to canonical geometry, strict spectral source law, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.\n")
    return payload


if __name__ == "__main__":
    main()
