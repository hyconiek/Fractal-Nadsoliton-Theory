#!/usr/bin/env python3
"""P2783/S1733: seven-class quotient integrity certificate.

The user asked whether the seven P2781 quotient classes might secretly collide
as isomorphic/nonisomorphic classes.  This audit does not extend the graph class.
Instead it verifies the integrity of the existing seven-class quotient by two
independent finite checks:

1. direct pairwise graph-isomorphism backtracking on the seven representatives;
2. pairwise full-Laplacian-spectrum signatures, which must agree for isomorphic
   graphs and are distinct here for all 21 representative pairs.

The certificate answers the local question positively: within the P2781 expanded
class, the seven representatives are pairwise nonisomorphic and pairwise
spectrally distinct.  This is still not canonical geometry or a strict spectral
source law.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2774_s1724_entropy_laplacian_trace_geometry_degeneracy import all_pairs_distances, distance_histogram
from p2779_s1729_16node_circulant_full_spectrum_quotient_audit import isomorphic, laplacian_spectrum, quotient_classes
from p2781_s1731_enumerated_two_layer_c8_spectrum_collision_audit import all_two_layer_c8_candidates
from p2778_s1728_max_symmetry_16node_geometry_source_audit import candidate_edge_sets

GEN = ROOT / "generated"
P2781 = GEN / "p2781_s1731_enumerated_two_layer_c8_spectrum_collision_audit.json"
P2782 = GEN / "p2782_s1732_bipartite_regular_enumerator_scale_obstruction.json"
OUT = GEN / "p2783_s1733_seven_class_quotient_integrity_certificate.json"
MD = GEN / "p2783_s1733_seven_class_quotient_integrity_certificate.md"
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


def representative_candidates() -> list[dict[str, Any]]:
    candidates = candidate_edge_sets() + all_two_layer_c8_candidates()
    by_name = {row["geometry"]: row for row in candidates}
    quotient = quotient_classes(candidates)
    reps: list[dict[str, Any]] = []
    for row in quotient:
        source = by_name[row["representative"]]
        spectrum = laplacian_spectrum(source["edges"])
        distances = all_pairs_distances(source["edges"])
        reps.append({
            "representative": row["representative"],
            "member_count": len(row["members"]),
            "members": row["members"],
            "edges": source["edges"],
            "laplacian_spectrum_rounded": spectrum,
            "spectrum_signature": json.dumps(spectrum),
            "distance_histogram": {str(k): v for k, v in distance_histogram(distances).items()},
        })
    return sorted(reps, key=lambda row: row["representative"])


def integrity_witness() -> dict[str, Any]:
    reps = representative_candidates()
    pair_rows: list[dict[str, Any]] = []
    for left, right in itertools.combinations(reps, 2):
        direct_iso = isomorphic(left["edges"], right["edges"])
        same_spectrum = left["spectrum_signature"] == right["spectrum_signature"]
        same_distance_histogram = left["distance_histogram"] == right["distance_histogram"]
        pair_rows.append({
            "left": left["representative"],
            "right": right["representative"],
            "direct_isomorphic": direct_iso,
            "same_full_laplacian_spectrum": same_spectrum,
            "same_distance_histogram": same_distance_histogram,
            "separation_certificate": "direct_nonisomorphism_and_distinct_spectrum" if (not direct_iso and not same_spectrum) else "collision_or_unresolved",
        })
    return {
        "source_class": "P2781 expanded class: P2779 base plus all two-C8-layer cross-matching shift pairs",
        "representative_count": len(reps),
        "pair_count": len(pair_rows),
        "representative_rows": [
            {key: value for key, value in row.items() if key != "edges"}
            for row in reps
        ],
        "pair_rows": pair_rows,
        "direct_isomorphism_collision_count": sum(1 for row in pair_rows if row["direct_isomorphic"]),
        "full_spectrum_collision_count": sum(1 for row in pair_rows if row["same_full_laplacian_spectrum"]),
        "distance_histogram_collision_count": sum(1 for row in pair_rows if row["same_distance_histogram"]),
        "all_representatives_pairwise_nonisomorphic": all(not row["direct_isomorphic"] for row in pair_rows),
        "all_representatives_pairwise_spectrally_distinct": all(not row["same_full_laplacian_spectrum"] for row in pair_rows),
        "finite_certificate_statement": "The seven P2781 representatives have zero direct isomorphism collisions and zero full-spectrum collisions across all 21 pairs.",
    }


def acceptance_matrix(witness: dict[str, Any], p2781: dict[str, Any], p2782: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2781_seven_class_witness_present": p2781.get("status") == "P2781_ENUMERATED_TWO_LAYER_C8_SPECTRUM_COLLISION_AUDIT_NO_CLOSURE",
        "p2782_enumerator_scale_boundary_present": p2782.get("status") == "P2782_BIPARTITE_REGULAR_ENUMERATOR_SCALE_OBSTRUCTION_NO_CLOSURE",
        "all_21_pairs_checked_by_direct_isomorphism": witness["pair_count"] == 21,
        "all_representatives_pairwise_nonisomorphic": witness["all_representatives_pairwise_nonisomorphic"],
        "all_representatives_pairwise_spectrally_distinct": witness["all_representatives_pairwise_spectrally_distinct"],
        "strict_nadsoliton_spectral_source_law_exported": False,
        "canonical_full_16node_4regular_generator_supplied": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_p2781_quotient_integrity_certificate": all(facts[key] for key in [
            "all_21_pairs_checked_by_direct_isomorphism",
            "all_representatives_pairwise_nonisomorphic",
            "all_representatives_pairwise_spectrally_distinct",
        ]),
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "The seven-class quotient is internally clean, but it remains the P2781 local expanded class.  No strict spectral source law, full canonical 16-node generator, or K/L_total coupling is exported.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    witness = payload["seven_class_integrity_witness"]
    lines = [
        "# P2783/S1733 seven-class quotient integrity certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Pairwise integrity result",
        f"- representative_count={witness['representative_count']}",
        f"- pair_count={witness['pair_count']}",
        f"- direct_isomorphism_collision_count={witness['direct_isomorphism_collision_count']}",
        f"- full_spectrum_collision_count={witness['full_spectrum_collision_count']}",
        f"- distance_histogram_collision_count={witness['distance_histogram_collision_count']}",
        f"- all_representatives_pairwise_nonisomorphic={witness['all_representatives_pairwise_nonisomorphic']}",
        f"- all_representatives_pairwise_spectrally_distinct={witness['all_representatives_pairwise_spectrally_distinct']}",
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
    p2781 = read_json(P2781)
    p2782 = read_json(P2782)
    witness = integrity_witness()
    acceptance = acceptance_matrix(witness, p2781, p2782)
    payload = {
        "status": "P2783_SEVEN_CLASS_QUOTIENT_INTEGRITY_CERTIFICATE_NO_CLOSURE",
        "input_hashes": {"P2781": sha(P2781), "P2782": sha(P2782)},
        "input_statuses": {"P2781": p2781.get("status"), "P2782": p2782.get("status")},
        "audited_question": "Do the seven P2781 quotient classes secretly collide by isomorphism or full-Laplacian spectrum?",
        "seven_class_integrity_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Use P2783 as a local quotient-integrity certificate only.  The next honest move is exactly one of: supply a canonical-generation theorem/tool certificate for the full connected 16-node 4-regular graph class and rerun this pairwise integrity plus full-spectrum collision audit; or export a strict nadsoliton spectral action/source law fixing the admissible class, target spectrum, and K/L_total coupling before testing.  Otherwise preserve the P2697-P2783 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2783/S1733 seven-class quotient integrity certificate", "## P2783/S1733 seven-class quotient integrity certificate\n\n`P2783/S1733` answers the local integrity question for the seven P2781 quotient classes.  It reruns the P2781 expanded class, extracts the seven representatives, and checks all 21 representative pairs by direct graph-isomorphism backtracking and full graph-Laplacian spectrum.  The result has zero direct isomorphism collisions and zero full-spectrum collisions, so the seven local classes are pairwise nonisomorphic and spectrally distinct.  This is only a local quotient-integrity certificate: no strict nadsoliton spectral source law, full canonical 16-node graph-class theorem, kernel full-expression theorem, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2783/S1733 seven-class quotient Ltotal guard", "## P2783/S1733 seven-class quotient Ltotal guard\n\n`P2783/S1733` adds no variational source term.  The seven P2781 representatives are pairwise nonisomorphic and spectrally distinct, but this remains a local quotient-integrity certificate rather than a sourced nonproxy `K`/`L_total` spectral action; therefore it cannot promote role-bearing `L_total` or canonical nadsoliton geometry.\n")
    append_once(AGENTS, "Current seven-class quotient integrity guardrail (P2783/S1733, 2026-06-15)", "## Current seven-class quotient integrity guardrail (P2783/S1733, 2026-06-15)\n\n- P2783 audits whether the seven P2781 quotient classes secretly collide by isomorphism or by full graph-Laplacian spectrum.\n- All 21 representative pairs are directly checked: there are zero direct isomorphism collisions and zero full-spectrum collisions, so the seven local classes are pairwise nonisomorphic and spectrally distinct.\n- Do not promote this local quotient-integrity certificate to canonical geometry, strict spectral source law, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  A next admissible move must supply a canonical-generation theorem/tool certificate for the full graph class, or export a strict spectral action/source law before testing.\n")
    return payload


if __name__ == "__main__":
    main()
