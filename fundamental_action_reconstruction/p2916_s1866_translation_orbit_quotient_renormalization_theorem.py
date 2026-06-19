#!/usr/bin/env python3
"""P2916/S1866: translation-orbit quotient renormalization theorem.

P2915 found three arithmetic quotient candidates for the P2914 12-vs-144
normalization mismatch.  P2916 supplies exactly one finite quotient theorem that
is stricter than a mere candidate: quotient directed Z12 edges by the diagonal
translation action (i,j) -> (i+a,j+a).  The quotient invariant is the relative
displacement d = j-i mod 12.

For a lay reader: the 12 is the number of nadsoliton phase/sites on the Z12
circle.  The 144 is the number of directed site-to-site relations because every
one of 12 starting sites can point to every one of 12 ending sites.  Counting all
144 edges overcounts the same relative jump at 12 translated basepoints.  The
translation-orbit quotient keeps only the relative jump d, leaving 12 quotient
classes.  This gives the exact renormalization factor 12.

This is a real finite quotient theorem, not L_total closure: it proves the
translation-orbit quotient and its 12-fold renormalization, but still does not
export a nonzero Gamma_9_5 source, continuum field-variable provenance, EOM,
Hamiltonian, role transfer, or ToE closure.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
import hashlib
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2915 = GEN / "p2915_s1865_gamma_edge_quotient_measure_bridge_candidate_gate.json"
OUT = GEN / "p2916_s1866_translation_orbit_quotient_renormalization_theorem.json"
MD = GEN / "p2916_s1866_translation_orbit_quotient_renormalization_theorem.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
N = 12
EDGE_COUNT = N * N


def f(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def edges() -> list[tuple[int, int]]:
    return [(i, j) for i in range(N) for j in range(N)]


def translate(edge: tuple[int, int], shift: int) -> tuple[int, int]:
    return ((edge[0] + shift) % N, (edge[1] + shift) % N)


def displacement(edge: tuple[int, int]) -> int:
    return (edge[1] - edge[0]) % N


def source_label(edge: tuple[int, int]) -> int:
    return edge[0]


def target_label(edge: tuple[int, int]) -> int:
    return edge[1]


def label_invariance_failures(label_name: str) -> list[dict[str, Any]]:
    label_fn = {"source": source_label, "target": target_label, "displacement": displacement}[label_name]
    failures: list[dict[str, Any]] = []
    for edge in edges():
        base = label_fn(edge)
        for shift in range(N):
            shifted = translate(edge, shift)
            if label_fn(shifted) != base:
                failures.append({"edge": [edge[0], edge[1]], "shift": shift, "base_label": base, "shifted_edge": [shifted[0], shifted[1]], "shifted_label": label_fn(shifted)})
                break
    return failures


def displacement_orbits() -> dict[int, list[tuple[int, int]]]:
    orbits: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for edge in edges():
        orbits[displacement(edge)].append(edge)
    return dict(sorted(orbits.items()))


def orbit_rows(site_weight: Fraction) -> list[dict[str, Any]]:
    per_edge_weight = site_weight / N
    rows: list[dict[str, Any]] = []
    for d, orbit in displacement_orbits().items():
        orbit_weight = len(orbit) * per_edge_weight
        rows.append({
            "displacement_class": d,
            "orbit_size": len(orbit),
            "per_edge_renormalized_weight": f(per_edge_weight),
            "orbit_total_weight": f(orbit_weight),
            "sample_edges": [[i, j] for i, j in orbit[:4]],
        })
    return rows


def build_payload(p2915: dict[str, Any]) -> dict[str, Any]:
    site_weight = Fraction(1, N)
    per_edge_weight = site_weight / N
    rows = orbit_rows(site_weight)
    source_failures = label_invariance_failures("source")
    target_failures = label_invariance_failures("target")
    displacement_failures = label_invariance_failures("displacement")
    quotient_total = sum(Fraction(row["orbit_size"], 1) * per_edge_weight for row in rows)
    return {
        "status": "P2916_TRANSLATION_ORBIT_QUOTIENT_RENORMALIZATION_THEOREM_FINITE_EXPORT_NO_LTOTAL",
        "input_hashes": {"P2915": hashlib.sha256(P2915.read_bytes()).hexdigest() if P2915.exists() else None},
        "lay_explanation": {
            "what_12_counts": "the 12 Z12 nadsoliton phase/sites around the finite circle",
            "what_144_counts": "the 12 x 12 directed relations from every start site to every end site",
            "why_mismatch_happens": "summing all directed relations counts each relative jump at 12 translated basepoints",
            "what_quotient_does": "identify translated copies of the same relative jump d=j-i, leaving 12 displacement classes",
        },
        "constructed_theoretical_objects": {
            "theorem_name": "Z12_Diagonal_Translation_Orbit_Edge_Quotient_Renormalization_Theorem",
            "group_action": "a·(i,j) = (i+a, j+a) mod 12",
            "quotient_map": "q(i,j) = j-i mod 12",
            "renormalization_rule": "with site-normalized m=1/12, each directed edge receives m/12=1/144 before quotient summation",
            "orbit_rows": rows,
            "candidate_selector_table": [
                {"candidate": "source-site quotient", "translation_invariant_label_failures": len(source_failures), "selected_by_translation_orbit_theorem": False},
                {"candidate": "target-site quotient", "translation_invariant_label_failures": len(target_failures), "selected_by_translation_orbit_theorem": False},
                {"candidate": "displacement quotient", "translation_invariant_label_failures": len(displacement_failures), "selected_by_translation_orbit_theorem": True},
            ],
            "acceptance_obligation_rows": [
                {"obligation": "diagonal Z12 translation action defined", "passed": True, "exported_strictly": True},
                {"obligation": "displacement quotient invariant under translations", "passed": len(displacement_failures) == 0, "exported_strictly": True},
                {"obligation": "12 orbits of size 12 computed", "passed": len(rows) == N and all(row["orbit_size"] == N for row in rows), "exported_strictly": True},
                {"obligation": "quotient renormalization total equals one", "passed": quotient_total == 1, "exported_strictly": True},
                {"obligation": "strict Gamma_9_5 action-unit source theorem", "passed": False, "exported_strictly": False},
                {"obligation": "continuum field-variable provenance theorem", "passed": False, "exported_strictly": False},
            ],
        },
        "acceptance_matrix": {
            "p2915_rechecked_three_candidate_boundary": p2915.get("acceptance_matrix", {}).get("quotient_candidate_count") == 3,
            "site_count": N,
            "directed_edge_count": EDGE_COUNT,
            "translation_orbit_count": len(rows),
            "all_translation_orbits_size_12": all(row["orbit_size"] == N for row in rows),
            "source_label_translation_failure_count": len(source_failures),
            "target_label_translation_failure_count": len(target_failures),
            "displacement_label_translation_failure_count": len(displacement_failures),
            "selected_quotient": "displacement_quotient_q(i,j)=j-i_mod12",
            "site_normalized_m": f(site_weight),
            "per_edge_renormalized_weight": f(per_edge_weight),
            "quotient_total_weight": f(quotient_total),
            "finite_translation_quotient_theorem_exported": True,
            "strict_gamma_9_5_source_exported": False,
            "continuum_field_variable_provenance_exported": False,
            "accepted_as_nonproxy_ltotal_measure_bridge": False,
        },
        "decision": {
            "positive_witnesses": {
                "one_quotient_selected_by_translation_invariance": True,
                "finite_orbit_quotient_theorem_exported": True,
                "renormalization_factor_12_computed": True,
            },
            "negative_export_flags": {
                "strict_gamma_9_5_source_exported": False,
                "continuum_field_variable_provenance_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2916 proves one finite quotient theorem: diagonal Z12 translation orbits select the displacement quotient q(i,j)=j-i mod 12.  This explains the 12 vs 144 mismatch: 144 directed edges are 12 translated copies of each of 12 relative jumps, so quotient integration uses a factor 1/12 and total weight one.  Source-site and target-site labels are not translation-invariant.  The theorem is finite quotient/renormalization progress only; Gamma_9_5 sourcehood and continuum field-variable provenance remain missing.",
            "next_honest_step": "The next proof-grade move should use this selected displacement quotient and prove exactly one remaining theorem: either a strict nonzero Gamma_9_5 action-unit source coupled to the quotient integral, or a continuum field-variable provenance theorem showing that q(i,j)=j-i quotient classes are the nonproxy variables integrated by L_total.  Without one of those, do not promote to EOM/Hamiltonian/ToE closure.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    acc = payload["acceptance_matrix"]
    lay = payload["lay_explanation"]
    lines = [
        "# P2916/S1866 translation-orbit quotient renormalization theorem",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Lay explanation of the 12 vs 144 mismatch",
        f"- 12 counts: {lay['what_12_counts']}.",
        f"- 144 counts: {lay['what_144_counts']}.",
        f"- mismatch: {lay['why_mismatch_happens']}.",
        f"- quotient repair: {lay['what_quotient_does']}.",
        "",
        "## Finite theorem gate",
        f"- selected quotient: `{acc['selected_quotient']}`",
        f"- translation orbit count: `{acc['translation_orbit_count']}`",
        f"- all orbits size 12: `{acc['all_translation_orbits_size_12']}`",
        f"- source-label translation failures: `{acc['source_label_translation_failure_count']}`",
        f"- target-label translation failures: `{acc['target_label_translation_failure_count']}`",
        f"- displacement-label translation failures: `{acc['displacement_label_translation_failure_count']}`",
        f"- per-edge renormalized weight: `{acc['per_edge_renormalized_weight']}`",
        f"- quotient total weight: `{acc['quotient_total_weight']}`",
        f"- accepted as nonproxy L_total measure bridge: `{acc['accepted_as_nonproxy_ltotal_measure_bridge']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2915))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2916/S1866 translation-orbit quotient renormalization theorem", "## P2916/S1866 translation-orbit quotient renormalization theorem\n\n`P2916/S1866` proves one finite quotient theorem for the P2914/P2915 normalization mismatch: quotient directed edges by the diagonal `Z12` translation action `(i,j) -> (i+a,j+a)`.  The invariant quotient is the displacement `q(i,j)=j-i mod 12`, giving `12` translation orbits of size `12`; source-site and target-site labels fail translation invariance.  Thus `144` directed edges are `12` translated copies of each of `12` relative jumps, and the quotient renormalization factor is `1/12`.  This is finite quotient/renormalization progress only: no strict `Gamma_9_5` action-unit source, continuum field-variable provenance, nonproxy `L_total`, EOM, Hamiltonian, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2916/S1866 displacement-quotient `L_total` guard", "## P2916/S1866 displacement-quotient `L_total` guard\n\n`P2916/S1866` selects the displacement quotient `q(i,j)=j-i mod 12` by diagonal translation invariance and explains the `12` vs `144` mismatch as twelve translated copies of each relative jump.  This supplies a finite renormalization theorem, but it still cannot define nonproxy `L_total` until a strict nonzero `Gamma_9_5` source theorem or continuum field-variable provenance theorem couples this quotient to the variational integral.\n")
    append_once(AGENTS, "Current translation-orbit quotient renormalization guardrail (P2916/S1866, 2026-06-19)", "## Current translation-orbit quotient renormalization guardrail (P2916/S1866, 2026-06-19)\n\n- P2916 proves a finite quotient theorem selecting the displacement quotient `q(i,j)=j-i mod 12` under diagonal `Z12` translations; this explains the `12` vs `144` mismatch as `12` translated copies of each relative jump.\n- Treat this as finite quotient/renormalization progress only: no strict nonzero `Gamma_9_5` source theorem and no continuum field-variable provenance theorem are exported.\n- Do not promote the displacement quotient, P2911 pullback readiness, P2912 Jacobian readiness, P2913 Gamma source candidates, or P2914/P2915 measure evidence to nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure without those remaining theorems.\n- A next admissible move must prove exactly one of the remaining theorems: strict `Gamma_9_5` action-unit source coupled to the quotient integral, or continuum field-variable provenance for the displacement quotient.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
