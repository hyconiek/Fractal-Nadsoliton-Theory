#!/usr/bin/env python3
"""P3012/S1962: strict-kernel phase-gradient selector-source obstruction.

P3011 separated the strict-kernel role from the selector role and allowed only a
new directed selector source or one closed kernel-to-coefficient atom.  The older
P2761-P2769 moment-provenance sequence already exhausts the generic
kernel-to-coefficient provenance replay, so this audit introduces exactly one new
candidate selector source: a strict-kernel phase/label-gradient arrow extracted
from sampled K_strict_gate values on Z/12Z unit labels.

The candidate is intentionally bounded.  It can produce a label-dependent arrow
(e.g. K(1) vs K(11)), but the finite Aut(Z12) equivariance test shows that this
arrow is not a nonpremise strict selector source: the score changes under unit
relabeling and therefore depends on an unexported chart/metric-label premise.
"""
from __future__ import annotations

import hashlib, json, math
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3011_s1961_toe_strict_kernel_selector_role_separation_matrix import OUT as P3011, strict_kernel

OUT = GEN / "p3012_s1962_strict_kernel_phase_gradient_selector_source_obstruction.json"
MD = GEN / "p3012_s1962_strict_kernel_phase_gradient_selector_source_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MODULUS = 12
UNITS = [1, 5, 7, 11]
PAIR = [1, 11]


def score(label: int) -> float:
    """Label-embedded strict-kernel score for one Z12 unit label."""
    return strict_kernel(label)


def phase_gradient_witness() -> dict[str, Any]:
    unit_scores = {str(u): score(u) for u in UNITS}
    pair_scores = {str(u): score(u) for u in PAIR}
    sorted_units = sorted(UNITS, key=lambda u: (-unit_scores[str(u)], u))
    pair_pick = max(PAIR, key=lambda u: pair_scores[str(u)])
    unit_pick = sorted_units[0]
    equivariance_rows = []
    for a, u in product(UNITS, UNITS):
        image = (a * u) % MODULUS
        equivariance_rows.append({
            "aut_unit": a,
            "source_unit": u,
            "image_unit": image,
            "source_score_rounded_12dp": round(unit_scores[str(u)], 12),
            "image_score_rounded_12dp": round(unit_scores[str(image)], 12),
            "score_preserved": abs(unit_scores[str(u)] - unit_scores[str(image)]) < 1e-12,
        })
    pair_rows = []
    for a, u in product(UNITS, PAIR):
        image = (a * u) % MODULUS
        pair_rows.append({
            "aut_unit": a,
            "source_directed_pair_unit": u,
            "image": image,
            "image_remains_in_pair": image in PAIR,
            "pair_score_available_after_action": image in PAIR,
        })
    invariant_directed_units = [u for u in UNITS if all((a * u) % MODULUS == u for a in UNITS)]
    return {
        "candidate_object": "StrictKernelPhaseGradientArrowSelectorSourceCandidate",
        "unit_scores_rounded_12dp": {k: round(v, 12) for k, v in unit_scores.items()},
        "pair_scores_rounded_12dp": {k: round(v, 12) for k, v in pair_scores.items()},
        "label_based_pair_pick": pair_pick,
        "label_based_unit_pick": unit_pick,
        "unit_score_order_desc": sorted_units,
        "pair_score_gap_K1_minus_K11": round(pair_scores["1"] - pair_scores["11"], 12),
        "unit_score_gap_top_minus_second": round(unit_scores[str(sorted_units[0])] - unit_scores[str(sorted_units[1])], 12),
        "equivariance_rows": equivariance_rows,
        "equivariance_row_count": len(equivariance_rows),
        "equivariance_failure_count": sum(1 for r in equivariance_rows if not r["score_preserved"]),
        "pair_action_rows": pair_rows,
        "pair_action_row_count": len(pair_rows),
        "pair_action_leaves_pair_count": sum(1 for r in pair_rows if not r["image_remains_in_pair"]),
        "aut_invariant_directed_units": invariant_directed_units,
        "accepted_as_nonpremise_selector_source": False,
    }


def obligation_rows(w: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"obligation": "finite_label_arrow_exists", "satisfied": w["pair_score_gap_K1_minus_K11"] != 0 and w["label_based_pair_pick"] == 1, "evidence": "K(1) and K(11) differ in the integer-label embedding"},
        {"obligation": "unit_orbit_score_table_exists", "satisfied": len(w["unit_scores_rounded_12dp"]) == 4, "evidence": "scores are computed for all Aut(Z12) unit labels"},
        {"obligation": "aut_equivariance", "satisfied": w["equivariance_failure_count"] == 0, "evidence": "a selector source must not depend on the arbitrary unit relabeling"},
        {"obligation": "pair_closed_under_full_aut_action", "satisfied": w["pair_action_leaves_pair_count"] == 0, "evidence": "the directed pair {1,11} is not preserved by all unit actions as a two-point carrier"},
        {"obligation": "nonpremise_chart_metric_source", "satisfied": False, "evidence": "no strict theorem exports the integer-label chart/metric that makes K(1)>K(11) canonical"},
        {"obligation": "accepted_current_selector_source", "satisfied": w["accepted_as_nonpremise_selector_source"], "evidence": "finite label arrow is not accepted because equivariance and chart-source obligations fail"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_label_arrow", "unit_score_table", "aut_equivariance", "pair_carrier_closed", "chart_metric_source", "torsor_coupling_theorem", "qw2191_boundary_respected"]
    rows = []
    for bits in product([False, True], repeat=len(names)):
        present = dict(zip(names, bits))
        rows.append({"present": present, "accepts_selector_source": all(bits)})
    return rows


def build_payload(p3011_path: Any) -> dict[str, Any]:
    witness = phase_gradient_witness()
    matrix = acceptance_matrix()
    return {
        "status": "P3012_STRICT_KERNEL_PHASE_GRADIENT_SELECTOR_SOURCE_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3011": hashlib.sha256(p3011_path.read_bytes()).hexdigest() if p3011_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "StrictKernelPhaseGradientSelectorSource_EquivarianceObstructionMatrix",
            "phase_gradient_witness": witness,
            "proof_obligation_rows": obligation_rows(witness),
            "finite_acceptance_matrix": matrix,
        },
        "selector_certificate": {
            "label_based_pair_pick": witness["label_based_pair_pick"],
            "label_based_unit_pick": witness["label_based_unit_pick"],
            "unit_score_order_desc": witness["unit_score_order_desc"],
            "pair_score_gap_K1_minus_K11": witness["pair_score_gap_K1_minus_K11"],
            "equivariance_row_count": witness["equivariance_row_count"],
            "equivariance_failure_count": witness["equivariance_failure_count"],
            "pair_action_leaves_pair_count": witness["pair_action_leaves_pair_count"],
            "aut_invariant_directed_unit_count": len(witness["aut_invariant_directed_units"]),
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_selector_source"]),
        },
        "decision": {
            "positive_progress": "P3012 constructs a new strict-kernel phase/label-gradient arrow candidate rather than replaying cube-map or generic selector lanes.",
            "breakthrough": "Bounded no-go: the label score can prefer +1 over -1, but the preference is not Aut(Z12)-equivariant and lacks a strict chart/metric source theorem, so it is not a nonpremise selector source.",
            "negative_export_flags": {k: False for k in ["strict_selector_source_exported", "qw2191_discharged", "directed_torsor_coupling_exported", "chart_metric_source_exported", "kernel_to_coefficient_provenance_exported", "unit_bearing_ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay kernel label-gradient selector attempts. A next proof-grade move must either supply a genuine strict chart/metric source theorem for a directed selector and rerun equivariance, or introduce a different new typed object outside cube-map, exhausted moment-provenance, selector replay, bridge, role-transfer, and L_total lanes.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["selector_certificate"]
    lines = [
        "# P3012/S1962 strict-kernel phase-gradient selector-source obstruction", "",
        f"Status: `{payload['status']}`", "", "## Selector certificate",
        f"- label-based pair pick: `{c['label_based_pair_pick']}`",
        f"- label-based unit pick: `{c['label_based_unit_pick']}`",
        f"- unit score order descending: `{c['unit_score_order_desc']}`",
        f"- K(1)-K(11): `{c['pair_score_gap_K1_minus_K11']}`",
        f"- Aut equivariance rows/failures: `{c['equivariance_row_count']}/{c['equivariance_failure_count']}`",
        f"- pair-action leaves-pair count: `{c['pair_action_leaves_pair_count']}`",
        f"- Aut-invariant directed unit count: `{c['aut_invariant_directed_unit_count']}`",
        f"- acceptance matrix rows/accepted: `{c['acceptance_matrix_rows']}/{c['accepted_rows']}`", "",
        "## Decision", payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P3011)
    payload = build_payload(P3011)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3012/S1962 strict-kernel phase-gradient selector-source obstruction", "## P3012/S1962 strict-kernel phase-gradient selector-source obstruction\n\n`P3012/S1962` constructs one new directed selector-source candidate after P3011: a strict-kernel phase/label-gradient arrow extracted from sampled `K_strict_gate` values on `Z/12Z` unit labels.  The finite score table can prefer `+1` over `-1` in the integer-label embedding (`K(1)-K(11)` is nonzero), but the preference is not `Aut(Z12)`-equivariant: the unit-score table has equivariance failures under unit relabeling, the pair `{+1,-1}` is not closed under all unit actions, and no strict chart/metric source theorem exports this label order as canonical.  Therefore no nonpremise selector source, `QW-2191` discharge, directed torsor coupling, kernel-to-coefficient provenance, unit-bearing `L_total`, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3012/S1962 phase-gradient selector `L_total` guard", "## P3012/S1962 phase-gradient selector `L_total` guard\n\n`P3012/S1962` adds no `L_total` term.  The strict-kernel label-gradient arrow is a finite selector-source candidate only; because it fails `Aut(Z12)` equivariance and lacks a strict chart/metric source theorem, it cannot install a selector source, directed torsor coupling, unit-bearing density, EOM/Hamiltonian term, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current strict-kernel phase-gradient selector-source obstruction guardrail (P3012/S1962, 2026-06-22)", "## Current strict-kernel phase-gradient selector-source obstruction guardrail (P3012/S1962, 2026-06-22)\n\n- P3012 introduces one new directed selector-source candidate after P3011: a strict-kernel phase/label-gradient arrow from sampled `K_strict_gate` values on `Z/12Z` unit labels.\n- The finite label score can prefer `+1` over `-1`, but this is not an `Aut(Z12)`-equivariant selector: unit relabeling changes scores, the pair `{+1,-1}` is not closed under all unit actions, and no strict chart/metric source theorem makes the label order canonical.\n- Do not promote strict-kernel phase-gradient arrows to nonpremise selector source, `QW-2191` discharge, directed torsor coupling, kernel-to-coefficient provenance, unit-bearing `L_total`, bridge closure, role transfer, observable generation, or ToE closure.\n- Do not replay kernel label-gradient selector attempts.  A next proof-grade move must supply a genuine strict chart/metric source theorem for a directed selector and rerun equivariance, or introduce a different new typed object outside cube-map, exhausted moment-provenance, selector replay, bridge, role-transfer, and `L_total` lanes.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
