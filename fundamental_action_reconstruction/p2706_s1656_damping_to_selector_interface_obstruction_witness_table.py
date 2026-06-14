#!/usr/bin/env python3
"""P2706/S1656: damping-to-selector interface obstruction/witness table.

Follow-up to P2705.  Test whether the P2377/P2378 damping-compression
transport primitive can become a non-premise Z12 directed-unit selector.
"""
from __future__ import annotations

import hashlib
import json
from itertools import combinations, product
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2706_s1656_damping_to_selector_interface_obstruction_witness_table.json"
MD = GEN / "p2706_s1656_damping_to_selector_interface_obstruction_witness_table.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2705": GEN / "p2705_s1655_release_9_3s_commit_boundary_alignment_audit.json",
    "P2377": GEN / "p2377_s1327_damping_compression_transport_primitive_uniform_coupling_theorem.json",
    "P2378": GEN / "p2378_s1328_unit_normalized_transport_coupling_insufficiency_theorem.json",
    "P2699": GEN / "p2699_s1649_z12_fractal_information_aut_invariant_selector_source_no_go.json",
    "P2700": GEN / "p2700_s1650_exhaustive_aut_invariant_selector_functional_enumeration_no_go.json",
    "P2701": GEN / "p2701_s1651_strict_sourced_symmetry_breaking_provider_inventory_matrix.json",
}

Z12_NODES = range(12)
SUPPORT_SIZE = 5
DIRECTED_UNITS = [1, 5, 7, 11]
SAMPLE_ETAS = [1.8, 1.9, 2.0]
SAMPLE_BETA_TORS = [0.0, 0.01, 0.1]
MASS_BUDGETS = [0.0, 1.0, 1.1757688202848233, 1.8435099395246706, 2.0]
STRICT_PARAMS = {"omega": 743.0 / 4000.0, "phi": 13.0 / 80.0, "beta": 1.0, "eta": 9.0 / 5.0}
NEGATIVE_EXPORT_FLAGS = [
    "qw2191_discharged",
    "non_premise_selector_provider_exported",
    "directed_unit_selected",
    "pair12_strict_core_upgrade_exported",
    "ltotal_promoted",
    "toe_closure_exported",
]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def k_strict(d: int) -> float:
    import math

    return math.cos(STRICT_PARAMS["omega"] * d + STRICT_PARAMS["phi"]) / (1.0 + STRICT_PARAMS["beta"] * d ** STRICT_PARAMS["eta"])


def compression_log_weight(d: int, eta: float, beta_tors: float) -> float:
    import math

    return math.log((1.0 + d**eta) / (1.0 + beta_tors * d))


def directed_edge_count(support: tuple[int, ...], step: int) -> int:
    support_set = set(support)
    return sum(1 for node in support_set if (node + step) % 12 in support_set)


def support_rows() -> list[dict[str, Any]]:
    rows = []
    for support in combinations(Z12_NODES, SUPPORT_SIZE):
        counts = {str(unit): directed_edge_count(support, unit) for unit in DIRECTED_UNITS}
        rows.append({"support": list(support), "directed_counts": counts})
    return rows


def orientation_blindness_scan() -> dict[str, Any]:
    rows = support_rows()
    worst_directed_count_delta = 0
    witness_examples = []
    for row in rows:
        counts = row["directed_counts"]
        for unit in [1, 5]:
            opposite = (-unit) % 12
            delta = counts[str(unit)] - counts[str(opposite)]
            worst_directed_count_delta = max(worst_directed_count_delta, abs(delta))
            if delta != 0 and len(witness_examples) < 5:
                witness_examples.append({"support": row["support"], "unit": unit, "opposite": opposite, "delta": delta})
    return {
        "support_count": len(rows),
        "directed_units": DIRECTED_UNITS,
        "worst_directed_count_delta_between_u_and_minus_u": worst_directed_count_delta,
        "nonzero_delta_witness_examples": witness_examples,
        "count_level_obstruction": "Every finite support has the same directed edge count for u and -u because edges are counted by membership pairs; damping weights multiply distances, not orientation.",
    }


def weighted_selector_scan() -> dict[str, Any]:
    rows = support_rows()
    max_abs_score_delta = 0.0
    total_cases = 0
    nonzero_cases = []
    samples = []
    for eta, beta_tors, mass in product(SAMPLE_ETAS, SAMPLE_BETA_TORS, MASS_BUDGETS):
        weights = {
            1: k_strict(1) + mass * compression_log_weight(1, eta, beta_tors),
            5: k_strict(5) + mass * compression_log_weight(5, eta, beta_tors),
        }
        for row in rows:
            counts = row["directed_counts"]
            plus_score = weights[1] * counts["1"] + weights[5] * counts["5"]
            minus_score = weights[1] * counts["11"] + weights[5] * counts["7"]
            delta = plus_score - minus_score
            total_cases += 1
            max_abs_score_delta = max(max_abs_score_delta, abs(delta))
            if abs(delta) > 1e-12 and len(nonzero_cases) < 5:
                nonzero_cases.append({"eta": eta, "beta_tors": beta_tors, "mass": mass, "support": row["support"], "delta": delta})
        if len(samples) < 9:
            samples.append({"eta": eta, "beta_tors": beta_tors, "mass": mass, "w1": weights[1], "w5": weights[5]})
    return {
        "total_weighted_cases": total_cases,
        "parameter_grid": {"etas": SAMPLE_ETAS, "beta_tors": SAMPLE_BETA_TORS, "mass_budgets": MASS_BUDGETS},
        "sample_weights": samples,
        "max_abs_score_delta_plus_direction_vs_minus_direction": max_abs_score_delta,
        "nonzero_score_delta_examples": nonzero_cases,
        "weighted_obstruction": "For every scanned eta/beta_tors/mass instance and every 5-node support, K(d)+M*C(d) scores +u and -u equally.  The damping primitive has distance weights only, so it cannot export a directed-unit selector.",
    }


def source_boundary() -> dict[str, Any]:
    p2705 = read_json(INPUTS["P2705"])
    p2377 = read_json(INPUTS["P2377"])
    p2378 = read_json(INPUTS["P2378"])
    p2699 = read_json(INPUTS["P2699"])
    p2700 = read_json(INPUTS["P2700"])
    p2701 = read_json(INPUTS["P2701"])
    return {
        "p2705_no_unlock": p2705.get("decision", {}).get("release_9_3s_pointer_unblocks_current_stage") is False,
        "p2377_status": p2377.get("status"),
        "p2378_status": p2378.get("status"),
        "p2699_bounded_no_go": p2699.get("decision", {}).get("bounded_no_go_now") is True,
        "p2700_bounded_no_go": p2700.get("decision", {}).get("bounded_no_go_now") is True,
        "p2701_bounded_no_go": p2701.get("decision", {}).get("bounded_no_go_now") is True,
    }


def obstruction_matrix(count_scan: dict[str, Any], weighted_scan: dict[str, Any], boundary: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "obligation": "directed_count_level_can_distinguish_plus_unit_from_minus_unit",
            "passes": count_scan["worst_directed_count_delta_between_u_and_minus_u"] == 0,
            "result": "obstruction",
            "evidence": {"support_count": count_scan["support_count"], "worst_delta": count_scan["worst_directed_count_delta_between_u_and_minus_u"]},
        },
        {
            "obligation": "p2377_p2378_weighted_scores_break_orientation_symmetry",
            "passes": weighted_scan["max_abs_score_delta_plus_direction_vs_minus_direction"] == 0.0,
            "result": "obstruction",
            "evidence": {"total_weighted_cases": weighted_scan["total_weighted_cases"], "max_abs_delta": weighted_scan["max_abs_score_delta_plus_direction_vs_minus_direction"]},
        },
        {
            "obligation": "prior_no_unlock_boundary_still_applies",
            "passes": boundary["p2705_no_unlock"] and boundary["p2699_bounded_no_go"] and boundary["p2700_bounded_no_go"] and boundary["p2701_bounded_no_go"],
            "result": "boundary_preserved",
            "evidence": boundary,
        },
    ]


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2706/S1656 damping-to-selector interface obstruction/witness table",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Obstruction matrix",
    ]
    for row in payload["obstruction_matrix"]:
        lines.append(f"- `{row['obligation']}`: passes={row['passes']}, result={row['result']}")
    lines.extend([
        "",
        "## Finite witness computation",
        f"- supports={payload['count_scan']['support_count']}",
        f"- weighted_cases={payload['weighted_scan']['total_weighted_cases']}",
        f"- max_abs_score_delta={payload['weighted_scan']['max_abs_score_delta_plus_direction_vs_minus_direction']}",
        "",
        "## Decision",
        payload["decision"]["current_unlock_reading"],
        "",
        "## Next honest step",
        payload["decision"]["next_honest_step"],
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    count_scan = orientation_blindness_scan()
    weighted_scan = weighted_selector_scan()
    boundary = source_boundary()
    matrix = obstruction_matrix(count_scan, weighted_scan, boundary)
    all_pass = all(row["passes"] for row in matrix)
    payload = {
        "status": "P2706_DAMPING_TO_SELECTOR_INTERFACE_OBSTRUCTION_NO_UNLOCK" if all_pass else "P2706_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "count_scan": count_scan,
        "weighted_scan": weighted_scan,
        "source_boundary": boundary,
        "obstruction_matrix": matrix,
        "decision": {
            "damping_transport_exports_directed_selector": False,
            "current_unlock_reading": "P2377/P2378 damping-compression transport is orientation-blind at the Z12 selector interface: across all 792 five-node supports and 45 eta/beta_tors/mass samples, +u and -u receive identical directed scores.  It therefore cannot discharge QW-2191 or export a non-premise directed-unit selector on current artifacts.",
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "next_honest_step": "No selector closure is unlocked by the damping interface.  The next admissible move requires a genuinely new strict-sourced symmetry-breaking provider or a new typed object outside the closed lanes; otherwise preserve the P2697-P2706 no-new-live-frontier/no-unlock certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2706/S1656 damping-to-selector interface obstruction", "## P2706/S1656 damping-to-selector interface obstruction\n\n`P2706/S1656` executes the P2705-recommended damping-to-selector witness table.  Across all 792 five-node Z12 supports and 45 sampled P2377/P2378 eta/beta_tors/mass settings, the distance-weighted damping score gives identical values to `+u` and `-u`.  Thus the damping-compression transport primitive is orientation-blind at the selector interface and does not export `QW-2191` discharge, a non-premise selector provider, pair12 strict-core upgrade, `L_total`, role transfer, bridge closure, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2706/S1656 damping interface Ltotal guard", "## P2706/S1656 damping interface Ltotal guard\n\n`P2706/S1656` proves a finite obstruction at the damping-to-selector interface: P2377/P2378 distance weights cannot distinguish `+u` from `-u` on the Z12 selector problem.  The obstruction is not a variational source and does not promote `L_total`, selector closure, pair12 strict-core upgrade, role transfer, bridge closure, or ToE closure.\n")
    append_once(AGENTS, "Current damping-to-selector interface obstruction guardrail (P2706/S1656, 2026-06-13)", "## Current damping-to-selector interface obstruction guardrail (P2706/S1656, 2026-06-13)\n\n- P2706 executes the P2705 damping-to-selector interface test and finds a finite orientation-blindness obstruction: all 792 five-node supports and 45 sampled P2377/P2378 parameter/mass settings give identical scores to `+u` and `-u`.\n- Do not use P2377/P2378 damping-compression transport to discharge `QW-2191`, export a non-premise selector provider, upgrade pair12 strict-core, promote `L_total`, transfer roles, close the bridge, or claim ToE closure.\n- A next admissible move requires a genuinely new strict-sourced symmetry-breaking provider or another new typed object outside closed lanes; otherwise preserve the P2697-P2706 no-new-live-frontier/no-unlock certificate.\n")
    return payload


if __name__ == "__main__":
    main()
