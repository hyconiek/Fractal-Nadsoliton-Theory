#!/usr/bin/env python3
"""P2782/S1732: bipartite 16-node 4-regular enumerator scale obstruction.

P2781 recommended leaving hand-picked/structured families and moving toward a
canonical enumerator for connected 16-node 4-regular graphs.  This bounded
follow-up tests the smallest natural exact enumerator subproblem: labeled
bipartite 4-regular graphs with fixed parts 8+8, represented as 8x8 binary
matrices with every row and column sum equal to 4.

A dynamic-programming count gives 116,963,796,250 labeled matrices before
connectivity filtering and before quotienting by row/column permutations.  Thus a
naive in-repo canonical enumerator is not an honest next proof move without a
new canonical-generation theorem/tool; this is an enumerator-scale obstruction,
not a spectral-source theorem or geometry closure.
"""
from __future__ import annotations

import functools
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
P2781 = GEN / "p2781_s1731_enumerated_two_layer_c8_spectrum_collision_audit.json"
OUT = GEN / "p2782_s1732_bipartite_regular_enumerator_scale_obstruction.json"
MD = GEN / "p2782_s1732_bipartite_regular_enumerator_scale_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

LEFT_SIZE = 8
RIGHT_SIZE = 8
DEGREE = 4
EXPECTED_LABELED_BIPARTITE_COUNT = 116_963_796_250
NEGATIVE_EXPORT_FLAGS = [
    "canonical_16node_4regular_enumerator_exported",
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


def row_masks(width: int = RIGHT_SIZE, weight: int = DEGREE) -> tuple[int, ...]:
    return tuple(sum(1 << bit for bit in combo) for combo in itertools.combinations(range(width), weight))


def decrement_caps(caps: tuple[int, ...], mask: int) -> tuple[int, ...] | None:
    next_caps = list(caps)
    for bit in range(len(caps)):
        if mask & (1 << bit):
            if next_caps[bit] == 0:
                return None
            next_caps[bit] -= 1
    return tuple(next_caps)


def count_labeled_bipartite_regular_matrices(left_size: int = LEFT_SIZE, right_size: int = RIGHT_SIZE, degree: int = DEGREE) -> dict[str, Any]:
    masks = row_masks(right_size, degree)

    @functools.lru_cache(maxsize=None)
    def dp(row_index: int, caps: tuple[int, ...]) -> int:
        if row_index == left_size:
            return int(all(cap == 0 for cap in caps))
        # Quick feasibility pruning: remaining rows cannot fill columns that need too much.
        remaining_rows = left_size - row_index
        if any(cap < 0 or cap > remaining_rows for cap in caps):
            return 0
        total = 0
        for mask in masks:
            next_caps = decrement_caps(caps, mask)
            if next_caps is not None:
                total += dp(row_index + 1, next_caps)
        return total

    count = dp(0, tuple([degree] * right_size))
    return {
        "left_size": left_size,
        "right_size": right_size,
        "degree": degree,
        "row_mask_count": len(masks),
        "labeled_bipartite_regular_matrix_count": count,
        "dp_cache_states": dp.cache_info().currsize,
        "matches_expected_regression_value": count == EXPECTED_LABELED_BIPARTITE_COUNT,
    }


def scale_witness() -> dict[str, Any]:
    count_row = count_labeled_bipartite_regular_matrices()
    count = count_row["labeled_bipartite_regular_matrix_count"]
    row_perm = 1
    for value in range(2, LEFT_SIZE + 1):
        row_perm *= value
    col_perm = 1
    for value in range(2, RIGHT_SIZE + 1):
        col_perm *= value
    return {
        "audited_subproblem": "fixed-bipartition labeled 8x8 0/1 matrices with all row/column sums equal to 4",
        "count_row": count_row,
        "row_permutation_count": row_perm,
        "column_permutation_count": col_perm,
        "row_column_relabeling_group_size": row_perm * col_perm,
        "labeled_count_exceeds_one_billion": count > 1_000_000_000,
        "labeled_count_exceeds_one_hundred_billion": count > 100_000_000_000,
        "connectivity_filter_not_yet_applied": True,
        "graph_isomorphism_quotient_not_yet_applied": True,
        "naive_enumerator_blocked_in_repo": count > 100_000_000_000,
        "finite_obstruction_statement": "Even the fixed-bipartition 16-node 4-regular subproblem has 116,963,796,250 labeled candidates before connectivity and isomorphism quotienting.",
    }


def acceptance_matrix(witness: dict[str, Any], p2781: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2781_local_family_witness_present": p2781.get("status") == "P2781_ENUMERATED_TWO_LAYER_C8_SPECTRUM_COLLISION_AUDIT_NO_CLOSURE",
        "exact_dp_count_performed": witness["count_row"]["matches_expected_regression_value"],
        "fixed_bipartition_subproblem_exceeds_naive_in_repo_scale": witness["naive_enumerator_blocked_in_repo"],
        "canonical_generation_algorithm_or_external_certificate_supplied": False,
        "connectivity_and_isomorphism_quotient_completed": False,
        "strict_nadsoliton_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_enumerator_scale_obstruction": facts["exact_dp_count_performed"] and facts["fixed_bipartition_subproblem_exceeds_naive_in_repo_scale"],
        "accepted_as_canonical_16node_4regular_enumerator": False,
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "The exact fixed-bipartition count is too large for a naive in-repo canonical enumerator, and no canonical-generation theorem/tool certificate is supplied.  No spectral source law or K/L_total coupling is exported.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    witness = payload["bipartite_regular_enumerator_scale_witness"]
    count_row = witness["count_row"]
    lines = [
        "# P2782/S1732 bipartite regular enumerator scale obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Exact DP count",
        f"- row_mask_count={count_row['row_mask_count']}",
        f"- dp_cache_states={count_row['dp_cache_states']}",
        f"- labeled_bipartite_regular_matrix_count={count_row['labeled_bipartite_regular_matrix_count']}",
        f"- row_column_relabeling_group_size={witness['row_column_relabeling_group_size']}",
        f"- naive_enumerator_blocked_in_repo={witness['naive_enumerator_blocked_in_repo']}",
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
    witness = scale_witness()
    acceptance = acceptance_matrix(witness, p2781)
    payload = {
        "status": "P2782_BIPARTITE_REGULAR_ENUMERATOR_SCALE_OBSTRUCTION_NO_CLOSURE",
        "input_hashes": {"P2781": sha(P2781)},
        "input_statuses": {"P2781": p2781.get("status")},
        "audited_question": "Is a naive canonical enumerator for 16-node 4-regular graphs honest after P2781?",
        "bipartite_regular_enumerator_scale_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not start a naive full 16-node 4-regular enumeration inside the repo.  The next honest move is exactly one of: supply a canonical-generation theorem/tool certificate (for example a checked nauty/geng-style graph6 import with reproducible hashes) and then run the full-spectrum collision audit; or export a strict nadsoliton spectral action/source law that fixes the admissible class and target before testing.  Otherwise preserve the P2697-P2782 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2782/S1732 bipartite regular enumerator scale obstruction", "## P2782/S1732 bipartite regular enumerator scale obstruction\n\n`P2782/S1732` audits the feasibility of moving from structured 16-node families to a naive canonical enumerator.  An exact dynamic-programming count for the fixed-bipartition bipartite subproblem—8x8 binary matrices with all row and column sums equal to 4—finds `116,963,796,250` labeled candidates before connectivity filtering and graph-isomorphism quotienting.  This blocks a naive in-repo full enumerator without a canonical-generation theorem/tool certificate.  No strict nadsoliton spectral source law, global graph-class theorem, kernel full-expression theorem, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2782/S1732 enumerator-scale Ltotal guard", "## P2782/S1732 enumerator-scale Ltotal guard\n\n`P2782/S1732` adds no variational source term.  The exact fixed-bipartition enumerator count is a scale obstruction to naive exhaustive graph search, not a sourced `K`/`L_total` spectral action; therefore it cannot promote role-bearing `L_total` or canonical nadsoliton geometry.\n")
    append_once(AGENTS, "Current bipartite regular enumerator scale obstruction guardrail (P2782/S1732, 2026-06-15)", "## Current bipartite regular enumerator scale obstruction guardrail (P2782/S1732, 2026-06-15)\n\n- P2782 tests the feasibility of the P2781-recommended canonical 16-node 4-regular enumeration by exactly counting the fixed-bipartition bipartite 8+8, degree-4 subproblem.\n- The dynamic-programming count gives `116,963,796,250` labeled 8x8 row/column-sum-4 matrices before connectivity filtering and graph-isomorphism quotienting, blocking a naive in-repo enumerator.\n- Do not promote enumerator scale evidence to canonical geometry, strict spectral source law, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  A next admissible move must supply a canonical-generation theorem/tool certificate for the full graph class, or export a strict spectral action/source law before testing.\n")
    return payload


if __name__ == "__main__":
    main()
