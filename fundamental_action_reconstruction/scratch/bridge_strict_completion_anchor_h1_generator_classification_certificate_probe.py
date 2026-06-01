#!/usr/bin/env python3
"""Scratch probe: anchor vs H^1 generator classification certificate.

This probe audits the claim that the remaining selector/anchor bit is already
"the generator" of H^1(Z12;GF(2)).  The previous path and cycle reports prove
that the open path has H^1=0, while the artificially closed 12-cycle has a
one-dimensional GF(2) cohomology quotient detected by total cycle parity.

The finite linear-algebra distinction is important:

* the left anchor b(0)=0 is a C^0 gauge fixing for an exact path cochain;
* the closed-cycle H^1 generator is a non-exact C^1 edge-cochain class with odd
  cycle parity, represented here by the single closing edge 11->0;
* the audited closed cochain has even parity and therefore represents the zero
  H^1 class.

So the AI opinion is accepted only in its weak form (the closed 12-cycle has a
Z2 obstruction/generator and the source remains open), but rejected as stated if
it identifies the d=0 anchor itself with the H^1 generator.  This is finite
GF(2) bookkeeping only; it is not a selector-source theorem, not beta_tors ->
chi_11, and not a bridge theorem.
"""
from __future__ import annotations

import json
from itertools import product
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_anchor_h1_generator_classification_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_anchor_h1_generator_classification_certificate_report.md"
PATH_REPORT = HERE / "bridge_strict_completion_phase_sign_path_cohomology_triviality_certificate_report.json"
CYCLE_REPORT = HERE / "bridge_strict_completion_phase_sign_cycle_closure_obstruction_certificate_report.json"
Z2_REPORT = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"
GF2_REPORT = HERE / "bridge_strict_completion_phase_sign_gf2_linear_system_certificate_report.json"

NODE_COUNT = 12
PATH_EDGE_COUNT = 11
CYCLE_EDGE_COUNT = 12
CLOSING_EDGE = "11->0"
CYCLE_EDGE_LABELS = [f"{d}->{d + 1}" for d in range(PATH_EDGE_COUNT)] + [CLOSING_EDGE]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def cycle_coboundary_matrix() -> list[list[int]]:
    rows = []
    for edge_index in range(CYCLE_EDGE_COUNT):
        row = [0] * NODE_COUNT
        row[edge_index % NODE_COUNT] = 1
        row[(edge_index + 1) % NODE_COUNT] = 1
        rows.append(row)
    return rows


def matvec_gf2(matrix: list[list[int]], vector: list[int]) -> list[int]:
    return [sum(value & vector[col] for col, value in enumerate(row)) % 2 for row in matrix]


def gf2_rref(matrix: list[list[int]]) -> dict[str, Any]:
    work = [row[:] for row in matrix]
    row_count = len(work)
    col_count = len(work[0]) if work else 0
    pivot_row = 0
    pivots = []
    for col in range(col_count):
        pivot = None
        for candidate in range(pivot_row, row_count):
            if work[candidate][col]:
                pivot = candidate
                break
        if pivot is None:
            continue
        if pivot != pivot_row:
            work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        for target in range(row_count):
            if target != pivot_row and work[target][col]:
                work[target] = [a ^ b for a, b in zip(work[target], work[pivot_row])]
        pivots.append({"row": pivot_row, "pivot_col": col})
        pivot_row += 1
        if pivot_row == row_count:
            break
    return {
        "rank": len(pivots),
        "nullity": col_count - len(pivots),
        "pivot_rows": pivots,
        "rref_matrix": work,
    }


def image_vectors(matrix: list[list[int]]) -> set[tuple[int, ...]]:
    return {tuple(matvec_gf2(matrix, list(candidate))) for candidate in product([0, 1], repeat=NODE_COUNT)}


def parity(vector: list[int]) -> int:
    return sum(vector) % 2


def build_payload() -> dict[str, Any]:
    path = load_json(PATH_REPORT)
    cycle = load_json(CYCLE_REPORT)
    z2 = load_json(Z2_REPORT)
    gf2 = load_json(GF2_REPORT)

    node_bits = [row["node_bit"] for row in z2["node_bit_rows"]]
    path_edge_bits = [row["edge_bit"] for row in z2["edge_bit_rows"]]
    forced_closing_bit = node_bits[-1] ^ node_bits[0]
    audited_cycle_bits = path_edge_bits + [forced_closing_bit]
    generator_bits = [0] * PATH_EDGE_COUNT + [1]
    alternate_generator_bits = [1] + [0] * (CYCLE_EDGE_COUNT - 1)
    even_pair_bits = [1, 1] + [0] * (CYCLE_EDGE_COUNT - 2)

    matrix = cycle_coboundary_matrix()
    rref = gf2_rref(matrix)
    image = image_vectors(matrix)
    image_size = len(image)

    generator_in_image = tuple(generator_bits) in image
    alternate_generator_same_class = tuple([a ^ b for a, b in zip(generator_bits, alternate_generator_bits)]) in image
    even_pair_in_image = tuple(even_pair_bits) in image
    audited_in_image = tuple(audited_cycle_bits) in image

    classification_rows = [
        {
            "claim": "H^1(Z12;GF(2)) is one-dimensional on the closed 12-cycle.",
            "verdict": "accepted_for_closed_cycle_only",
            "finite_witness": "cycle report has rank(delta)=11, dim C1=12, so dim H^1=1",
        },
        {
            "claim": "The d=0 left anchor is itself the H^1 generator.",
            "verdict": "rejected_type_mismatch",
            "finite_witness": "b(0) is a C0 gauge-fixing value; an H^1 generator is a C1/im(delta) edge-cochain class with odd cycle parity",
        },
        {
            "claim": "The audited phase pattern represents the nonzero H^1 generator.",
            "verdict": "rejected_for_audited_zero_closure",
            "finite_witness": "audited path plus forced closing bit has total parity 0 and is exact on the closed cycle",
        },
        {
            "claim": "The missing selector/source is still open.",
            "verdict": "accepted",
            "finite_witness": "the generator classification gives the address of the obstruction, not a source theorem or QW-2191 discharge",
        },
    ]

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_ANCHOR_H1_GENERATOR_CLASSIFICATION_CERTIFICATE__GF2_TYPE_AUDIT",
        "status": "anchor-is-c0-gauge-fix-not-h1-generator__closed-cycle-generator-is-odd-parity-c1-class",
        "source_reports": {
            "phase_sign_path_cohomology_triviality_certificate": str(PATH_REPORT.relative_to(ROOT)),
            "phase_sign_cycle_closure_obstruction_certificate": str(CYCLE_REPORT.relative_to(ROOT)),
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
            "phase_sign_gf2_linear_system_certificate": str(GF2_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "H^1",
                "cohomology generator",
                "anchor",
                "left anchor",
                "cycle closure",
                "Z12",
                "selector source",
            ],
            "conclusion": "Existing reports certify path H1=0 and cycle H1=1, but no report explicitly separated the C0 anchor from the closed-cycle C1/im(delta) generator before this file.",
        },
        "cohomology_type_ledger": {
            "field": "GF(2)",
            "path_dim_C0": NODE_COUNT,
            "path_dim_C1": PATH_EDGE_COUNT,
            "path_h1_dimension": path["path_cohomology_summary"]["h1_dimension_dim_C1_minus_rank_delta"],
            "path_anchor_object": "b(0)=0 in C0 / node-potential gauge fixing",
            "cycle_dim_C0": NODE_COUNT,
            "cycle_dim_C1": CYCLE_EDGE_COUNT,
            "cycle_rank_delta": rref["rank"],
            "cycle_h1_dimension_dim_C1_minus_rank_delta": CYCLE_EDGE_COUNT - rref["rank"],
            "cycle_generator_object": "odd-parity C1 edge-cochain class modulo im(delta)",
            "cycle_parity_functional": "sum of all 12 cycle-edge bits mod 2",
            "image_size": image_size,
            "expected_image_size_2_rank": 2 ** rref["rank"],
        },
        "representative_rows": [
            {
                "name": "audited_closed_phase_cochain",
                "edge_bits": audited_cycle_bits,
                "parity": parity(audited_cycle_bits),
                "in_image_delta": audited_in_image,
                "h1_class": "zero",
            },
            {
                "name": "closing_edge_generator_candidate",
                "edge_bits": generator_bits,
                "support": [CLOSING_EDGE],
                "parity": parity(generator_bits),
                "in_image_delta": generator_in_image,
                "h1_class": "nonzero_generator",
            },
            {
                "name": "first_edge_generator_candidate",
                "edge_bits": alternate_generator_bits,
                "support": ["0->1"],
                "parity": parity(alternate_generator_bits),
                "same_h1_class_as_closing_edge_generator": alternate_generator_same_class,
                "h1_class": "nonzero_generator",
            },
            {
                "name": "even_two_edge_boundary_candidate",
                "edge_bits": even_pair_bits,
                "support": ["0->1", "1->2"],
                "parity": parity(even_pair_bits),
                "in_image_delta": even_pair_in_image,
                "h1_class": "zero",
            },
        ],
        "opinion_audit_rows": classification_rows,
        "classification_summary": {
            "path_h1_zero": path["path_cohomology_summary"]["h1_dimension_dim_C1_minus_rank_delta"] == 0,
            "cycle_h1_one": CYCLE_EDGE_COUNT - rref["rank"] == 1,
            "cycle_image_equals_even_parity_kernel": image_size == 2 ** rref["rank"] and all(parity(list(vector)) == 0 for vector in image),
            "closing_edge_odd_parity_generator": parity(generator_bits) == 1 and not generator_in_image,
            "alternate_odd_edge_same_generator_class": alternate_generator_same_class,
            "audited_closed_cochain_is_exact_zero_h1_class": audited_in_image and parity(audited_cycle_bits) == 0,
            "left_anchor_is_c0_gauge_fix_not_c1_generator": path["cochain_complex_definition"]["left_anchor_node_bit_b0"] == 0,
            "ai_opinion_as_stated_rejected": True,
            "ai_opinion_weak_form_accepted": True,
            "selector_source_remains_open": True,
            "gf2_path_solution_unique_inherited": gf2["gf2_linear_system_summary"]["unique_solution"],
        },
        "proof_certificate": {
            "path_step": "On the open path, rank(delta)=11=dim C1, so H^1(path;GF(2))=0; b(0)=0 only fixes the constant C0 kernel.",
            "cycle_step": "On the closed 12-cycle, rank(delta)=11 and dim C1=12, so C1/im(delta) has dimension 1.",
            "parity_step": "Every coboundary has even total cycle parity; the image has 2^11 vectors, matching the even-parity kernel, so odd parity represents the nonzero H^1 class.",
            "generator_step": "The single closing-edge cochain 11->0 has odd parity, is not a coboundary, and is therefore a valid generator representative; any other single edge differs from it by an exact even-parity cochain.",
            "anchor_step": "The left anchor b(0)=0 is a node-potential gauge choice in C0, not a C1 edge cochain and not an H^1 quotient class.",
            "audited_pattern_step": "The audited path plus forced closing bit has parity 0 and is exact, so it represents the zero H^1 class rather than the generator.",
            "theoretical_limit": "This locates the type and representative of the closed-cycle obstruction only; it does not derive the selector/source bit from nadsoliton dynamics, beta_tors, SSB, observer readout, or any bridge theorem.",
        },
        "hard_limits": [
            "No beta_tors -> chi_11 theorem is claimed.",
            "No strict-core selector source is exported.",
            "No QW-2191 selector discharge is claimed.",
            "No legacy physical-role transfer to K_strict_gate is claimed.",
            "No unqualified identity K_legacy_ont == K_strict_gate is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    summary = payload["classification_summary"]
    lines = [
        "# Anchor vs H1 generator classification certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "The d=0 left anchor is classified as a C0 gauge-fixing datum, not as the",
        "closed-cycle H1 generator.  The closed-cycle generator is represented by an",
        "odd-parity C1 edge cochain such as the single closing edge `11->0`.",
        "",
        "## Summary",
        "",
    ]
    for key, value in summary.items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Opinion audit", ""])
    for row in payload["opinion_audit_rows"]:
        lines.append(f"- `{row['verdict']}`: {row['claim']} ({row['finite_witness']})")
    lines.extend(["", "## Hard limits", ""])
    for limit in payload["hard_limits"]:
        lines.append(f"- {limit}")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(payload)
    print(json.dumps(payload["classification_summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
