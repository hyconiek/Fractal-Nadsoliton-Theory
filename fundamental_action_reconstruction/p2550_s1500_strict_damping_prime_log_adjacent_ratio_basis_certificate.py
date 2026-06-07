#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from typing import Any

import sympy as sp

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2550_s1500_strict_damping_prime_log_adjacent_ratio_basis_certificate.json"
MD = GEN / "p2550_s1500_strict_damping_prime_log_adjacent_ratio_basis_certificate.md"

SOURCE_FILES = {
    "P2527_PRIME_LOG_SLOPE_LINE": GEN / "p2527_s1477_strict_damping_prime_log_proportionality_slope_line_certificate.json",
    "P2542_PRIME_LOG_OBSTRUCTION": GEN / "p2542_s1492_strict_damping_prime_log_proportionality_current_premise_obstruction_certificate.json",
    "P2547_POST_IDENTITY_RESIDUAL_TRIKEY": GEN / "p2547_s1497_strict_damping_post_identity_residual_trikey_certificate.json",
    "P2549_TRACE_SOURCE_OBSTRUCTION": GEN / "p2549_s1499_strict_damping_quadruple_trace_current_premise_obstruction_certificate.json",
}

PRIMES = [2, 3, 5, 7, 11]
RATIO_SYMBOLS = [f"r_{p}" for p in PRIMES]
ADJACENT_EDGES = [(0, 1), (1, 2), (2, 3), (3, 4)]
NEGATIVE_EXPORT_FLAGS = [
    "prime_log_proportionality_source_exported",
    "slope_value_or_prime_anchor_source_exported",
    "m2_operator_signature_source_exported",
    "strict_quadruple_trace_source_exported",
    "beta_eta_numeric_source_exported",
    "strict_damping_beta_eta_source_exported",
    "source_obligation_discharge_exported",
    "damping_compression_bridge_component_ready",
    "full_bridge_theorem_exported",
    "role_transfer_theorem_exported",
    "selector_closure_exported",
    "qw2191_discharged_by_this_certificate",
    "role_bearing_ltotal_exported",
    "toe_closure_claimed",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2550|S1500|adjacent ratio basis|prime-log adjacent|ratio-basis certificate",
        "intended_research_nonduplication": "adjacent.*prime.*ratio|prime.*adjacent.*ratio|ratio equality basis|minimal ratio basis|prime-log.*basis|ratio-basis",
        "prime_log_precursors": "P2527|S1477|P2542|S1492|prime_log_proportionality_source|prime-log proportionality|normalized ratio",
        "post_identity_frontier": "P2547|S1497|P2549|S1499|residual.*prime-log|quadruple trace.*obstruction",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def load_theorems(sources: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {
        "P2527_PRIME_LOG_SLOPE_LINE": theorem(sources["P2527_PRIME_LOG_SLOPE_LINE"], "strict_damping_prime_log_proportionality_slope_line_certificate"),
        "P2542_PRIME_LOG_OBSTRUCTION": theorem(sources["P2542_PRIME_LOG_OBSTRUCTION"], "strict_damping_prime_log_proportionality_current_premise_obstruction_certificate"),
        "P2547_POST_IDENTITY_RESIDUAL_TRIKEY": theorem(sources["P2547_POST_IDENTITY_RESIDUAL_TRIKEY"], "strict_damping_post_identity_residual_trikey_certificate"),
        "P2549_TRACE_SOURCE_OBSTRUCTION": theorem(sources["P2549_TRACE_SOURCE_OBSTRUCTION"], "strict_damping_quadruple_trace_current_premise_obstruction_certificate"),
    }


def adjacent_constraint_matrix() -> sp.Matrix:
    rows = []
    for left, right in ADJACENT_EDGES:
        row = [0] * len(PRIMES)
        row[left] = 1
        row[right] = -1
        rows.append(row)
    return sp.Matrix(rows)


def matrix_rows_as_ints(matrix: sp.Matrix) -> list[list[int]]:
    return [[int(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def nullspace_as_ints(nullspace: list[sp.Matrix]) -> list[list[int]]:
    return [[int(value) for value in vector] for vector in nullspace]


def omission_witness(edge_index: int) -> dict[str, Any]:
    omitted = ADJACENT_EDGES[edge_index]
    ratios = [0 if i <= edge_index else 1 for i in range(len(PRIMES))]
    satisfied_edges = []
    failed_edges = []
    for current_index, (left, right) in enumerate(ADJACENT_EDGES):
        defect = ratios[left] - ratios[right]
        row = {
            "edge_index": current_index,
            "edge_primes": [PRIMES[left], PRIMES[right]],
            "defect": defect,
            "satisfied": defect == 0,
        }
        if defect == 0:
            satisfied_edges.append(row)
        else:
            failed_edges.append(row)
    prime_values_symbolic = [f"{ratio}*log({prime})" for ratio, prime in zip(ratios, PRIMES)]
    prime_values_numeric = [ratio * math.log(prime) for ratio, prime in zip(ratios, PRIMES)]
    return {
        "omitted_edge_index": edge_index,
        "omitted_edge_primes": [PRIMES[omitted[0]], PRIMES[omitted[1]]],
        "normalized_ratios_r_p": dict(zip(RATIO_SYMBOLS, ratios)),
        "prime_values_v_p_symbolic": dict(zip([f"v_{p}" for p in PRIMES], prime_values_symbolic)),
        "prime_values_v_p_numeric": dict(zip([f"v_{p}" for p in PRIMES], prime_values_numeric)),
        "satisfies_all_non_omitted_adjacent_equalities": len(satisfied_edges) == len(ADJACENT_EDGES) - 1,
        "violates_omitted_adjacent_equality": len(failed_edges) == 1 and failed_edges[0]["edge_index"] == edge_index,
        "ratio_spread": max(ratios) - min(ratios),
        "prime_log_proportionality_accepts": max(ratios) == min(ratios),
        "satisfied_edges": satisfied_edges,
        "failed_edges": failed_edges,
    }


def proper_subset_audit() -> dict[str, Any]:
    witnesses = [omission_witness(edge_index) for edge_index in range(len(ADJACENT_EDGES))]
    return {
        "single_omission_witnesses": witnesses,
        "single_omission_witness_count": len(witnesses),
        "all_single_omissions_have_nonproportional_witness": all(not witness["prime_log_proportionality_accepts"] for witness in witnesses),
        "all_single_omissions_satisfy_remaining_edges": all(witness["satisfies_all_non_omitted_adjacent_equalities"] for witness in witnesses),
        "basis_irredundant_under_single_omission": all(
            witness["satisfies_all_non_omitted_adjacent_equalities"] and not witness["prime_log_proportionality_accepts"]
            for witness in witnesses
        ),
    }


def basis_audit() -> dict[str, Any]:
    matrix = adjacent_constraint_matrix()
    nullspace = matrix.nullspace()
    omission = proper_subset_audit()
    return {
        "prime_order": PRIMES,
        "ratio_symbols": RATIO_SYMBOLS,
        "adjacent_edges_by_prime": [[PRIMES[left], PRIMES[right]] for left, right in ADJACENT_EDGES],
        "constraint_matrix_rows": matrix_rows_as_ints(matrix),
        "constraint_matrix_rank": int(matrix.rank()),
        "constraint_matrix_nullity": len(PRIMES) - int(matrix.rank()),
        "nullspace_basis": nullspace_as_ints(nullspace),
        "nullspace_is_constant_ratio_line": nullspace_as_ints(nullspace) == [[1, 1, 1, 1, 1]],
        "full_adjacent_basis_equivalent_to_prime_log_proportionality": int(matrix.rank()) == 4 and nullspace_as_ints(nullspace) == [[1, 1, 1, 1, 1]],
        "proper_subset_audit": omission,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    source_theorems = load_theorems(sources)
    basis = basis_audit()
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2550_T1_strict_damping_prime_log_adjacent_ratio_basis_certificate",
        "audited_chain": ["P2527/S1477", "P2542/S1492", "P2547/S1497", "P2549/S1499"],
        "frontier_source_under_attack": "prime_log_proportionality_source",
        "p2527_prime_log_slope_line_inherited": bool(source_theorems["P2527_PRIME_LOG_SLOPE_LINE"]),
        "p2542_prime_log_current_premise_obstruction_inherited": source_theorems["P2542_PRIME_LOG_OBSTRUCTION"].get("current_unital_multiplicative_premises_do_not_entail_prime_log_proportionality") is True,
        "p2547_post_identity_residual_trikey_inherited": source_theorems["P2547_POST_IDENTITY_RESIDUAL_TRIKEY"].get("post_identity_residual_trikey_irredundancy_exported") is True,
        "p2549_trace_source_obstruction_inherited": source_theorems["P2549_TRACE_SOURCE_OBSTRUCTION"].get("current_premise_quadruple_trace_nonentailment_exported") is True,
        "prime_log_adjacent_ratio_basis_audit": basis,
        "minimal_constraint_count_for_five_prime_ratios": len(ADJACENT_EDGES),
        "full_adjacent_ratio_basis_suffices_for_prime_log_proportionality": basis["full_adjacent_basis_equivalent_to_prime_log_proportionality"],
        "single_omission_countermodels_exported": basis["proper_subset_audit"]["basis_irredundant_under_single_omission"],
        "adjacent_ratio_basis_is_source_obligation_not_source": True,
        "prime_log_proportionality_source_exported": False,
        "residual_after_hypothetical_identity_m2_and_prime_log": ["slope_value_or_prime_anchor_source"],
        "recommended_next_honest_step": (
            "Do not treat the adjacent ratio basis as a strict source: it is only the exact four-constraint obligation needed to collapse five prime ratios. "
            "The next honest step is to seek a strict nadsoliton mechanism that exports all four adjacent prime-ratio equalities at once; "
            "if such a mechanism is unavailable, pivot to the remaining slope_value_or_prime_anchor_source with the prime-log basis kept conditional."
        ),
        "not_licensed": [
            "No strict mechanism exporting the four adjacent prime-ratio equalities is supplied.",
            "No prime-log source, slope-value source, m2 source, or strict damping beta/eta source is exported.",
            "No bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
    }
    for flag in NEGATIVE_EXPORT_FLAGS:
        theorem_export[flag] = False
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2542_obstruction_inherited": theorem_export["p2542_prime_log_current_premise_obstruction_inherited"],
        "p2547_residual_frontier_inherited": theorem_export["p2547_post_identity_residual_trikey_inherited"],
        "p2549_trace_obstruction_inherited": theorem_export["p2549_trace_source_obstruction_inherited"],
        "adjacent_matrix_rank_four": basis["constraint_matrix_rank"] == 4,
        "adjacent_matrix_nullity_one": basis["constraint_matrix_nullity"] == 1,
        "single_omission_countermodels_exist": theorem_export["single_omission_countermodels_exported"],
        "no_false_prime_log_source_export": theorem_export["prime_log_proportionality_source_exported"] is False,
        "no_false_bridge_or_role_transfer": theorem_export["full_bridge_theorem_exported"] is False and theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_or_toe_claim": theorem_export["qw2191_discharged_by_this_certificate"] is False and theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2550",
        "stage_id": "S1500",
        "status": "STRICT_DAMPING_PRIME_LOG_ADJACENT_RATIO_BASIS_CERTIFICATE_SOURCE_OBLIGATION_ONLY_NO_PRIME_LOG_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_prime_log_adjacent_ratio_basis_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_prime_log_adjacent_ratio_basis_certificate"]["theorem_export"]
    basis = t["prime_log_adjacent_ratio_basis_audit"]
    lines = [
        "# P2550/S1500 strict damping prime-log adjacent ratio basis certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier source under attack: `{t['frontier_source_under_attack']}`.",
        f"- P2542 prime-log obstruction inherited: `{t['p2542_prime_log_current_premise_obstruction_inherited']}`.",
        f"- P2547 residual tri-key inherited: `{t['p2547_post_identity_residual_trikey_inherited']}`.",
        f"- P2549 trace-source obstruction inherited: `{t['p2549_trace_source_obstruction_inherited']}`.",
        f"- Prime order: `{basis['prime_order']}`.",
        f"- Adjacent ratio edges: `{basis['adjacent_edges_by_prime']}`.",
        f"- Constraint matrix rank/nullity: `{basis['constraint_matrix_rank']}/{basis['constraint_matrix_nullity']}`.",
        f"- Nullspace basis: `{basis['nullspace_basis']}`.",
        f"- Full adjacent basis equivalent to prime-log proportionality: `{basis['full_adjacent_basis_equivalent_to_prime_log_proportionality']}`.",
        f"- Single-omission countermodels exported: `{t['single_omission_countermodels_exported']}`.",
        f"- Prime-log source exported: `{t['prime_log_proportionality_source_exported']}`.",
        "",
        "## Interpretation",
        "",
        "The four adjacent equalities `r_2=r_3`, `r_3=r_5`, `r_5=r_7`, and `r_7=r_11` are an exact finite basis for collapsing the five normalized prime ratios `r_p=v_p/log(p)` to one line.  Each single omitted equality has an explicit nonproportional witness satisfying all remaining adjacent equalities, so the four constraints are irredundant.",
        "",
        "This is a source-obligation basis, not a strict source theorem: P2550 does not explain why strict nadsoliton dynamics must export those four ratio equalities.",
        "",
        "## Recommended next honest step",
        "",
        t["recommended_next_honest_step"],
        "",
        "## Negative controls",
        "",
        "No prime-log source, residual slope source, `m2_operator_signature_source`, exact trace source, source-obligation discharge, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, physical-value generator, or ToE closure is exported.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_prime_log_adjacent_ratio_basis_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2550/S1500` attacks the residual `prime_log_proportionality_source` by isolating its exact finite source-obligation basis.  For normalized prime ratios `r_p=v_p/log(p)` on primes `2,3,5,7,11`, the four adjacent equations `r_2=r_3`, `r_3=r_5`, `r_5=r_7`, and `r_7=r_11` have rank/nullity `4/1` and nullspace the constant ratio line.  Every single omitted adjacent equality has an explicit nonproportional witness satisfying all other adjacent equalities, so the basis is irredundant.  This still exports no prime-log source: a strict nadsoliton mechanism must supply all four ratio equalities rather than merely list them.
""".strip()
    lag_section = """
`P2550/S1500` keeps the prime-log part of damping non-role-bearing in `L_total`.  The adjacent-ratio audit gives an exact four-constraint obligation for prime-log proportionality, but it does not derive those constraints from strict dynamics.  Therefore the damping/compression term still lacks prime-log source export, slope/anchor source export, bridge completion, and role-transfer license.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2550/S1500 prime-log adjacent-ratio basis guard", "## P2550/S1500 prime-log adjacent-ratio basis guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2550/S1500 prime-log adjacent-ratio Ltotal guard", "## P2550/S1500 prime-log adjacent-ratio Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
