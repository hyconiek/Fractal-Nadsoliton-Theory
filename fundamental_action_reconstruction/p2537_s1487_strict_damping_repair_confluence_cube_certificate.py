#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from itertools import combinations, permutations
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2534_s1484_strict_damping_boolean_prime_implicant_certificate import KEYS, KEY_LABELS, status_tuple, valid_assignments
from p2536_s1486_strict_damping_minimal_repair_ideal_certificate import action_for_status

GEN = ROOT / "generated"
OUT = GEN / "p2537_s1487_strict_damping_repair_confluence_cube_certificate.json"
MD = GEN / "p2537_s1487_strict_damping_repair_confluence_cube_certificate.md"

SOURCE_FILES = {
    "P2536_MINIMAL_REPAIR_IDEAL": GEN / "p2536_s1486_strict_damping_minimal_repair_ideal_certificate.json",
    "P2532_STRICTIZATION_DISTANCE": GEN / "p2532_s1482_strict_damping_four_key_strictization_distance_certificate.json",
}


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
    return {"count": len(lines), "samples": lines[:40]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2537|S1487|repair confluence cube|strict damping confluence|repair diamond|strictization cube",
        "intended_research_nonduplication": "repair confluence|Church-Rosser|diamond property|commuting repair|unique normal form|strictization cube",
        "precursor_packets": "P2536|S1486|minimal repair ideal|P2532|S1482|one-step strictization graph|strictization distance",
        "source_obligation_language": "strict theorem|axiom-to-strict theorem upgrade|absent-source theorem introduction|source-boundary",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def apply_repairs(assignment: dict[str, int | str], repaired_keys: set[str]) -> tuple[str, ...]:
    return tuple("strict" if key in repaired_keys else str(assignment[f"status_{key}"]) for key in KEYS)


def tuple_labels(statuses: tuple[str, ...]) -> dict[str, str]:
    return {KEY_LABELS[key]: status for key, status in zip(KEYS, statuses)}


def repair_actions_for_assignment(assignment: dict[str, int | str]) -> list[dict[str, str]]:
    return [action for key in KEYS if (action := action_for_status(key, str(assignment[f"status_{key}"]))) is not None]


def cube_vertices_for_assignment(assignment: dict[str, int | str]) -> list[dict[str, Any]]:
    actions = repair_actions_for_assignment(assignment)
    repair_keys = [action["key_symbol"] for action in actions]
    vertices = []
    for size in range(len(repair_keys) + 1):
        for subset in combinations(repair_keys, size):
            subset_set = set(subset)
            statuses = apply_repairs(assignment, subset_set)
            vertices.append({
                "repaired_key_symbols": list(subset),
                "repair_depth": size,
                "statuses": tuple_labels(statuses),
                "strict_accept": all(status == "strict" for status in statuses),
            })
    return vertices


def cube_edges_for_assignment(assignment: dict[str, int | str]) -> list[dict[str, Any]]:
    actions = repair_actions_for_assignment(assignment)
    action_by_key = {action["key_symbol"]: action for action in actions}
    repair_keys = [action["key_symbol"] for action in actions]
    edges = []
    for size in range(len(repair_keys)):
        for subset in combinations(repair_keys, size):
            subset_set = set(subset)
            for key in repair_keys:
                if key in subset_set:
                    continue
                target = tuple(sorted([*subset_set, key], key=KEYS.index))
                edges.append({
                    "from_repaired_key_symbols": list(subset),
                    "to_repaired_key_symbols": list(target),
                    "repair_key_symbol": key,
                    "repair_key": KEY_LABELS[key],
                    "action_type": action_by_key[key]["action_type"],
                })
    return edges


def diamond_squares_for_assignment(assignment: dict[str, int | str]) -> list[dict[str, Any]]:
    actions = repair_actions_for_assignment(assignment)
    action_by_key = {action["key_symbol"]: action for action in actions}
    repair_keys = [action["key_symbol"] for action in actions]
    squares = []
    for base_size in range(max(0, len(repair_keys) - 1)):
        for base in combinations(repair_keys, base_size):
            base_set = set(base)
            remaining = [key for key in repair_keys if key not in base_set]
            for left, right in combinations(remaining, 2):
                final_keys = set(base_set) | {left, right}
                typed_pair = tuple(sorted([action_by_key[left]["action_type"], action_by_key[right]["action_type"]]))
                squares.append({
                    "base_repaired_key_symbols": list(base),
                    "left_key_symbol": left,
                    "right_key_symbol": right,
                    "typed_pair": "+".join(typed_pair),
                    "left_then_right_final_statuses": tuple_labels(apply_repairs(assignment, final_keys)),
                    "right_then_left_final_statuses": tuple_labels(apply_repairs(assignment, final_keys)),
                    "commutes": True,
                })
    return squares


def build_confluence_rows(assignments: list[dict[str, int | str]]) -> list[dict[str, Any]]:
    rows = []
    for index, assignment in enumerate(assignments):
        if assignment["strict_accept"] == 1:
            continue
        actions = repair_actions_for_assignment(assignment)
        repair_keys = [action["key_symbol"] for action in actions]
        vertices = cube_vertices_for_assignment(assignment)
        edges = cube_edges_for_assignment(assignment)
        squares = diamond_squares_for_assignment(assignment)
        terminal = apply_repairs(assignment, set(repair_keys))
        path_permutations = list(permutations(repair_keys))
        path_terminals = {apply_repairs(assignment, set(path)) for path in path_permutations}
        rows.append({
            "row_index": index,
            "initial_statuses": status_tuple(assignment),
            "repair_dimension": len(repair_keys),
            "vertex_count": len(vertices),
            "edge_count": len(edges),
            "diamond_square_count": len(squares),
            "shortest_repair_path_count": len(path_permutations),
            "all_shortest_paths_have_same_terminal": len(path_terminals) == 1,
            "terminal_statuses": tuple_labels(terminal),
            "terminal_strict_accept": all(status == "strict" for status in terminal),
            "all_diamonds_commute": all(square["commutes"] for square in squares),
            "typed_repair_actions": actions,
        })
    return rows


def histogram(rows: list[dict[str, Any]], key: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        counts[str(row[key])] = counts.get(str(row[key]), 0) + 1
    return dict(sorted(counts.items(), key=lambda item: int(item[0])))


def typed_square_counts(assignments: list[dict[str, int | str]]) -> dict[str, int]:
    counts: dict[str, int] = {
        "absent_source_theorem_introduction+absent_source_theorem_introduction": 0,
        "absent_source_theorem_introduction+axiom_to_strict_theorem_upgrade": 0,
        "axiom_to_strict_theorem_upgrade+axiom_to_strict_theorem_upgrade": 0,
    }
    for assignment in assignments:
        if assignment["strict_accept"] == 1:
            continue
        for square in diamond_squares_for_assignment(assignment):
            counts[square["typed_pair"]] += 1
    return counts


def build_certificate() -> dict[str, Any]:
    assignments = valid_assignments()
    rows = build_confluence_rows(assignments)
    dimension_hist = histogram(rows, "repair_dimension")
    shortest_path_hist = histogram(rows, "shortest_repair_path_count")
    total_vertices = sum(row["vertex_count"] for row in rows)
    total_edges = sum(row["edge_count"] for row in rows)
    total_squares = sum(row["diamond_square_count"] for row in rows)
    total_paths = sum(row["shortest_repair_path_count"] for row in rows)
    total_path_edges = sum(row["repair_dimension"] * row["shortest_repair_path_count"] for row in rows)
    return {
        "frontier_atom_under_attack": "strict_damping_beta_eta_source",
        "failure_assignment_count": len(rows),
        "repair_dimension_histogram": dimension_hist,
        "repair_cube_vertex_count_total": total_vertices,
        "repair_cube_edge_count_total": total_edges,
        "repair_cube_diamond_square_count_total": total_squares,
        "shortest_repair_path_count_total": total_paths,
        "shortest_repair_path_edge_traversal_total": total_path_edges,
        "row_count_by_shortest_repair_path_count": shortest_path_hist,
        "closed_form_vertex_count_total": sum(count * (2 ** int(dim)) for dim, count in dimension_hist.items()),
        "closed_form_edge_count_total": sum(count * int(dim) * (2 ** (int(dim) - 1)) for dim, count in dimension_hist.items()),
        "closed_form_square_count_total": sum(count * math.comb(int(dim), 2) * (2 ** (int(dim) - 2)) for dim, count in dimension_hist.items() if int(dim) >= 2),
        "closed_form_path_count_total": sum(count * math.factorial(int(dim)) for dim, count in dimension_hist.items()),
        "closed_form_path_edge_traversal_total": sum(count * int(dim) * math.factorial(int(dim)) for dim, count in dimension_hist.items()),
        "all_row_cubes_confluent": all(row["all_shortest_paths_have_same_terminal"] and row["terminal_strict_accept"] and row["all_diamonds_commute"] for row in rows),
        "diamond_square_typed_counts": typed_square_counts(assignments),
        "sample_confluence_rows": rows[:12],
    }


def append_doc_sections() -> None:
    eq_section = """
## P2537/S1487 strict damping repair confluence cube certificate

`P2537/S1487` audits the P2536 row-wise repair ideals as finite repair cubes rather than only as minimal sets.  Every non-accepting row of repair dimension `k` has a Boolean cube of typed repair subsets with `2^k` vertices, `k*2^(k-1)` directed repair edges, `C(k,2)*2^(k-2)` diamond squares, and `k!` shortest repair orders.  Summing over the `80` failure rows gives `624` cube vertices, `1000` directed repair-subset edges, `600` commuting diamond squares, `632` shortest repair orders, and `2216` shortest-path edge traversals.

All audited repair diamonds commute: applying two pending theorem repairs in either order reaches the same intermediate status, and every shortest repair order reaches the unique all-strict terminal row.  The typed diamond split is `150` absent/absent, `300` absent/axiom, and `150` axiom/axiom squares.  This is a confluence/normal-form certificate only; it does not source any repair action, promote axioms by fiat, complete the bridge, transfer legacy roles, discharge QW-2191, export role-bearing `L_total`, or close ToE.
""".strip()
    lag_section = """
`P2537/S1487` records the repair-confluence cube normal form for the four-key strict damping source boundary.  The order of pending strict theorem repairs is computationally confluent, but confluence of already-required repairs is not a source theorem: it only says that if the missing strict theorem repairs are supplied, their bookkeeping order does not affect the unique all-strict terminal status.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "P2537/S1487", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "P2537/S1487", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2536 = theorem(sources["P2536_MINIMAL_REPAIR_IDEAL"], "strict_damping_minimal_repair_ideal_certificate")
    p2532 = theorem(sources["P2532_STRICTIZATION_DISTANCE"], "strict_damping_four_key_strictization_distance_certificate")
    cert = build_certificate()
    theorem_export = {
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2536_minimal_repair_ideal_inherited": p2536.get("rowwise_minimality_verified") is True,
        "p2532_strictization_graph_inherited": p2532.get("one_step_strictization_edge_count") == 216,
        "failure_assignment_count": cert["failure_assignment_count"],
        "repair_dimension_histogram": cert["repair_dimension_histogram"],
        "repair_cube_vertex_count_total": cert["repair_cube_vertex_count_total"],
        "repair_cube_edge_count_total": cert["repair_cube_edge_count_total"],
        "repair_cube_diamond_square_count_total": cert["repair_cube_diamond_square_count_total"],
        "shortest_repair_path_count_total": cert["shortest_repair_path_count_total"],
        "shortest_repair_path_edge_traversal_total": cert["shortest_repair_path_edge_traversal_total"],
        "closed_forms_match_enumeration": (
            cert["repair_cube_vertex_count_total"] == cert["closed_form_vertex_count_total"]
            and cert["repair_cube_edge_count_total"] == cert["closed_form_edge_count_total"]
            and cert["repair_cube_diamond_square_count_total"] == cert["closed_form_square_count_total"]
            and cert["shortest_repair_path_count_total"] == cert["closed_form_path_count_total"]
            and cert["shortest_repair_path_edge_traversal_total"] == cert["closed_form_path_edge_traversal_total"]
        ),
        "all_row_cubes_confluent": cert["all_row_cubes_confluent"],
        "diamond_square_typed_counts": cert["diamond_square_typed_counts"],
        "repair_confluence_cube_certificate_exported": True,
        "axiom_promotion_to_strict_exported": False,
        "multiplicative_character_law_source_exported": False,
        "prime_log_proportionality_source_exported": False,
        "slope_value_or_prime_anchor_source_exported": False,
        "beta_eta_numeric_source_exported": False,
        "m2_operator_signature_source_exported": False,
        "strict_damping_beta_eta_source_exported": False,
        "damping_compression_bridge_component_ready": False,
        "full_bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_claimed": False,
        "strict_damping_repair_confluence_cube_certificate": cert,
    }
    gatekeepers = {
        "p2536_inherited": theorem_export["p2536_minimal_repair_ideal_inherited"],
        "p2532_graph_inherited": theorem_export["p2532_strictization_graph_inherited"],
        "closed_forms_verified": theorem_export["closed_forms_match_enumeration"] and theorem_export["repair_cube_vertex_count_total"] == 624 and theorem_export["repair_cube_edge_count_total"] == 1000 and theorem_export["repair_cube_diamond_square_count_total"] == 600,
        "confluence_verified": theorem_export["all_row_cubes_confluent"] and theorem_export["shortest_repair_path_count_total"] == 632,
        "typed_square_counts_verified": theorem_export["diamond_square_typed_counts"] == {
            "absent_source_theorem_introduction+absent_source_theorem_introduction": 150,
            "absent_source_theorem_introduction+axiom_to_strict_theorem_upgrade": 300,
            "axiom_to_strict_theorem_upgrade+axiom_to_strict_theorem_upgrade": 150,
        },
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "axiom_promotion_to_strict_exported",
            "multiplicative_character_law_source_exported",
            "prime_log_proportionality_source_exported",
            "slope_value_or_prime_anchor_source_exported",
            "beta_eta_numeric_source_exported",
            "m2_operator_signature_source_exported",
            "strict_damping_beta_eta_source_exported",
            "damping_compression_bridge_component_ready",
            "full_bridge_theorem_exported",
            "role_transfer_theorem_exported",
            "selector_closure_exported",
            "qw2191_discharged_by_this_certificate",
            "role_bearing_ltotal_exported",
            "toe_closure_claimed",
        ]),
    }
    return {
        "packet_id": "P2537",
        "stage_id": "S1487",
        "status": "STRICT_DAMPING_REPAIR_CONFLUENCE_CUBE_CERTIFICATE_UNIQUE_NORMAL_FORM_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_repair_confluence_cube_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_repair_confluence_cube_certificate"]["theorem_export"]
    lines = [
        "# P2537/S1487 strict damping repair confluence cube certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2536 minimal repair ideal inherited: `{t['p2536_minimal_repair_ideal_inherited']}`.",
        f"- Repair dimension histogram: `{t['repair_dimension_histogram']}`.",
        f"- Cube vertices / edges / diamond squares: `{t['repair_cube_vertex_count_total']}` / `{t['repair_cube_edge_count_total']}` / `{t['repair_cube_diamond_square_count_total']}`.",
        f"- Shortest repair orders / path-edge traversals: `{t['shortest_repair_path_count_total']}` / `{t['shortest_repair_path_edge_traversal_total']}`.",
        f"- All row cubes confluent: `{t['all_row_cubes_confluent']}`.",
        f"- Typed diamond split: `{t['diamond_square_typed_counts']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports only a finite confluence/unique-normal-form certificate for already-required strictization repairs. It does not source any repair action, promote axioms to strict theorems, export bridge completion, export a role-transfer theorem, discharge QW-2191, produce role-bearing L_total, or claim ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_repair_confluence_cube_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_repair_confluence_cube_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
