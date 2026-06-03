#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from collections import Counter
from math import comb, factorial
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2532_s1482_strict_damping_four_key_strictization_distance_certificate.json"
MD = GEN / "p2532_s1482_strict_damping_four_key_strictization_distance_certificate.md"

SOURCE_FILES = {
    "P2531_FOUR_KEY_AXIOM_BOUNDARY": GEN / "p2531_s1481_strict_damping_four_key_axiom_boundary_certificate.json",
}

SOURCE_KEYS = [
    "M_multiplicative_character_law_source",
    "P_prime_log_proportionality_source",
    "A_slope_value_or_prime_anchor_source",
    "O_m2_operator_signature_source",
]
STATUSES = ["absent", "axiom", "strict"]


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
        "new_packet": "P2532|S1482|four-key strictization distance|strictization-distance|axiom upgrade distance|theorem deficit stratification|one-step strictization graph",
        "intended_research_nonduplication": "strictization distance|axiom upgrade|theorem deficit|upgrade distance|strict theorem deficit|deficit stratification|one-step strictization|strictization graph",
        "precursor_packets": "P2531|S1481|four-key axiom boundary|P2530|four-key irredundancy|P2529|numeric subkey rank lattice",
        "source_obligation_language": "axiom-augmented|non-strict|strict theorem|absent key|strict_damping_beta_eta_source",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def row_counts(row: dict[str, Any]) -> dict[str, int]:
    statuses = row["status_by_key"].values()
    counts = Counter(statuses)
    return {status: counts.get(status, 0) for status in STATUSES}


def strictization_profile(row: dict[str, Any]) -> dict[str, Any]:
    counts = row_counts(row)
    absent_count = counts["absent"]
    axiom_count = counts["axiom"]
    strict_count = counts["strict"]
    theorem_deficit = absent_count + axiom_count
    if row["strict_damping_source_strict_accepts"]:
        class_name = "strict_accept"
        strict_upgrade_distance_if_present = 0
        absent_source_gap = 0
    elif row["strict_damping_source_axiom_augmented_only"]:
        class_name = "present_non_strict_axiom_augmented"
        strict_upgrade_distance_if_present = axiom_count
        absent_source_gap = 0
    else:
        class_name = "blocked_missing_key"
        strict_upgrade_distance_if_present = None
        absent_source_gap = absent_count
    return {
        "status_by_key": row["status_by_key"],
        "strict_count": strict_count,
        "axiom_count": axiom_count,
        "absent_count": absent_count,
        "theorem_deficit_to_all_strict": theorem_deficit,
        "absent_source_gap": absent_source_gap,
        "strict_upgrade_distance_if_all_keys_present": strict_upgrade_distance_if_present,
        "strict_damping_source_verdict": row["strict_damping_source_verdict"],
        "strictization_class": class_name,
    }


def stratify(lattice: list[dict[str, Any]]) -> dict[str, Any]:
    profiles = [strictization_profile(row) for row in lattice]
    by_class = Counter(profile["strictization_class"] for profile in profiles)
    by_theorem_deficit = Counter(profile["theorem_deficit_to_all_strict"] for profile in profiles)
    present_axiom_by_upgrade_distance = Counter(
        profile["strict_upgrade_distance_if_all_keys_present"]
        for profile in profiles
        if profile["strictization_class"] == "present_non_strict_axiom_augmented"
    )
    blocked_by_absent_gap = Counter(
        profile["absent_source_gap"]
        for profile in profiles
        if profile["strictization_class"] == "blocked_missing_key"
    )
    expected_deficit_counts = {}
    for k in range(len(SOURCE_KEYS) + 1):
        # Pick k non-strict coordinates and choose absent/axiom for each.
        expected_deficit_counts[str(k)] = comb(len(SOURCE_KEYS), k) * (2 ** k)
    minimal_present_axiom_upgrade_rows = [
        profile for profile in profiles
        if profile["strictization_class"] == "present_non_strict_axiom_augmented"
        and profile["strict_upgrade_distance_if_all_keys_present"] == 1
    ]
    maximal_absent_blocker_rows = [
        profile for profile in profiles
        if profile["strictization_class"] == "blocked_missing_key"
        and profile["absent_count"] == len(SOURCE_KEYS)
    ]
    return {
        "strictization_profiles": profiles,
        "class_counts": dict(sorted(by_class.items())),
        "theorem_deficit_counts": {str(k): by_theorem_deficit[k] for k in range(len(SOURCE_KEYS) + 1)},
        "expected_theorem_deficit_counts": expected_deficit_counts,
        "present_axiom_upgrade_distance_counts": {str(k): present_axiom_by_upgrade_distance[k] for k in range(1, len(SOURCE_KEYS) + 1)},
        "blocked_absent_source_gap_counts": {str(k): blocked_by_absent_gap[k] for k in range(1, len(SOURCE_KEYS) + 1)},
        "minimal_present_axiom_upgrade_rows": minimal_present_axiom_upgrade_rows,
        "maximal_absent_blocker_rows": maximal_absent_blocker_rows,
    }


def strictization_edges(profile: dict[str, Any]) -> list[dict[str, Any]]:
    edges = []
    for key in SOURCE_KEYS:
        source_status = profile["status_by_key"][key]
        if source_status == "strict":
            continue
        target_status_by_key = dict(profile["status_by_key"])
        target_status_by_key[key] = "strict"
        edge_type = "absent_source_theorem_introduction" if source_status == "absent" else "axiom_to_strict_theorem_upgrade"
        edges.append({
            "upgraded_key": key,
            "source_status": source_status,
            "target_status": "strict",
            "edge_type": edge_type,
            "deficit_before": profile["theorem_deficit_to_all_strict"],
            "deficit_after": profile["theorem_deficit_to_all_strict"] - 1,
            "source_strictization_class": profile["strictization_class"],
            "target_status_by_key": target_status_by_key,
        })
    return edges


def strictization_graph_summary(profiles: list[dict[str, Any]]) -> dict[str, Any]:
    all_edges = []
    path_count_by_deficit: Counter[int] = Counter()
    for profile in profiles:
        deficit = profile["theorem_deficit_to_all_strict"]
        path_count_by_deficit[deficit] += factorial(deficit)
        all_edges.extend(strictization_edges(profile))

    edge_type_counts = Counter(edge["edge_type"] for edge in all_edges)
    edge_source_class_counts = Counter(edge["source_strictization_class"] for edge in all_edges)
    upgraded_key_counts = Counter(edge["upgraded_key"] for edge in all_edges)
    return {
        "one_step_edge_count": len(all_edges),
        "one_step_edge_type_counts": dict(sorted(edge_type_counts.items())),
        "one_step_edge_source_class_counts": dict(sorted(edge_source_class_counts.items())),
        "upgraded_key_counts": {key: upgraded_key_counts[key] for key in SOURCE_KEYS},
        "shortest_strictization_path_count_by_deficit": {str(k): path_count_by_deficit[k] for k in range(len(SOURCE_KEYS) + 1)},
        "total_shortest_strictization_path_count": sum(path_count_by_deficit.values()),
        "edge_samples": all_edges[:16],
    }


def append_doc_sections() -> None:
    eq_section = """
`P2532/S1482` refines the P2531 ternary axiom boundary by measuring strictization distance on the full `3^4=81` four-key table.  The unique all-strict row has theorem deficit `0`; every other row has positive theorem deficit.  The 15 all-present axiom-augmented rows stratify by axiom-upgrade distance `1,2,3,4` with counts `4,6,4,1`, while the 65 blocked rows stratify by absent source gap `1,2,3,4` with counts `32,24,8,1`.  A one-step strictization graph has `216` directed upgrade edges: `108` absent-source theorem introductions and `108` axiom-to-strict theorem upgrades.  Thus the nearest non-strict rows are exactly the four one-axiom rows, and even these require one strict theorem upgrade rather than an axiom promotion.
""".strip()
    lag_section = """
`P2532/S1482` adds a strictization-distance audit for the four-key strict damping source boundary.  It does not reduce the source burden: all non-accepting rows have positive theorem deficit, the one-step strictization graph only records required theorem-upgrade edges, axiom-present rows need strict theorem upgrades, and missing-key rows need absent source theorems before strict damping can even be considered.  No nonlinear compression-flow source, role-bearing `L_total`, bridge completion, role-transfer theorem, QW-2191 discharge, or ToE closure is licensed.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "P2532/S1482", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "P2532/S1482", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2531 = theorem(sources["P2531_FOUR_KEY_AXIOM_BOUNDARY"], "strict_damping_four_key_axiom_boundary_certificate")
    lattice = p2531["strict_damping_four_key_axiom_boundary_certificate"]["ternary_status_lattice"]
    strata = stratify(lattice)
    graph = strictization_graph_summary(strata["strictization_profiles"])
    theorem_export = {
        "frontier_atom_under_attack": "four_key_strictization_distance_and_theorem_deficit_stratification",
        "p2531_axiom_boundary_inherited": bool(p2531.get("four_key_axiom_boundary_exported", False)),
        "source_keys": SOURCE_KEYS,
        "status_values": STATUSES,
        "ternary_status_row_count": len(lattice),
        "strict_accept_row_count": strata["class_counts"].get("strict_accept", 0),
        "present_non_strict_axiom_augmented_row_count": strata["class_counts"].get("present_non_strict_axiom_augmented", 0),
        "blocked_missing_key_row_count": strata["class_counts"].get("blocked_missing_key", 0),
        "theorem_deficit_counts": strata["theorem_deficit_counts"],
        "theorem_deficit_counts_match_closed_form": strata["theorem_deficit_counts"] == strata["expected_theorem_deficit_counts"],
        "present_axiom_upgrade_distance_counts": strata["present_axiom_upgrade_distance_counts"],
        "present_axiom_upgrade_distance_counts_match_binomial": strata["present_axiom_upgrade_distance_counts"] == {"1": 4, "2": 6, "3": 4, "4": 1},
        "blocked_absent_source_gap_counts": strata["blocked_absent_source_gap_counts"],
        "blocked_absent_source_gap_counts_match_enumeration": strata["blocked_absent_source_gap_counts"] == {"1": 32, "2": 24, "3": 8, "4": 1},
        "minimal_non_strict_present_rows_are_one_axiom_rows": len(strata["minimal_present_axiom_upgrade_rows"]) == 4 and all(profile["axiom_count"] == 1 for profile in strata["minimal_present_axiom_upgrade_rows"]),
        "maximal_absent_blocker_row_count": len(strata["maximal_absent_blocker_rows"]),
        "one_step_strictization_edge_count": graph["one_step_edge_count"],
        "one_step_strictization_edge_count_matches_deficit_sum": graph["one_step_edge_count"] == 216,
        "one_step_edge_type_counts": graph["one_step_edge_type_counts"],
        "one_step_edge_type_counts_match_closed_form": graph["one_step_edge_type_counts"] == {"absent_source_theorem_introduction": 108, "axiom_to_strict_theorem_upgrade": 108},
        "one_step_edge_source_class_counts": graph["one_step_edge_source_class_counts"],
        "one_step_edge_source_class_counts_match_strata": graph["one_step_edge_source_class_counts"] == {"blocked_missing_key": 184, "present_non_strict_axiom_augmented": 32},
        "upgraded_key_counts": graph["upgraded_key_counts"],
        "upgraded_key_counts_uniform": all(count == 54 for count in graph["upgraded_key_counts"].values()),
        "shortest_strictization_path_count_by_deficit": graph["shortest_strictization_path_count_by_deficit"],
        "shortest_strictization_path_counts_match_factorial_closed_form": graph["shortest_strictization_path_count_by_deficit"] == {"0": 1, "1": 8, "2": 48, "3": 192, "4": 384},
        "all_non_accepting_rows_have_positive_theorem_deficit": all(profile["theorem_deficit_to_all_strict"] > 0 for profile in strata["strictization_profiles"] if profile["strictization_class"] != "strict_accept"),
        "strictization_distance_certificate_exported": True,
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
        "strict_damping_four_key_strictization_distance_certificate": {
            "class_counts": strata["class_counts"],
            "theorem_deficit_counts": strata["theorem_deficit_counts"],
            "present_axiom_upgrade_distance_counts": strata["present_axiom_upgrade_distance_counts"],
            "blocked_absent_source_gap_counts": strata["blocked_absent_source_gap_counts"],
            "minimal_present_axiom_upgrade_rows": strata["minimal_present_axiom_upgrade_rows"],
            "maximal_absent_blocker_rows": strata["maximal_absent_blocker_rows"],
            "one_step_strictization_graph_summary": graph,
            "source_obligation_normal_form": "strict damping source acceptance is at strictization distance zero only; axiom rows require strict theorem upgrades and absent rows require missing source theorems",
        },
    }
    gatekeepers = {
        "p2531_inherited": theorem_export["p2531_axiom_boundary_inherited"],
        "strictization_counts_verified": theorem_export["theorem_deficit_counts_match_closed_form"] and theorem_export["present_axiom_upgrade_distance_counts_match_binomial"] and theorem_export["blocked_absent_source_gap_counts_match_enumeration"],
        "one_step_graph_verified": theorem_export["one_step_strictization_edge_count_matches_deficit_sum"] and theorem_export["one_step_edge_type_counts_match_closed_form"] and theorem_export["one_step_edge_source_class_counts_match_strata"] and theorem_export["upgraded_key_counts_uniform"] and theorem_export["shortest_strictization_path_counts_match_factorial_closed_form"],
        "positive_deficit_for_all_non_accepting_rows": theorem_export["all_non_accepting_rows_have_positive_theorem_deficit"],
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
        "packet_id": "P2532",
        "stage_id": "S1482",
        "status": "STRICT_DAMPING_FOUR_KEY_STRICTIZATION_DISTANCE_CERTIFICATE_DEFICIT_STRATIFICATION_AND_GRAPH_ONLY_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_four_key_strictization_distance_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_four_key_strictization_distance_certificate"]["theorem_export"]
    lines = [
        "# P2532/S1482 strict damping four-key strictization distance certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2531 axiom boundary inherited: `{t['p2531_axiom_boundary_inherited']}`.",
        f"- Ternary status rows: `{t['ternary_status_row_count']}`.",
        f"- Strict accept rows: `{t['strict_accept_row_count']}`.",
        f"- Present non-strict axiom-augmented rows: `{t['present_non_strict_axiom_augmented_row_count']}`.",
        f"- Blocked missing-key rows: `{t['blocked_missing_key_row_count']}`.",
        f"- Theorem-deficit counts: `{t['theorem_deficit_counts']}`.",
        f"- Present axiom upgrade-distance counts: `{t['present_axiom_upgrade_distance_counts']}`.",
        f"- Blocked absent source-gap counts: `{t['blocked_absent_source_gap_counts']}`.",
        f"- All non-accepting rows have positive theorem deficit: `{t['all_non_accepting_rows_have_positive_theorem_deficit']}`.",
        f"- One-step strictization edges: `{t['one_step_strictization_edge_count']}`.",
        f"- One-step edge type counts: `{t['one_step_edge_type_counts']}`.",
        f"- Shortest strictization path counts by deficit: `{t['shortest_strictization_path_count_by_deficit']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports only a strictization-distance stratification of the conditional four-key source table. It does not promote axioms to strict theorems, source any key, export bridge completion, export a role-transfer theorem, discharge QW-2191, produce role-bearing L_total, or claim ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_four_key_strictization_distance_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_four_key_strictization_distance_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
