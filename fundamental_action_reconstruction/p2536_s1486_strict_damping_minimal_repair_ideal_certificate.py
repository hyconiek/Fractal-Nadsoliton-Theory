#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from itertools import combinations
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

GEN = ROOT / "generated"
OUT = GEN / "p2536_s1486_strict_damping_minimal_repair_ideal_certificate.json"
MD = GEN / "p2536_s1486_strict_damping_minimal_repair_ideal_certificate.md"

SOURCE_FILES = {
    "P2535_DUAL_FAILURE_COVER": GEN / "p2535_s1485_strict_damping_dual_failure_cover_certificate.json",
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
        "new_packet": "P2536|S1486|minimal repair ideal|strict damping repair ideal|minimal strictization repair|repair-set ideal",
        "intended_research_nonduplication": "minimal repair ideal|strictization repair|repair action|proper repair subset|minimal repair set|repair subset lattice",
        "precursor_packets": "P2535|S1485|dual failure cover|P2532|S1482|strictization distance|P2534|Boolean prime implicant",
        "source_obligation_language": "strict theorem|axiom-to-strict theorem upgrade|absent-source theorem introduction|source-boundary",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def action_for_status(key: str, status: str) -> dict[str, str] | None:
    if status == "strict":
        return None
    if status == "axiom":
        action_type = "axiom_to_strict_theorem_upgrade"
    elif status == "absent":
        action_type = "absent_source_theorem_introduction"
    else:
        raise ValueError(f"unknown status {status!r}")
    return {
        "key": KEY_LABELS[key],
        "key_symbol": key,
        "from_status": status,
        "to_status": "strict",
        "action_type": action_type,
    }


def repaired_statuses(assignment: dict[str, int | str], repaired_keys: set[str]) -> dict[str, str]:
    statuses: dict[str, str] = {}
    for key in KEYS:
        status = str(assignment[f"status_{key}"])
        statuses[key] = "strict" if key in repaired_keys else status
    return statuses


def is_strict_accept_statuses(statuses: dict[str, str]) -> bool:
    return all(statuses[key] == "strict" for key in KEYS)


def powerset(items: list[str]) -> list[tuple[str, ...]]:
    return [subset for size in range(len(items) + 1) for subset in combinations(items, size)]


def multinomial4(absent: int, axiom: int, strict: int) -> int:
    if absent + axiom + strict != 4:
        raise ValueError("expected four keys")
    return math.factorial(4) // (math.factorial(absent) * math.factorial(axiom) * math.factorial(strict))


def build_repair_rows(assignments: list[dict[str, int | str]]) -> list[dict[str, Any]]:
    rows = []
    for index, assignment in enumerate(assignments):
        if assignment["strict_accept"] == 1:
            continue
        actions = [action for key in KEYS if (action := action_for_status(key, str(assignment[f"status_{key}"]))) is not None]
        repair_keys = [action["key_symbol"] for action in actions]
        proper_subsets = [set(subset) for subset in powerset(repair_keys) if set(subset) != set(repair_keys)]
        accepting_proper_subsets = [sorted(subset) for subset in proper_subsets if is_strict_accept_statuses(repaired_statuses(assignment, subset))]
        full_repair_accepts = is_strict_accept_statuses(repaired_statuses(assignment, set(repair_keys)))
        rows.append({
            "row_index": index,
            "initial_statuses": status_tuple(assignment),
            "absent_count": sum(assignment[f"status_{key}"] == "absent" for key in KEYS),
            "axiom_count": sum(assignment[f"status_{key}"] == "axiom" for key in KEYS),
            "strict_count": sum(assignment[f"status_{key}"] == "strict" for key in KEYS),
            "theorem_deficit": len(actions),
            "minimal_repair_actions": actions,
            "minimal_repair_action_types": [action["action_type"] for action in actions],
            "proper_repair_subset_count_including_empty": len(proper_subsets),
            "accepting_proper_repair_subset_count": len(accepting_proper_subsets),
            "full_repair_accepts": full_repair_accepts,
        })
    return rows


def histogram(rows: list[dict[str, Any]], key: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        counts[str(row[key])] = counts.get(str(row[key]), 0) + 1
    return dict(sorted(counts.items(), key=lambda item: int(item[0])))


def repair_bigrade_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result = []
    for absent in range(5):
        for axiom in range(5 - absent):
            strict = 4 - absent - axiom
            if absent == 0 and axiom == 0:
                continue
            matching = [row for row in rows if row["absent_count"] == absent and row["axiom_count"] == axiom and row["strict_count"] == strict]
            result.append({
                "absent_source_theorem_introductions": absent,
                "axiom_to_strict_theorem_upgrades": axiom,
                "already_strict_keys": strict,
                "minimal_repair_set_size": absent + axiom,
                "row_count": len(matching),
                "closed_form_multinomial_count": multinomial4(absent, axiom, strict),
                "total_absent_introduction_actions": absent * len(matching),
                "total_axiom_upgrade_actions": axiom * len(matching),
            })
    return result


def key_action_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result = []
    for key in KEYS:
        intro_rows = [row for row in rows if any(action["key_symbol"] == key and action["action_type"] == "absent_source_theorem_introduction" for action in row["minimal_repair_actions"])]
        upgrade_rows = [row for row in rows if any(action["key_symbol"] == key and action["action_type"] == "axiom_to_strict_theorem_upgrade" for action in row["minimal_repair_actions"])]
        result.append({
            "key": KEY_LABELS[key],
            "key_symbol": key,
            "absent_source_theorem_introduction_rows": len(intro_rows),
            "axiom_to_strict_theorem_upgrade_rows": len(upgrade_rows),
            "total_minimal_repair_action_incidence": len(intro_rows) + len(upgrade_rows),
        })
    return result


def build_certificate() -> dict[str, Any]:
    assignments = valid_assignments()
    repair_rows = build_repair_rows(assignments)
    bigrades = repair_bigrade_rows(repair_rows)
    key_rows = key_action_rows(repair_rows)
    minimal_size_hist = histogram(repair_rows, "theorem_deficit")
    proper_subset_total = sum(row["proper_repair_subset_count_including_empty"] for row in repair_rows)
    return {
        "frontier_atom_under_attack": "strict_damping_beta_eta_source",
        "valid_ternary_assignment_count": len(assignments),
        "strict_accept_assignment_count": sum(assignment["strict_accept"] == 1 for assignment in assignments),
        "failure_assignment_count": len(repair_rows),
        "minimal_repair_set_count": len(repair_rows),
        "minimal_repair_set_size_histogram": minimal_size_hist,
        "minimal_repair_bigrade_rows": bigrades,
        "repair_bigrade_closed_form": "[(r_abs+r_ax+1)^4-1] with coefficients 4!/(r_abs! r_ax! s!) for r_abs+r_ax+s=4, excluding s^4",
        "rowwise_minimality_verified": all(row["full_repair_accepts"] and row["accepting_proper_repair_subset_count"] == 0 for row in repair_rows),
        "proper_repair_subset_count_including_empty": proper_subset_total,
        "proper_repair_subset_count_excluding_empty": sum(max(0, row["proper_repair_subset_count_including_empty"] - 1) for row in repair_rows),
        "full_minimal_repair_subset_count": len(repair_rows),
        "all_candidate_repair_subsets_count": proper_subset_total + len(repair_rows),
        "key_action_rows": key_rows,
        "total_absent_source_theorem_introduction_actions": sum(row["absent_source_theorem_introduction_rows"] for row in key_rows),
        "total_axiom_to_strict_theorem_upgrade_actions": sum(row["axiom_to_strict_theorem_upgrade_rows"] for row in key_rows),
        "total_minimal_repair_action_incidence": sum(row["total_minimal_repair_action_incidence"] for row in key_rows),
        "uniform_key_repair_action_incidence": len({row["total_minimal_repair_action_incidence"] for row in key_rows}) == 1,
        "uniform_absent_and_axiom_action_incidence_per_key": all(row["absent_source_theorem_introduction_rows"] == 27 and row["axiom_to_strict_theorem_upgrade_rows"] == 27 for row in key_rows),
        "sample_minimal_repair_rows": repair_rows[:12],
    }


def append_doc_sections() -> None:
    eq_section = """
## P2536/S1486 strict damping minimal repair ideal certificate

`P2536/S1486` turns the P2535 failure cover into a row-wise minimal strictization-repair theorem.  For every one of the `80` non-accepting ternary rows, the unique minimal repair set is exactly the set of non-strict keys in that row: absent keys require `absent_source_theorem_introduction`, and axiom keys require `axiom_to_strict_theorem_upgrade`.  Exhaustive subset enumeration verifies that all `544` proper repair subsets, including empty subsets, still fail, while the `80` full row-wise repair sets accept.

The repair bigrade has closed form `[(r_abs+r_ax+1)^4-1]`: its theorem-deficit histogram is `1:8, 2:24, 3:32, 4:16`, and its action incidence totals are `108` absent-source theorem introductions plus `108` axiom-to-strict theorem upgrades, uniformly `54` repair incidences per source key.  This is a repair-obligation certificate only; it does not source any key, promote axioms by fiat, complete the bridge, transfer legacy roles, discharge QW-2191, export role-bearing `L_total`, or close ToE.
""".strip()
    lag_section = """
`P2536/S1486` records the minimal repair ideal for the four-key strict damping source boundary.  The normal form is row-wise exact: a row reaches the strict source only after every non-strict key is converted to a strict theorem by its typed repair action.  Proper repair subsets remain outside the accepted `s^4` atom, so partial source repairs cannot be read as a nonlinear compression-flow source or a role-bearing Lagrangian term.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "P2536/S1486", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "P2536/S1486", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2535 = theorem(sources["P2535_DUAL_FAILURE_COVER"], "strict_damping_dual_failure_cover_certificate")
    p2532 = theorem(sources["P2532_STRICTIZATION_DISTANCE"], "strict_damping_four_key_strictization_distance_certificate")
    cert = build_certificate()
    theorem_export = {
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2535_dual_failure_cover_inherited": p2535.get("strict_theorem_failure_cover_exhausts_all_failures") is True,
        "p2532_strictization_graph_inherited": p2532.get("one_step_strictization_edge_count") == 216,
        "valid_ternary_assignment_count": cert["valid_ternary_assignment_count"],
        "strict_accept_assignment_count": cert["strict_accept_assignment_count"],
        "failure_assignment_count": cert["failure_assignment_count"],
        "minimal_repair_set_count": cert["minimal_repair_set_count"],
        "minimal_repair_set_size_histogram": cert["minimal_repair_set_size_histogram"],
        "repair_bigrade_closed_form": cert["repair_bigrade_closed_form"],
        "repair_bigrade_matches_closed_form": all(row["row_count"] == row["closed_form_multinomial_count"] for row in cert["minimal_repair_bigrade_rows"]),
        "rowwise_minimality_verified": cert["rowwise_minimality_verified"],
        "proper_repair_subset_count_including_empty": cert["proper_repair_subset_count_including_empty"],
        "proper_repair_subset_count_excluding_empty": cert["proper_repair_subset_count_excluding_empty"],
        "full_minimal_repair_subset_count": cert["full_minimal_repair_subset_count"],
        "all_candidate_repair_subsets_count": cert["all_candidate_repair_subsets_count"],
        "total_absent_source_theorem_introduction_actions": cert["total_absent_source_theorem_introduction_actions"],
        "total_axiom_to_strict_theorem_upgrade_actions": cert["total_axiom_to_strict_theorem_upgrade_actions"],
        "total_minimal_repair_action_incidence": cert["total_minimal_repair_action_incidence"],
        "uniform_key_repair_action_incidence": cert["uniform_key_repair_action_incidence"],
        "uniform_absent_and_axiom_action_incidence_per_key": cert["uniform_absent_and_axiom_action_incidence_per_key"],
        "minimal_repair_ideal_certificate_exported": True,
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
        "strict_damping_minimal_repair_ideal_certificate": cert,
    }
    gatekeepers = {
        "p2535_inherited": theorem_export["p2535_dual_failure_cover_inherited"],
        "p2532_graph_inherited": theorem_export["p2532_strictization_graph_inherited"],
        "rowwise_minimality_verified": theorem_export["rowwise_minimality_verified"] and theorem_export["proper_repair_subset_count_including_empty"] == 544,
        "repair_bigrade_verified": theorem_export["repair_bigrade_matches_closed_form"] and theorem_export["minimal_repair_set_size_histogram"] == {"1": 8, "2": 24, "3": 32, "4": 16},
        "action_incidence_matches_p2532": theorem_export["total_absent_source_theorem_introduction_actions"] == 108 and theorem_export["total_axiom_to_strict_theorem_upgrade_actions"] == 108 and theorem_export["total_minimal_repair_action_incidence"] == 216,
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
        "packet_id": "P2536",
        "stage_id": "S1486",
        "status": "STRICT_DAMPING_MINIMAL_REPAIR_IDEAL_CERTIFICATE_ROWWISE_EXACT_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_minimal_repair_ideal_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_minimal_repair_ideal_certificate"]["theorem_export"]
    lines = [
        "# P2536/S1486 strict damping minimal repair ideal certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2535 dual failure cover inherited: `{t['p2535_dual_failure_cover_inherited']}`.",
        f"- P2532 strictization graph inherited: `{t['p2532_strictization_graph_inherited']}`.",
        f"- Failure rows / minimal repair sets: `{t['failure_assignment_count']}` / `{t['minimal_repair_set_count']}`.",
        f"- Minimal repair size histogram: `{t['minimal_repair_set_size_histogram']}`.",
        f"- Proper repair subsets still failing, including empty subsets: `{t['proper_repair_subset_count_including_empty']}`.",
        f"- Full row-wise repair subsets accepting: `{t['full_minimal_repair_subset_count']}`.",
        f"- Typed repair action incidence: `{t['total_absent_source_theorem_introduction_actions']}` absent-source introductions + `{t['total_axiom_to_strict_theorem_upgrade_actions']}` axiom-to-strict upgrades = `{t['total_minimal_repair_action_incidence']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports only a row-wise minimal strictization-repair ideal for the conditional four-key source boundary. It does not source any key, promote axioms to strict theorems, export bridge completion, export a role-transfer theorem, discharge QW-2191, produce role-bearing L_total, or claim ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_minimal_repair_ideal_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_minimal_repair_ideal_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
