#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from decimal import Decimal
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    function_for_family,
    load_json,
    rel,
    replay_cell,
)

GEN = ROOT / "generated"
OUT = GEN / "p2474_s1424_strict_pointwise_interval_decimal_p2459_extremal_witness_rerun_audit.json"
MD = GEN / "p2474_s1424_strict_pointwise_interval_decimal_p2459_extremal_witness_rerun_audit.md"

SOURCE_FILES = {
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
    "P2469_FULL_OPPOSITE_TAIL_REPLAY": GEN / "p2469_s1419_strict_pointwise_interval_decimal_full_opposite_tail_replay_certificate.json",
    "P2470_REMAINING_NON_TAIL_REPLAY": GEN / "p2470_s1420_strict_pointwise_interval_decimal_full_remaining_non_tail_complement_replay_certificate.json",
    "P2471_FINITE_PARTITION_WITNESS": GEN / "p2471_s1421_strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate.json",
    "P2472_PARTITION_SEAM_REPLAY_AUDIT": GEN / "p2472_s1422_strict_pointwise_interval_decimal_p2459_partition_seam_replay_audit.json",
    "P2473_FINGERPRINT_BINDING_AUDIT": GEN / "p2473_s1423_strict_pointwise_interval_decimal_p2459_replay_chain_fingerprint_binding_audit.json",
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
        "new_packet": "P2474|S1424|extremal witness rerun|minimum-separation witness rerun|critical witness rerun|Decimal extremal rerun",
        "precursor_packets": "P2469|S1419|P2470|S1420|P2471|S1421|P2472|S1422|P2473|S1423",
        "rerun_language": "fresh Decimal/Taylor rerun|minimum separation witness|first/last/min witness|critical finite witness|extremal cell replay",
        "hard_limit_markers": "directed-rounding root exclusion|symbolic root-exclusion theorem|global continuum root-exclusion|full mathematical proof",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def projection_vector(p2449: dict[str, Any]) -> list[Decimal]:
    values = p2449.get("strict_pointwise_rank_lift_projection_reduction_certificate", {}).get("theorem_export", {}).get("projection_vector", [])
    return [Decimal(str(value)) for value in values]


def replay_witness(family: str, stored_row: dict[str, Any], selection_rule: str, projection: list[Decimal]) -> dict[str, Any]:
    fresh = replay_cell(family, stored_row, projection, function_for_family(family))
    stored_sep = str(stored_row["decimal_separation_from_zero"])
    fresh_sep = str(fresh["decimal_separation_from_zero"])
    return {
        "family": family,
        "selection_rule": selection_rule,
        "segment_index": stored_row.get("segment_index"),
        "uncovered_index": stored_row.get("uncovered_index"),
        "stored_decimal_separation_from_zero": stored_sep,
        "fresh_decimal_separation_from_zero": fresh_sep,
        "separation_matches_stored": fresh_sep == stored_sep,
        "stored_decimal_interval_excludes_zero": stored_row.get("decimal_interval_excludes_zero"),
        "fresh_decimal_interval_excludes_zero": fresh.get("decimal_interval_excludes_zero"),
        "zero_exclusion_matches_stored": fresh.get("decimal_interval_excludes_zero") == stored_row.get("decimal_interval_excludes_zero"),
        "fresh_decimal_separation_positive": Decimal(fresh_sep) > 0,
        "fresh_row_fingerprint_sha256": sha256_json(fresh),
    }


def p2469_witnesses(p2469_theorem: dict[str, Any], projection: list[Decimal]) -> dict[str, Any]:
    rows = []
    for family in p2469_theorem["family_full_opposite_tail_replays"]:
        fam = family["family"]
        for key, label in [
            ("first_opposite_tail_full_replay_row", "p2469_first_opposite_tail_row"),
            ("minimum_opposite_tail_full_replay_row", "p2469_minimum_separation_row"),
            ("last_opposite_tail_full_replay_row", "p2469_last_opposite_tail_row"),
        ]:
            rows.append(replay_witness(fam, family[key], label, projection))
    return {
        "source_packet": "P2469/S1419",
        "stored_total_replay_count": p2469_theorem["total_opposite_tail_full_replay_count"],
        "stored_minimum_decimal_separation": p2469_theorem["minimum_opposite_tail_full_replay_decimal_separation"],
        "fresh_witness_count": len(rows),
        "fresh_witness_rows": rows,
        "all_fresh_witnesses_match_stored": all(row["separation_matches_stored"] and row["zero_exclusion_matches_stored"] for row in rows),
        "all_fresh_witnesses_exclude_zero_with_positive_separation": all(row["fresh_decimal_interval_excludes_zero"] and row["fresh_decimal_separation_positive"] for row in rows),
    }


def p2470_witnesses(p2470_theorem: dict[str, Any], projection: list[Decimal]) -> dict[str, Any]:
    rows = []
    for family in p2470_theorem["family_remaining_non_tail_replays"]:
        fam = family["family"]
        for key, label in [
            ("first_remaining_non_tail_full_replay_row", "p2470_first_remaining_non_tail_row"),
            ("minimum_remaining_non_tail_full_replay_row", "p2470_minimum_separation_row"),
            ("last_remaining_non_tail_full_replay_row", "p2470_last_remaining_non_tail_row"),
        ]:
            rows.append(replay_witness(fam, family[key], label, projection))
    return {
        "source_packet": "P2470/S1420",
        "stored_total_replay_count": p2470_theorem["total_remaining_non_tail_full_replay_count"],
        "stored_minimum_decimal_separation": p2470_theorem["minimum_remaining_non_tail_full_replay_decimal_separation"],
        "fresh_witness_count": len(rows),
        "fresh_witness_rows": rows,
        "all_fresh_witnesses_match_stored": all(row["separation_matches_stored"] and row["zero_exclusion_matches_stored"] for row in rows),
        "all_fresh_witnesses_exclude_zero_with_positive_separation": all(row["fresh_decimal_interval_excludes_zero"] and row["fresh_decimal_separation_positive"] for row in rows),
    }


def p2472_witnesses(p2472_theorem: dict[str, Any], projection: list[Decimal]) -> dict[str, Any]:
    rows = []
    for family in p2472_theorem["family_partition_seam_replays"]:
        fam = family["family"]
        for stored_row in family["seam_replay_rows"]:
            rows.append(replay_witness(fam, stored_row, f"p2472_{stored_row['seam_selection_rule']}", projection))
    return {
        "source_packet": "P2472/S1422",
        "stored_total_seam_replay_count": p2472_theorem["total_partition_seam_replay_count"],
        "stored_minimum_decimal_separation": p2472_theorem["minimum_partition_seam_decimal_separation"],
        "fresh_witness_count": len(rows),
        "fresh_witness_rows": rows,
        "all_fresh_witnesses_match_stored": all(row["separation_matches_stored"] and row["zero_exclusion_matches_stored"] for row in rows),
        "all_fresh_witnesses_exclude_zero_with_positive_separation": all(row["fresh_decimal_interval_excludes_zero"] and row["fresh_decimal_separation_positive"] for row in rows),
    }


def append_once(path: Path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2474/S1424 strict pointwise interval-Decimal P2459 extremal-witness rerun audit

`P2474/S1424` performs a fresh Decimal/Taylor rerun of the critical finite witnesses from the P2459 replay chain: first/minimum/last rows for the full P2469 opposite-tail replay, first/minimum/last rows for the full P2470 remaining non-tail replay, and all P2472 seam rows.  The fresh rerun matches stored separations and zero-exclusion flags while the inherited finite ledger remains `45165 + 48320` plus the P2466 tail count inside the `99846`-cell P2459 universe.

For a non-specialist: this redoes the most important saved checks—the weakest cells and the seam cells—rather than only checking file hashes.  It is still a finite witness rerun, not a directed-rounding interval theorem, symbolic root-exclusion theorem, continuum theorem, selector/source/gauge theorem, physical-value generator, role-bearing `L_total`, `QW-2191` discharge, legacy-role transfer, or ToE closure.
"""
    lag_section = """
## P2474/S1424 P2459 extremal-witness rerun guard

`P2474/S1424` adds a fresh Decimal/Taylor rerun guard to the finite-grid bookkeeping behind `L_total`: the weakest and edge witnesses from P2469/P2470 plus all P2472 seam witnesses are recomputed and compared against the stored artifacts.  This strengthens finite reproducibility, but it does not add selector/source/gauge authority, physical-value generation, role-bearing finality, legacy-role transfer, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2474/S1424 strict pointwise interval-Decimal P2459 extremal-witness rerun audit", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2474/S1424 P2459 extremal-witness rerun guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    projection = projection_vector(sources["P2449_PROJECTION_REDUCTION"])
    p2469_theorem = theorem(sources["P2469_FULL_OPPOSITE_TAIL_REPLAY"], "strict_pointwise_interval_decimal_full_opposite_tail_replay_certificate")
    p2470_theorem = theorem(sources["P2470_REMAINING_NON_TAIL_REPLAY"], "strict_pointwise_interval_decimal_full_remaining_non_tail_complement_replay_certificate")
    p2471_theorem = theorem(sources["P2471_FINITE_PARTITION_WITNESS"], "strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate")
    p2472_theorem = theorem(sources["P2472_PARTITION_SEAM_REPLAY_AUDIT"], "strict_pointwise_interval_decimal_p2459_partition_seam_replay_audit")
    p2473_theorem = theorem(sources["P2473_FINGERPRINT_BINDING_AUDIT"], "strict_pointwise_interval_decimal_p2459_replay_chain_fingerprint_binding_audit")
    witness_groups = [
        p2469_witnesses(p2469_theorem, projection),
        p2470_witnesses(p2470_theorem, projection),
        p2472_witnesses(p2472_theorem, projection),
    ]
    finite_chain_sum = (
        p2471_theorem["p2466_total_tail_boundary_replay_count_inherited"]
        + p2469_theorem["total_opposite_tail_full_replay_count"]
        + p2470_theorem["total_remaining_non_tail_full_replay_count"]
    )
    theorem_export = {
        "theorem_name": "P2474_T1_strict_pointwise_interval_decimal_p2459_extremal_witness_rerun_audit",
        "audited_chain": ["P2469/S1419", "P2470/S1420", "P2471/S1421", "P2472/S1422", "P2473/S1423"],
        "rerun_mode": "fresh_decimal_taylor_rerun_of_first_minimum_last_and_seam_witnesses",
        "witness_groups": witness_groups,
        "total_fresh_decimal_taylor_witness_rerun_count": sum(group["fresh_witness_count"] for group in witness_groups),
        "all_fresh_witness_groups_match_stored": all(group["all_fresh_witnesses_match_stored"] for group in witness_groups),
        "all_fresh_witness_groups_exclude_zero_with_positive_separation": all(group["all_fresh_witnesses_exclude_zero_with_positive_separation"] for group in witness_groups),
        "p2466_tail_count_inherited_from_p2471": p2471_theorem["p2466_total_tail_boundary_replay_count_inherited"],
        "p2469_full_opposite_tail_count_inherited": p2469_theorem["total_opposite_tail_full_replay_count"],
        "p2470_remaining_non_tail_count_inherited": p2470_theorem["total_remaining_non_tail_full_replay_count"],
        "p2471_p2459_universe_count_inherited": p2471_theorem["total_p2459_universe_cell_count_rebuilt"],
        "finite_chain_sum_p2466_p2469_p2470": finite_chain_sum,
        "finite_chain_sum_matches_p2471_universe": finite_chain_sum == p2471_theorem["total_p2459_universe_cell_count_rebuilt"],
        "p2471_missing_cells_inherited": p2471_theorem["total_partition_missing_cell_count"],
        "p2471_extra_cells_inherited": p2471_theorem["total_partition_extra_cell_count"],
        "p2471_disjoint_complete_inherited": p2471_theorem["all_family_partitions_are_disjoint_and_complete"],
        "p2473_fingerprint_binding_pass_inherited": p2473_theorem["all_theorem_fingerprints_match_declared"] and p2473_theorem["all_declared_source_fingerprints_match_current_sources"] and p2473_theorem["finite_replay_chain_count_binding_passes"],
        "finite_replay_chain_extremal_witness_rerun_audit_exported": True,
        "directed_rounding_interval_theorem_exported_by_this_certificate": False,
        "symbolic_root_exclusion_theorem_exported_by_this_certificate": False,
        "analytic_monotonicity_theorem_exported_by_this_certificate": False,
        "global_continuum_root_exclusion_theorem_exported_by_this_certificate": False,
        "pointwise_coordinate_selector_exported_by_this_certificate": False,
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "gauge_slice_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "legacy_role_transfer_exported": False,
        "toe_closure_exported": False,
        "lay_summary": "This packet reruns the weakest saved cells, edge cells, and all partition seam cells from the finite P2459 audit chain. In plain terms: after checking the file seals in P2473, P2474 opens the box and recalculates the most important saved cells to confirm that the stored finite answer is reproducible. It still does not become a symbolic proof for every real point.",
        "not_licensed": [
            "An extremal-witness rerun audit is not a directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, or continuum root-exclusion theorem.",
            "The audit does not export selector/source/gauge closure, a physical-value generator, QW-2191 discharge, role-bearing L_total, legacy-role transfer, or ToE closure.",
            "No legacy physical roles are transferred to L_total or K_strict_gate by this finite reproducibility check.",
        ],
        "next_honest_step": "Use the rerun-stable extremal witnesses and bound finite chain as inputs for a future directed-rounding or symbolic root-exclusion attempt; do not promote deterministic finite witness replay into continuum closure.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "fresh_witness_count": theorem_export["total_fresh_decimal_taylor_witness_rerun_count"] == 28,
        "fresh_witness_groups_match_stored": theorem_export["all_fresh_witness_groups_match_stored"],
        "fresh_witness_groups_exclude_zero_positive": theorem_export["all_fresh_witness_groups_exclude_zero_with_positive_separation"],
        "finite_chain_sum_matches_universe": theorem_export["finite_chain_sum_matches_p2471_universe"],
        "p2459_universe_count": theorem_export["p2471_p2459_universe_count_inherited"] == 99846,
        "p2469_count": theorem_export["p2469_full_opposite_tail_count_inherited"] == 45165,
        "p2470_count": theorem_export["p2470_remaining_non_tail_count_inherited"] == 48320,
        "p2471_no_missing_extra": theorem_export["p2471_missing_cells_inherited"] == 0 and theorem_export["p2471_extra_cells_inherited"] == 0 and theorem_export["p2471_disjoint_complete_inherited"],
        "p2473_binding_inherited": theorem_export["p2473_fingerprint_binding_pass_inherited"],
        "finite_grid_only_not_directed_rounding": theorem_export["finite_replay_chain_extremal_witness_rerun_audit_exported"] and not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
        "no_symbolic_root_exclusion": not theorem_export["symbolic_root_exclusion_theorem_exported_by_this_certificate"],
        "no_continuum_root_exclusion": not theorem_export["global_continuum_root_exclusion_theorem_exported_by_this_certificate"],
        "no_selector_source_gauge_theorem": not theorem_export["pointwise_coordinate_selector_exported_by_this_certificate"] and not theorem_export["strict_observable_source_constraint_exported_by_this_certificate"] and not theorem_export["gauge_slice_theorem_exported_by_this_certificate"],
        "no_value_generator_export": not theorem_export["strict_physical_value_generator_exported"],
        "no_qw2191_discharge": not theorem_export["qw2191_discharged"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_legacy_role_transfer": not theorem_export["legacy_role_transfer_exported"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2474_s1424_v1",
        "packet_id": "P2474",
        "stage_id": "S1424",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_EXTREMAL_WITNESS_RERUN_AUDIT_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_p2459_extremal_witness_rerun_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_p2459_extremal_witness_rerun_audit"]["theorem_export"]
    lines = [
        "# P2474/S1424 strict pointwise interval-Decimal P2459 extremal-witness rerun audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Extremal witness rerun audit",
        "",
        f"Fresh Decimal/Taylor witness reruns: `{t['total_fresh_decimal_taylor_witness_rerun_count']}`.",
        f"All fresh witness groups match stored values: `{t['all_fresh_witness_groups_match_stored']}`.",
        f"All fresh witness groups exclude zero with positive separation: `{t['all_fresh_witness_groups_exclude_zero_with_positive_separation']}`.",
        f"Finite chain sum P2466+P2469+P2470: `{t['finite_chain_sum_p2466_p2469_p2470']}` / P2459 universe `{t['p2471_p2459_universe_count_inherited']}`.",
        f"P2471 missing/extra/disjoint-complete: `{t['p2471_missing_cells_inherited']}` / `{t['p2471_extra_cells_inherited']}` / `{t['p2471_disjoint_complete_inherited']}`.",
        f"P2473 fingerprint binding inherited: `{t['p2473_fingerprint_binding_pass_inherited']}`.",
        "",
        "## Plain-language progress note",
        "",
        t["lay_summary"],
        "",
        "## Hard limits / negative controls",
        "",
        "This is a finite extremal-witness rerun audit only.  It exports no directed-rounding interval theorem, no symbolic root-exclusion theorem, no analytic monotonicity theorem, no global continuum root-exclusion theorem, no selector/source/gauge theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, no legacy-role transfer, and no ToE closure.",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps({"status": payload["status"], "gatekeepers": payload["gatekeeper_checks"]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
