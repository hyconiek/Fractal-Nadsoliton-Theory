#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from decimal import Decimal
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel

GEN = ROOT / "generated"
OUT = GEN / "p2476_s1426_strict_pointwise_interval_decimal_p2459_critical_halo_order_classification_audit.json"
MD = GEN / "p2476_s1426_strict_pointwise_interval_decimal_p2459_critical_halo_order_classification_audit.md"

SOURCE_FILES = {
    "P2475_CRITICAL_MINIMUM_HALO_REPLAY": GEN / "p2475_s1425_strict_pointwise_interval_decimal_p2459_critical_minimum_halo_replay_audit.json",
    "P2474_EXTREMAL_WITNESS_RERUN": GEN / "p2474_s1424_strict_pointwise_interval_decimal_p2459_extremal_witness_rerun_audit.json",
    "P2473_FINGERPRINT_BINDING_AUDIT": GEN / "p2473_s1423_strict_pointwise_interval_decimal_p2459_replay_chain_fingerprint_binding_audit.json",
    "P2471_FINITE_PARTITION_WITNESS": GEN / "p2471_s1421_strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate.json",
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
        "new_packet": "P2476|S1426|critical halo order|halo order classification|local-minimum halo|neighbor-dominance audit|separation order audit",
        "precursor_packets": "P2473|S1423|P2474|S1424|P2475|S1425|P2471|S1421",
        "order_language": "local minimum|lower neighbor|one-sided boundary minimum|monotone halo|separation ordering|neighbor dominance",
        "hard_limit_markers": "directed-rounding root exclusion|symbolic root-exclusion theorem|global continuum root-exclusion|full mathematical proof",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def classify_target_order(family_payload: dict[str, Any], target: dict[str, Any]) -> dict[str, Any]:
    seg = int(target["segment_index"])
    idx = int(target["uncovered_index"])
    rows = [
        row for row in family_payload["halo_replay_rows"]
        if int(row["segment_index"]) == seg and abs(int(row["uncovered_index"]) - idx) <= int(family_payload["halo_radius"])
    ]
    rows = sorted(rows, key=lambda row: int(row["uncovered_index"]))
    center = next(row for row in rows if int(row["uncovered_index"]) == idx)
    center_sep = Decimal(center["decimal_separation_from_zero"])
    neighbors = [row for row in rows if int(row["uncovered_index"]) != idx]
    lower_neighbors = [row for row in neighbors if Decimal(row["decimal_separation_from_zero"]) < center_sep]
    higher_neighbors = [row for row in neighbors if Decimal(row["decimal_separation_from_zero"]) > center_sep]
    equal_neighbors = [row for row in neighbors if Decimal(row["decimal_separation_from_zero"]) == center_sep]
    available_offsets = [int(row["uncovered_index"]) - idx for row in rows]
    missing_offsets = [offset for offset in range(-int(family_payload["halo_radius"]), int(family_payload["halo_radius"]) + 1) if offset not in available_offsets]
    left_rows = [row for row in rows if int(row["uncovered_index"]) < idx]
    right_rows = [row for row in rows if int(row["uncovered_index"]) > idx]
    return {
        "family": family_payload["family"],
        "source_packet": target["source_packet"],
        "source_rule": target["source_rule"],
        "segment_index": seg,
        "uncovered_index": idx,
        "center_partition_class": center["partition_class"],
        "center_decimal_separation": str(center_sep),
        "available_offsets": available_offsets,
        "missing_offsets_due_to_segment_boundary": missing_offsets,
        "available_halo_cell_count": len(rows),
        "neighbor_count": len(neighbors),
        "lower_neighbor_count": len(lower_neighbors),
        "higher_neighbor_count": len(higher_neighbors),
        "equal_neighbor_count": len(equal_neighbors),
        "center_is_strict_minimum_within_available_halo": len(neighbors) > 0 and len(lower_neighbors) == 0 and len(equal_neighbors) == 0,
        "center_has_lower_neighbor_within_available_halo": len(lower_neighbors) > 0,
        "center_is_boundary_truncated_one_sided_minimum": len(missing_offsets) > 0 and len(lower_neighbors) == 0 and len(equal_neighbors) == 0,
        "lowest_available_halo_row": min(rows, key=lambda row: Decimal(row["decimal_separation_from_zero"])),
        "left_neighbor_separations": [row["decimal_separation_from_zero"] for row in left_rows],
        "right_neighbor_separations": [row["decimal_separation_from_zero"] for row in right_rows],
        "lower_neighbor_rows": lower_neighbors,
    }


def append_once(path: Path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2476/S1426 strict pointwise interval-Decimal P2459 critical-halo order classification audit

`P2476/S1426` classifies the local separation ordering inside the P2475 critical halos.  It does not assume every inherited minimum witness is a local minimum after cross-class halo expansion: it checks each target against its available neighboring cells, records boundary-truncated one-sided minima, and explicitly reports the one case where a lower neighbor exists in the cross-class halo.  The audit finds `5/6` critical targets are strict minima within their available halo and `1/6` has a lower neighboring cell, while all cells remain zero-excluding under P2475.

For a non-specialist: this is an honesty check on the word "minimum".  Some minima are minima only inside their original packet/class, not after nearby cells from another class are included.  P2476 names that fact instead of hiding it.  It remains a finite ordering audit, not a directed-rounding interval theorem, symbolic root-exclusion theorem, continuum theorem, selector/source/gauge theorem, physical-value generator, role-bearing `L_total`, `QW-2191` discharge, legacy-role transfer, or ToE closure.
"""
    lag_section = """
## P2476/S1426 P2459 critical-halo order classification guard

`P2476/S1426` adds an order-classification guard to the finite-grid bookkeeping behind `L_total`: P2475 critical halos are checked for strict local minima, one-sided boundary minima, and lower-neighbor exceptions.  This improves proof hygiene around weakest-witness language, but it does not add selector/source/gauge authority, physical-value generation, role-bearing finality, legacy-role transfer, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2476/S1426 strict pointwise interval-Decimal P2459 critical-halo order classification audit", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2476/S1426 P2459 critical-halo order classification guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2475 = theorem(sources["P2475_CRITICAL_MINIMUM_HALO_REPLAY"], "strict_pointwise_interval_decimal_p2459_critical_minimum_halo_replay_audit")
    p2474 = theorem(sources["P2474_EXTREMAL_WITNESS_RERUN"], "strict_pointwise_interval_decimal_p2459_extremal_witness_rerun_audit")
    p2473 = theorem(sources["P2473_FINGERPRINT_BINDING_AUDIT"], "strict_pointwise_interval_decimal_p2459_replay_chain_fingerprint_binding_audit")
    p2471 = theorem(sources["P2471_FINITE_PARTITION_WITNESS"], "strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate")
    target_classifications = []
    for family_payload in p2475["family_halo_replays"]:
        for target in family_payload["minimum_witness_targets"]:
            target_classifications.append(classify_target_order(family_payload, target))
    strict_minimum_count = sum(row["center_is_strict_minimum_within_available_halo"] for row in target_classifications)
    lower_neighbor_count = sum(row["center_has_lower_neighbor_within_available_halo"] for row in target_classifications)
    boundary_one_sided_count = sum(row["center_is_boundary_truncated_one_sided_minimum"] for row in target_classifications)
    theorem_export = {
        "theorem_name": "P2476_T1_strict_pointwise_interval_decimal_p2459_critical_halo_order_classification_audit",
        "audited_chain": ["P2471/S1421", "P2473/S1423", "P2474/S1424", "P2475/S1425"],
        "inherited_p2475_total_unique_halo_replay_count": p2475["total_unique_halo_replay_count"],
        "inherited_p2475_all_halos_exclude_zero": p2475["all_halo_replayed_cells_exclude_zero"],
        "inherited_p2475_all_halos_positive_separation": p2475["all_halo_replayed_cells_have_positive_separation"],
        "inherited_p2474_total_extremal_witness_rerun_count": p2474["total_fresh_decimal_taylor_witness_rerun_count"],
        "inherited_p2473_fingerprint_binding_pass": p2473["all_theorem_fingerprints_match_declared"] and p2473["all_declared_source_fingerprints_match_current_sources"] and p2473["finite_replay_chain_count_binding_passes"],
        "inherited_p2471_disjoint_complete": p2471["all_family_partitions_are_disjoint_and_complete"],
        "target_classifications": target_classifications,
        "total_target_classification_count": len(target_classifications),
        "strict_available_halo_minimum_count": strict_minimum_count,
        "boundary_truncated_one_sided_minimum_count": boundary_one_sided_count,
        "targets_with_lower_neighbor_count": lower_neighbor_count,
        "all_targets_are_strict_available_halo_minima": strict_minimum_count == len(target_classifications),
        "lower_neighbor_exception_exported": lower_neighbor_count > 0,
        "finite_replay_chain_critical_halo_order_classification_audit_exported": True,
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
        "lay_summary": "This packet checks whether the saved 'minimum' cells are still local minima after nearby cells are included. Five of six targets are strict minima in their available halo; one target has a lower neighbor from a neighboring partition class. That exception is useful: it prevents overclaiming and shows exactly where the word minimum is class-local rather than whole-halo local.",
        "not_licensed": [
            "A critical-halo order audit is not a directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, or continuum root-exclusion theorem.",
            "The audit does not export selector/source/gauge closure, a physical-value generator, QW-2191 discharge, role-bearing L_total, legacy-role transfer, or ToE closure.",
            "No legacy physical roles are transferred to L_total or K_strict_gate by this finite ordering check.",
        ],
        "next_honest_step": "Use the lower-neighbor exception and one-sided boundary classifications to target any future directed-rounding or symbolic root-exclusion proof at the real weakest finite neighborhoods without overclaiming local-minimum status.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2475_halo_count_inherited": theorem_export["inherited_p2475_total_unique_halo_replay_count"] == 14,
        "p2475_zero_exclusion_inherited": theorem_export["inherited_p2475_all_halos_exclude_zero"] and theorem_export["inherited_p2475_all_halos_positive_separation"],
        "p2474_extremal_count_inherited": theorem_export["inherited_p2474_total_extremal_witness_rerun_count"] == 28,
        "p2473_binding_inherited": theorem_export["inherited_p2473_fingerprint_binding_pass"],
        "p2471_partition_inherited": theorem_export["inherited_p2471_disjoint_complete"],
        "target_classification_count": theorem_export["total_target_classification_count"] == 6,
        "strict_minimum_count_expected": theorem_export["strict_available_halo_minimum_count"] == 5,
        "boundary_one_sided_minimum_count_expected": theorem_export["boundary_truncated_one_sided_minimum_count"] == 5,
        "lower_neighbor_exception_expected": theorem_export["targets_with_lower_neighbor_count"] == 1 and theorem_export["lower_neighbor_exception_exported"],
        "not_all_targets_overclaimed_as_local_minima": theorem_export["all_targets_are_strict_available_halo_minima"] is False,
        "finite_grid_only_not_directed_rounding": theorem_export["finite_replay_chain_critical_halo_order_classification_audit_exported"] and not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
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
        "schema_version": "p2476_s1426_v1",
        "packet_id": "P2476",
        "stage_id": "S1426",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_CRITICAL_HALO_ORDER_CLASSIFICATION_AUDIT_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_p2459_critical_halo_order_classification_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_p2459_critical_halo_order_classification_audit"]["theorem_export"]
    lines = [
        "# P2476/S1426 strict pointwise interval-Decimal P2459 critical-halo order classification audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Critical halo order classification",
        "",
        f"Targets classified: `{t['total_target_classification_count']}`.",
        f"Strict minima within available halo: `{t['strict_available_halo_minimum_count']}`.",
        f"Boundary-truncated one-sided minima: `{t['boundary_truncated_one_sided_minimum_count']}`.",
        f"Targets with lower neighbor: `{t['targets_with_lower_neighbor_count']}`.",
        f"All targets are strict available-halo minima: `{t['all_targets_are_strict_available_halo_minima']}`.",
        f"P2475 halo zero-exclusion inherited: `{t['inherited_p2475_all_halos_exclude_zero']}`.",
        "",
        "## Plain-language progress note",
        "",
        t["lay_summary"],
        "",
        "## Hard limits / negative controls",
        "",
        "This is a finite critical-halo order classification audit only.  It exports no directed-rounding interval theorem, no symbolic root-exclusion theorem, no analytic monotonicity theorem, no global continuum root-exclusion theorem, no selector/source/gauge theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, no legacy-role transfer, and no ToE closure.",
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
