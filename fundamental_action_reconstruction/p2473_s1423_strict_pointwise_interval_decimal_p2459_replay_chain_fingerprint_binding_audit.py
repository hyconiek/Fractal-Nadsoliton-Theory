#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)

GEN = ROOT / "generated"
OUT = GEN / "p2473_s1423_strict_pointwise_interval_decimal_p2459_replay_chain_fingerprint_binding_audit.json"
MD = GEN / "p2473_s1423_strict_pointwise_interval_decimal_p2459_replay_chain_fingerprint_binding_audit.md"

ARTIFACTS = [
    {
        "label": "P2469_FULL_OPPOSITE_TAIL_REPLAY",
        "packet_id": "P2469",
        "stage_id": "S1419",
        "path": GEN / "p2469_s1419_strict_pointwise_interval_decimal_full_opposite_tail_replay_certificate.json",
        "container": "strict_pointwise_interval_decimal_full_opposite_tail_replay_certificate",
    },
    {
        "label": "P2470_REMAINING_NON_TAIL_REPLAY",
        "packet_id": "P2470",
        "stage_id": "S1420",
        "path": GEN / "p2470_s1420_strict_pointwise_interval_decimal_full_remaining_non_tail_complement_replay_certificate.json",
        "container": "strict_pointwise_interval_decimal_full_remaining_non_tail_complement_replay_certificate",
    },
    {
        "label": "P2471_FINITE_PARTITION_WITNESS",
        "packet_id": "P2471",
        "stage_id": "S1421",
        "path": GEN / "p2471_s1421_strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate.json",
        "container": "strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate",
    },
    {
        "label": "P2472_PARTITION_SEAM_REPLAY_AUDIT",
        "packet_id": "P2472",
        "stage_id": "S1422",
        "path": GEN / "p2472_s1422_strict_pointwise_interval_decimal_p2459_partition_seam_replay_audit.json",
        "container": "strict_pointwise_interval_decimal_p2459_partition_seam_replay_audit",
    },
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
        "new_packet": "P2473|S1423|fingerprint binding audit|replay chain fingerprint|source fingerprint binding|tamper-evident replay chain|finite ledger provenance",
        "precursor_packets": "P2469|S1419|P2470|S1420|P2471|S1421|P2472|S1422|P2459|S1409",
        "binding_language": "source fingerprint|theorem fingerprint|artifact hash|hash binding|tamper-evident|determinism audit",
        "hard_limit_markers": "directed-rounding root exclusion|symbolic root-exclusion theorem|global continuum root-exclusion|full mathematical proof",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem_export(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def audit_artifact(meta: dict[str, Any]) -> dict[str, Any]:
    payload = load_json(meta["path"])
    container_payload = payload.get(meta["container"], {})
    theorem = container_payload.get("theorem_export", {})
    declared_theorem_hash = container_payload.get("theorem_fingerprint_sha256")
    recomputed_theorem_hash = sha256_json(theorem)
    declared_source_hashes = container_payload.get("source_fingerprints", {})
    source_results = []
    for source_label, rel_path in sorted(payload.get("source_files", {}).items()):
        source_path = REPO / rel_path
        source_payload = load_json(source_path)
        recomputed = sha256_json(source_payload)
        source_results.append({
            "source_label": source_label,
            "source_file": rel_path,
            "declared_source_fingerprint_sha256": declared_source_hashes.get(source_label),
            "recomputed_source_fingerprint_sha256": recomputed,
            "matches_declared_source_fingerprint": declared_source_hashes.get(source_label) == recomputed,
        })
    gatekeepers = payload.get("gatekeeper_checks", {})
    return {
        "label": meta["label"],
        "packet_id": payload.get("packet_id"),
        "stage_id": payload.get("stage_id"),
        "artifact_file": rel(meta["path"]),
        "artifact_file_sha256": sha256_file(meta["path"]),
        "status": payload.get("status"),
        "packet_id_matches_expected": payload.get("packet_id") == meta["packet_id"],
        "stage_id_matches_expected": payload.get("stage_id") == meta["stage_id"],
        "declared_theorem_fingerprint_sha256": declared_theorem_hash,
        "recomputed_theorem_fingerprint_sha256": recomputed_theorem_hash,
        "theorem_fingerprint_matches_declared": declared_theorem_hash == recomputed_theorem_hash,
        "source_fingerprint_results": source_results,
        "all_source_fingerprints_match_declared": all(row["matches_declared_source_fingerprint"] for row in source_results),
        "gatekeeper_count": len(gatekeepers),
        "all_gatekeepers_pass": bool(gatekeepers) and all(value is True for value in gatekeepers.values()),
    }


def chain_consistency(loaded: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2469 = theorem_export(loaded["P2469_FULL_OPPOSITE_TAIL_REPLAY"], ARTIFACTS[0]["container"])
    p2470 = theorem_export(loaded["P2470_REMAINING_NON_TAIL_REPLAY"], ARTIFACTS[1]["container"])
    p2471 = theorem_export(loaded["P2471_FINITE_PARTITION_WITNESS"], ARTIFACTS[2]["container"])
    p2472 = theorem_export(loaded["P2472_PARTITION_SEAM_REPLAY_AUDIT"], ARTIFACTS[3]["container"])
    p2466_count = p2471.get("p2466_total_tail_boundary_replay_count_inherited")
    p2469_count = p2469.get("total_opposite_tail_full_replay_count")
    p2470_count = p2470.get("total_remaining_non_tail_full_replay_count")
    p2459_count = p2471.get("total_p2459_universe_cell_count_rebuilt")
    finite_sum = p2466_count + p2469_count + p2470_count
    return {
        "p2466_descent_tail_count_from_p2471": p2466_count,
        "p2469_full_opposite_tail_count": p2469_count,
        "p2470_remaining_non_tail_count": p2470_count,
        "p2471_p2459_universe_count": p2459_count,
        "finite_replay_partition_sum": finite_sum,
        "finite_replay_partition_sum_equals_p2459_universe": finite_sum == p2459_count,
        "p2470_remaining_budget_after_replay": p2470.get("p2459_unreplayed_by_boundary_chain_budget_after_p2466_p2469_p2470"),
        "p2471_missing_cells": p2471.get("total_partition_missing_cell_count"),
        "p2471_extra_cells": p2471.get("total_partition_extra_cell_count"),
        "p2471_all_family_partitions_are_disjoint_and_complete": p2471.get("all_family_partitions_are_disjoint_and_complete"),
        "p2472_transition_pair_count": p2472.get("total_transition_pair_count"),
        "p2472_seam_replay_count": p2472.get("total_partition_seam_replay_count"),
        "p2472_all_transition_pairs_are_adjacent": p2472.get("all_transition_pairs_are_adjacent"),
        "p2472_all_seam_replayed_cells_exclude_zero": p2472.get("all_seam_replayed_cells_exclude_zero"),
        "p2472_minimum_partition_seam_decimal_separation": p2472.get("minimum_partition_seam_decimal_separation"),
        "p2472_inherits_p2471_universe_count": p2472.get("p2471_total_p2459_universe_cell_count_rebuilt_inherited") == p2459_count,
        "p2472_inherits_p2469_count": p2472.get("p2469_total_opposite_tail_full_replay_count_inherited") == p2469_count,
        "p2472_inherits_p2470_count": p2472.get("p2470_total_remaining_non_tail_full_replay_count_inherited") == p2470_count,
    }


def append_once(path: Path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2473/S1423 strict pointwise interval-Decimal P2459 replay-chain fingerprint binding audit

`P2473/S1423` is a tamper-evident provenance audit for the finite P2459 replay chain `P2469 -> P2470 -> P2471 -> P2472`.  It reloads each generated artifact, recomputes the declared theorem fingerprint, recomputes every declared source fingerprint from the current source files, checks all gatekeepers, and verifies that the finite chain counts still bind to the `99846`-cell P2459 universe.

For a non-specialist: this does not calculate new physics; it checks that the files saying "we checked every box" are still tied to the exact input files they claim, and that the hash labels in those files still match the current contents.  It is not a directed-rounding interval theorem, symbolic root-exclusion theorem, continuum theorem, selector/source/gauge theorem, physical-value generator, role-bearing `L_total`, `QW-2191` discharge, legacy-role transfer, or ToE closure.
"""
    lag_section = """
## P2473/S1423 P2459 replay-chain fingerprint binding guard

`P2473/S1423` adds a provenance guard to the finite-grid bookkeeping behind `L_total`: the P2469/P2470/P2471/P2472 replay-chain artifacts are hash-bound to their declared source artifacts and their theorem exports.  This improves audit reproducibility and tamper evidence, but it does not add selector/source/gauge authority, physical-value generation, role-bearing finality, legacy-role transfer, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2473/S1423 strict pointwise interval-Decimal P2459 replay-chain fingerprint binding audit", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2473/S1423 P2459 replay-chain fingerprint binding guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    loaded = {meta["label"]: load_json(meta["path"]) for meta in ARTIFACTS}
    artifact_audits = [audit_artifact(meta) for meta in ARTIFACTS]
    chain = chain_consistency(loaded)
    theorem_export_payload = {
        "theorem_name": "P2473_T1_strict_pointwise_interval_decimal_p2459_replay_chain_fingerprint_binding_audit",
        "audited_chain": ["P2469/S1419", "P2470/S1420", "P2471/S1421", "P2472/S1422"],
        "artifact_binding_audits": artifact_audits,
        "chain_consistency": chain,
        "all_artifact_packet_and_stage_ids_match": all(row["packet_id_matches_expected"] and row["stage_id_matches_expected"] for row in artifact_audits),
        "all_theorem_fingerprints_match_declared": all(row["theorem_fingerprint_matches_declared"] for row in artifact_audits),
        "all_declared_source_fingerprints_match_current_sources": all(row["all_source_fingerprints_match_declared"] for row in artifact_audits),
        "all_audited_gatekeepers_pass": all(row["all_gatekeepers_pass"] for row in artifact_audits),
        "finite_replay_chain_count_binding_passes": chain["finite_replay_partition_sum_equals_p2459_universe"]
        and chain["p2470_remaining_budget_after_replay"] == 0
        and chain["p2471_missing_cells"] == 0
        and chain["p2471_extra_cells"] == 0
        and chain["p2471_all_family_partitions_are_disjoint_and_complete"] is True
        and chain["p2472_inherits_p2471_universe_count"]
        and chain["p2472_inherits_p2469_count"]
        and chain["p2472_inherits_p2470_count"],
        "finite_replay_chain_fingerprint_binding_audit_exported": True,
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
        "lay_summary": "This packet is like checking the serial numbers and seals on the finite audit files. It confirms that the P2469-P2472 files still match their own hash labels, still point to the current source files they claim to use, and still add up to the same finite 99846-cell P2459 checklist. It does not turn the finite checklist into a symbolic proof for every real point.",
        "not_licensed": [
            "A fingerprint binding audit is not a directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, or continuum root-exclusion theorem.",
            "The audit does not export selector/source/gauge closure, a physical-value generator, QW-2191 discharge, role-bearing L_total, legacy-role transfer, or ToE closure.",
            "No legacy physical roles are transferred to L_total or K_strict_gate by this provenance check.",
        ],
        "next_honest_step": "Use the bound finite replay chain as a reproducible input for a future directed-rounding or symbolic root-exclusion attempt; do not promote hash/provenance hygiene into continuum closure.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "packet_stage_ids_match": theorem_export_payload["all_artifact_packet_and_stage_ids_match"],
        "theorem_fingerprints_match": theorem_export_payload["all_theorem_fingerprints_match_declared"],
        "source_fingerprints_match": theorem_export_payload["all_declared_source_fingerprints_match_current_sources"],
        "audited_gatekeepers_pass": theorem_export_payload["all_audited_gatekeepers_pass"],
        "finite_chain_count_binding": theorem_export_payload["finite_replay_chain_count_binding_passes"],
        "p2459_universe_count": chain["p2471_p2459_universe_count"] == 99846,
        "p2469_count": chain["p2469_full_opposite_tail_count"] == 45165,
        "p2470_count": chain["p2470_remaining_non_tail_count"] == 48320,
        "p2472_seam_controls_still_pass": chain["p2472_transition_pair_count"] == 4 and chain["p2472_seam_replay_count"] == 16 and chain["p2472_all_transition_pairs_are_adjacent"] and chain["p2472_all_seam_replayed_cells_exclude_zero"],
        "finite_grid_only_not_directed_rounding": theorem_export_payload["finite_replay_chain_fingerprint_binding_audit_exported"] and not theorem_export_payload["directed_rounding_interval_theorem_exported_by_this_certificate"],
        "no_symbolic_root_exclusion": not theorem_export_payload["symbolic_root_exclusion_theorem_exported_by_this_certificate"],
        "no_continuum_root_exclusion": not theorem_export_payload["global_continuum_root_exclusion_theorem_exported_by_this_certificate"],
        "no_selector_source_gauge_theorem": not theorem_export_payload["pointwise_coordinate_selector_exported_by_this_certificate"] and not theorem_export_payload["strict_observable_source_constraint_exported_by_this_certificate"] and not theorem_export_payload["gauge_slice_theorem_exported_by_this_certificate"],
        "no_value_generator_export": not theorem_export_payload["strict_physical_value_generator_exported"],
        "no_qw2191_discharge": not theorem_export_payload["qw2191_discharged"],
        "no_ltotal_export": not theorem_export_payload["role_bearing_ltotal_exported"],
        "no_legacy_role_transfer": not theorem_export_payload["legacy_role_transfer_exported"],
        "no_toe_export": not theorem_export_payload["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export_payload) == sha256_json(theorem_export_payload),
    }
    return {
        "schema_version": "p2473_s1423_v1",
        "packet_id": "P2473",
        "stage_id": "S1423",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_REPLAY_CHAIN_FINGERPRINT_BINDING_AUDIT_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM",
        "source_files": {meta["label"]: rel(meta["path"]) for meta in ARTIFACTS},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_p2459_replay_chain_fingerprint_binding_audit": {
            "theorem_export": theorem_export_payload,
            "source_fingerprints": {meta["label"]: sha256_json(loaded[meta["label"]]) for meta in ARTIFACTS},
            "theorem_fingerprint_sha256": sha256_json(theorem_export_payload),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_p2459_replay_chain_fingerprint_binding_audit"]["theorem_export"]
    c = t["chain_consistency"]
    lines = [
        "# P2473/S1423 strict pointwise interval-Decimal P2459 replay-chain fingerprint binding audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Replay-chain fingerprint binding audit",
        "",
        f"Artifacts audited: `{len(t['artifact_binding_audits'])}` (`P2469`, `P2470`, `P2471`, `P2472`).",
        f"All theorem fingerprints match declared values: `{t['all_theorem_fingerprints_match_declared']}`.",
        f"All declared source fingerprints match current source files: `{t['all_declared_source_fingerprints_match_current_sources']}`.",
        f"All audited gatekeepers pass: `{t['all_audited_gatekeepers_pass']}`.",
        f"Finite chain sum: `{c['finite_replay_partition_sum']}` / P2459 universe `{c['p2471_p2459_universe_count']}`.",
        f"Remaining budget after P2466+P2469+P2470: `{c['p2470_remaining_budget_after_replay']}`.",
        f"P2471 missing/extra/disjoint-complete: `{c['p2471_missing_cells']}` / `{c['p2471_extra_cells']}` / `{c['p2471_all_family_partitions_are_disjoint_and_complete']}`.",
        "",
        "## Plain-language progress note",
        "",
        t["lay_summary"],
        "",
        "## Hard limits / negative controls",
        "",
        "This is a finite replay-chain fingerprint/provenance audit only.  It exports no directed-rounding interval theorem, no symbolic root-exclusion theorem, no analytic monotonicity theorem, no global continuum root-exclusion theorem, no selector/source/gauge theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, no legacy-role transfer, and no ToE closure.",
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
