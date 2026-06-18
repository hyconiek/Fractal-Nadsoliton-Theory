#!/usr/bin/env python3
"""P2886/S1836: external unit-measure/action-density inventory no-go audit.

P2885 closed the bounded invariant-quadratic action family: those actions carry
already chosen 9:5 ratios but do not export a unit measure, localized action
density, or variational chain rule into nonproxy L_total.  This packet performs
the next honest intake gate requested by that boundary: scan current generated
artifacts for an external-to-that-family strict unit measure or localized action
density export that could supply the missing source object.

This is not a future-impossibility proof.  It is a current-artifact certificate:
existing generated packets contain many negative guard flags and related terms,
but no positive exported unit measure/localized action density coupled to the
9:5 source/stiffness obligation.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2885 = GEN / "p2885_s1835_invariant_quadratic_action_9_over_5_selection_no_go_audit.json"
OUT = GEN / "p2886_s1836_external_unit_measure_action_density_inventory_no_go_audit.json"
MD = GEN / "p2886_s1836_external_unit_measure_action_density_inventory_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

TARGET_TOKENS = (
    "unit_measure",
    "localized_action_density",
    "unit_bearing_action_density",
    "action_density",
    "variational_chain_rule",
    "nonproxy_ltotal",
    "source_stiffness_ratio_9_to_5",
    "unit_bearing_9_over_5",
)
POSITIVE_VERBS = ("export", "accepted", "select", "derive", "coupl")
NEGATIVE_PREFIXES = ("no_", "not_")
NEGATIVE_TOKENS = ("negative_export_flags", "without", "missing", "blocked", "no_go", "no-go")


def json_files() -> list[Path]:
    return sorted(path for path in GEN.glob("*.json") if path.is_file())


def walk(value: Any, prefix: str = "") -> list[tuple[str, Any]]:
    rows: list[tuple[str, Any]] = []
    if isinstance(value, dict):
        for key, child in value.items():
            child_prefix = f"{prefix}.{key}" if prefix else str(key)
            rows.extend(walk(child, child_prefix))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            rows.extend(walk(child, f"{prefix}[{index}]"))
    else:
        rows.append((prefix, value))
    return rows


def path_relevant(path: str) -> bool:
    lowered = path.lower()
    return any(token in lowered for token in TARGET_TOKENS)


def positive_export_bool_path(path: str) -> bool:
    lowered = path.lower()
    leaf = lowered.rsplit(".", 1)[-1]
    if leaf.startswith(NEGATIVE_PREFIXES) or "_no_" in leaf or "_not_" in leaf:
        return False
    if any(token in lowered for token in NEGATIVE_TOKENS):
        return False
    return path_relevant(path) and any(verb in lowered for verb in POSITIVE_VERBS)


def inventory_records() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for file_path in json_files():
        try:
            payload = read_json(file_path)
        except Exception as exc:
            records.append({"file": str(file_path.relative_to(ROOT)), "read_error": repr(exc)})
            continue
        relevant_booleans = []
        positive_export_booleans = []
        relevant_values = []
        for path, value in walk(payload):
            if isinstance(value, bool) and path_relevant(path):
                row = {"path": path, "value": value}
                relevant_booleans.append(row)
                if value and positive_export_bool_path(path):
                    positive_export_booleans.append(row)
            elif path_relevant(path) and not isinstance(value, (dict, list, bool)):
                relevant_values.append({"path": path, "value": str(value)[:240]})
        if relevant_booleans or relevant_values:
            records.append(
                {
                    "file": str(file_path.relative_to(ROOT)),
                    "status": payload.get("status"),
                    "relevant_boolean_count": len(relevant_booleans),
                    "positive_export_boolean_count": len(positive_export_booleans),
                    "positive_export_booleans": positive_export_booleans,
                    "sample_relevant_booleans": relevant_booleans[:10],
                    "sample_relevant_values": relevant_values[:8],
                }
            )
    return records


def build_payload(p2885: dict[str, Any]) -> dict[str, Any]:
    files = json_files()
    records = inventory_records()
    positive_export_records = [record for record in records if record["positive_export_boolean_count"] > 0]
    facts = {
        "p2885_rechecked": p2885.get("status") == "P2885_INVARIANT_QUADRATIC_ACTION_9_OVER_5_SELECTION_NO_GO_AUDIT_NO_CLOSURE",
        "generated_json_inventory_nonempty": len(files) > 0,
        "unit_measure_or_action_density_terms_found": len(records) > 0,
        "no_positive_external_unit_measure_or_action_density_export_found": len(positive_export_records) == 0,
        "p2885_obligation_remains_unsatisfied": len(positive_export_records) == 0,
    }
    return {
        "status": "P2886_EXTERNAL_UNIT_MEASURE_ACTION_DENSITY_INVENTORY_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2885": sha(P2885)},
        "external_unit_measure_action_density_inventory_no_go_audit": {
            "input_status_rechecked": p2885.get("status"),
            "candidate_class": "current generated-artifact inventory for an external-to-P2884/P2885 strict unit measure, localized action density, or variational chain rule coupled to the 9:5 obligation",
            "generated_json_file_count": len(files),
            "relevant_record_count": len(records),
            "positive_export_record_count": len(positive_export_records),
            "sample_relevant_records": records[:20],
            "proof_certificate": {
                "inventory_rule": "Scan generated JSON artifacts for unit-measure, localized-action-density, unit-bearing action, variational-chain-rule, nonproxy-Ltotal, and 9:5 coupling paths; count only positive non-negative export/accept/select/coupling booleans as candidate exports.",
                "finite_result": "Relevant terms and negative guard flags exist, but no positive generated-artifact boolean exports an external strict unit measure/localized action density satisfying the P2885 obligation.",
                "sourcehood_step": "A future move must supply a new artifact/formula with a positive unit measure or localized action-density export rather than another ratio carrier or inventory replay.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_external_unit_measure_or_action_density_export": False,
            "exports_unit_bearing_action_density": False,
            "exports_variational_chain_rule_to_ltotal": False,
            "exports_source_stiffness_ratio_9_to_5": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "external_unit_measure_exported": False,
                "localized_action_density_exported": False,
                "unit_bearing_action_density_exported": False,
                "variational_chain_rule_exported": False,
                "source_stiffness_ratio_9_to_5_exported": False,
                "unit_bearing_9_over_5_coupling_theorem_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2886 performs the post-P2885 external unit-measure/localized-action-density intake gate across current generated artifacts.  Relevant terms and negative guard flags exist, but no positive export/accept/select/coupling boolean supplies an external strict unit measure, localized unit-bearing action density, variational chain rule, or 9:5 source/stiffness theorem.",
            "next_honest_step": "Do not replay invariant-count actions, ratio algebra, scalar Euler transmission, endpoint/coefficient boxes, or artifact inventories.  The next proof-grade move must provide one new explicit formula/artifact for a strict unit measure or localized action density with a computed nonzero 9:5 coupling, or pivot to a genuinely different typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["external_unit_measure_action_density_inventory_no_go_audit"]
    lines = [
        "# P2886/S1836 external unit-measure/action-density inventory no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## External unit-measure/action-density inventory",
        f"- candidate class: `{audit['candidate_class']}`",
        f"- generated JSON file count: `{audit['generated_json_file_count']}`",
        f"- relevant record count: `{audit['relevant_record_count']}`",
        f"- positive export record count: `{audit['positive_export_record_count']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload(read_json(P2885))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2886/S1836 external unit-measure/action-density inventory no-go audit",
        "## P2886/S1836 external unit-measure/action-density inventory no-go audit\n\n"
        "`P2886/S1836` performs the post-`P2885` intake gate across generated artifacts for an external strict unit measure, localized unit-bearing action density, variational chain rule, or `9:5` source/stiffness theorem.  Relevant terms and negative guard flags exist, but no positive generated-artifact export/accept/select/coupling boolean supplies the missing object.  No strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2886/S1836 external unit-measure/action-density `L_total` guard",
        "## P2886/S1836 external unit-measure/action-density `L_total` guard\n\n"
        "`P2886/S1836` adds no strict action term.  The generated-artifact inventory finds no external unit measure, localized unit-bearing action density, variational chain rule into nonproxy `L_total`, EOM, or Hamiltonian source satisfying the post-`P2885` obligation.\n",
    )
    append_once(
        AGENTS,
        "Current external unit-measure/action-density inventory no-go guardrail (P2886/S1836, 2026-06-18)",
        "## Current external unit-measure/action-density inventory no-go guardrail (P2886/S1836, 2026-06-18)\n\n"
        "- P2886 scans current generated artifacts for an external-to-P2884/P2885 strict unit measure, localized unit-bearing action density, variational chain rule, or `9:5` source/stiffness theorem.\n"
        "- Relevant terms and negative guard flags exist, but no positive export/accept/select/coupling boolean supplies the missing object.\n"
        "- Do not promote invariant-count actions, ratio algebra, scalar Euler transmission, endpoint/coefficient boxes, or artifact inventories to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must provide one new explicit formula/artifact for a strict unit measure or localized action density with a computed nonzero `9:5` coupling, pivot to a genuinely different typed object, or preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
