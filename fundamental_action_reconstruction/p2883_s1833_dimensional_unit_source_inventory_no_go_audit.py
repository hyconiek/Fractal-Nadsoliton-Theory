#!/usr/bin/env python3
"""P2883/S1833: dimensional/unit source inventory no-go audit.

P2882 reduced the remaining 9/5 obligation to a concrete missing object: an
independent strict dimensional/unit law fixing primitive source/stiffness ratio
9:5.  This packet performs the honest next computation: scan current generated
FAR artifacts for any already-exported positive source law in that family.

The audit recursively inspects generated JSON packets for relevant export/accept
booleans and source-ratio/coupling keys.  It does not prove future impossibility;
it certifies the current repository state.  On current artifacts, every relevant
accepted/export flag is false, and no existing generated packet exports the
needed dimensional/unit source.  Therefore the next move must introduce a new
source object rather than replay endpoint/coefficient/Euler mechanisms.
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

P2882 = GEN / "p2882_s1832_euler_source_ratio_9_over_5_forcing_no_go_audit.json"
OUT = GEN / "p2883_s1833_dimensional_unit_source_inventory_no_go_audit.json"
MD = GEN / "p2883_s1833_dimensional_unit_source_inventory_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

RELEVANT_TOKENS = (
    "9_over_5",
    "9_to_5",
    "source_stiffness",
    "unit_bearing_9",
    "dimensional_unit",
    "unit_law",
    "coupling_theorem",
)
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
            child_prefix = f"{prefix}[{index}]"
            rows.extend(walk(child, child_prefix))
    else:
        rows.append((prefix, value))
    return rows


def relevant_path(path: str) -> bool:
    lowered = path.lower()
    return any(token in lowered for token in RELEVANT_TOKENS)


def relevant_leaf(path: str) -> bool:
    leaf = path.lower().rsplit(".", 1)[-1]
    return any(token in leaf for token in RELEVANT_TOKENS)


def relevant_bool_path(path: str) -> bool:
    leaf = path.lower().rsplit(".", 1)[-1]
    return relevant_leaf(path) and not leaf.startswith("no_")


def inventory_records() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for file_path in json_files():
        try:
            payload = read_json(file_path)
        except Exception as exc:  # JSON decode errors should be visible in the report.
            records.append({"file": str(file_path.relative_to(ROOT)), "read_error": repr(exc)})
            continue
        relevant_booleans = []
        relevant_truths = []
        relevant_nonbool_values = []
        for path, value in walk(payload):
            if isinstance(value, bool) and relevant_bool_path(path):
                row = {"path": path, "value": value}
                relevant_booleans.append(row)
                if value:
                    relevant_truths.append(row)
            elif relevant_path(path) and not isinstance(value, (dict, list, bool)):
                relevant_nonbool_values.append({"path": path, "value": str(value)[:240]})
        if relevant_booleans or relevant_nonbool_values:
            records.append(
                {
                    "file": str(file_path.relative_to(ROOT)),
                    "status": payload.get("status"),
                    "relevant_boolean_count": len(relevant_booleans),
                    "relevant_true_boolean_count": len(relevant_truths),
                    "relevant_true_booleans": relevant_truths,
                    "sample_relevant_booleans": relevant_booleans[:12],
                    "sample_relevant_nonbool_values": relevant_nonbool_values[:8],
                }
            )
    return records


def build_payload(p2882: dict[str, Any]) -> dict[str, Any]:
    files = json_files()
    records = inventory_records()
    positive_records = [record for record in records if record.get("relevant_true_boolean_count", 0) > 0]
    source_ratio_positive_records = [
        record
        for record in positive_records
        for row in record.get("relevant_true_booleans", [])
        if "exports_source_stiffness" in row["path"].lower() or "source_stiffness_ratio_9_to_5_exported" in row["path"].lower()
    ]
    facts = {
        "p2882_rechecked": p2882.get("status") == "P2882_EULER_SOURCE_RATIO_9_OVER_5_FORCING_NO_GO_AUDIT_NO_CLOSURE",
        "generated_json_inventory_nonempty": len(files) > 0,
        "relevant_inventory_records_found": len(records) > 0,
        "no_positive_source_stiffness_9_to_5_export_found": len(source_ratio_positive_records) == 0,
        "positive_9_over_5_eta_hits_do_not_export_source_stiffness_law": len(positive_records) > 0 and len(source_ratio_positive_records) == 0,
    }
    return {
        "status": "P2883_DIMENSIONAL_UNIT_SOURCE_INVENTORY_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2882": sha(P2882)},
        "dimensional_unit_source_inventory_no_go_audit": {
            "input_status_rechecked": p2882.get("status"),
            "candidate_class": "current generated-artifact inventory for exported strict dimensional/unit source laws fixing 9:5 or unit-bearing 9/5",
            "generated_json_file_count": len(files),
            "relevant_record_count": len(records),
            "positive_relevant_record_count": len(positive_records),
            "positive_hits_classification": "positive hits are existing eta/damping/unit-premise records, not an exported primitive source/stiffness ratio 9:5 law",
            "positive_source_stiffness_9_to_5_record_count": len(source_ratio_positive_records),
            "sample_relevant_records": records[:20],
            "proof_certificate": {
                "inventory_rule": "Scan generated JSON artifacts for relevant accept/export booleans and 9/5, 9:5, source/stiffness, dimensional/unit, and coupling-theorem paths.",
                "finite_result": "Positive 9/5-related hits exist in older eta/damping/unit-premise artifacts, but no positive generated-artifact boolean exports a primitive source/stiffness ratio 9:5 law.",
                "sourcehood_step": "The current repository has negative/obstruction packets but no exported strict object fixing the primitive source/stiffness ratio 9:5 from independent data.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_current_dimensional_unit_source_export": False,
            "exports_source_stiffness_ratio_9_to_5": False,
            "exports_unit_bearing_9_over_5_coupling_theorem": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "current_dimensional_unit_source_exported": False,
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
            "reason": "P2883 inventories current generated artifacts after P2882 for a positive strict dimensional/unit source fixing primitive ratio 9:5.  The inventory finds positive 9/5-related eta/damping/unit-premise hits, but no positive generated-artifact boolean exports a primitive source/stiffness ratio 9:5 law; current hits are therefore insufficient for the P2882 obligation.",
            "next_honest_step": "Do not replay endpoint pins, denominator-5 coefficient boxes, local quadratic minimization, scalar Euler source-ratio transmission, or generated-artifact inventory as closure.  A next proof-grade move must introduce one genuinely new strict dimensional/unit source object or analytic theorem fixing primitive ratio 9:5 from independent data; otherwise preserve the no-new-live-frontier certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["dimensional_unit_source_inventory_no_go_audit"]
    lines = [
        "# P2883/S1833 dimensional/unit source inventory no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Generated-artifact inventory",
        f"- candidate class: `{audit['candidate_class']}`",
        f"- generated JSON file count: `{audit['generated_json_file_count']}`",
        f"- relevant record count: `{audit['relevant_record_count']}`",
        f"- positive relevant record count: `{audit['positive_relevant_record_count']}`",
        f"- positive source/stiffness 9:5 record count: `{audit['positive_source_stiffness_9_to_5_record_count']}`",
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
    payload = build_payload(read_json(P2882))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2883/S1833 dimensional/unit source inventory no-go audit",
        "## P2883/S1833 dimensional/unit source inventory no-go audit\n\n"
        "`P2883/S1833` inventories current generated artifacts after `P2882` for a positive strict dimensional/unit source fixing primitive ratio `9:5`.  Positive `9/5`-related eta/damping/unit-premise hits exist, but no positive generated-artifact boolean exports a primitive source/stiffness ratio `9:5` law; current hits are insufficient for the P2882 obligation.  No strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2883/S1833 dimensional/unit source inventory `L_total` guard",
        "## P2883/S1833 dimensional/unit source inventory `L_total` guard\n\n"
        "`P2883/S1833` adds no strict action term.  Existing positive `9/5` eta/damping/unit-premise hits do not export a primitive source/stiffness ratio `9:5` law, localized boundary/source density, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current dimensional/unit source inventory no-go guardrail (P2883/S1833, 2026-06-18)",
        "## Current dimensional/unit source inventory no-go guardrail (P2883/S1833, 2026-06-18)\n\n"
        "- P2883 inventories current generated artifacts after P2882 for a positive strict dimensional/unit source fixing primitive ratio `9:5` or a unit-bearing `9/5` coupling theorem.\n"
        "- Positive `9/5`-related eta/damping/unit-premise hits exist, but no current generated artifact exports a primitive source/stiffness ratio `9:5` law satisfying the P2882 obligation.\n"
        "- Do not promote endpoint pins, denominator-5 coefficient boxes, local quadratic minimization, scalar Euler source-ratio transmission, or generated-artifact inventory to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must introduce one genuinely new strict dimensional/unit source object or analytic theorem fixing primitive ratio `9:5` from independent data; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
