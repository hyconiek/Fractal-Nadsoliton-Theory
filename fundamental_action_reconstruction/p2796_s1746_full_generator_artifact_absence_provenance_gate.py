#!/usr/bin/env python3
"""P2796/S1746: full-generator artifact absence / provenance gate.

P2795 showed that the recent named-subclass lane is saturated against P2791.
The next honest computational move is therefore not another named-subclass
certificate.  P2796 audits the repository for the object that the guardrails now
require: an actual certified full connected 16-node 4-regular generator artifact
or toolchain with graph6/hash provenance.

This script scans generated JSON certificates for full-generator acceptance
flags, explicit canonical-generator booleans, graph6/provenance traces, and known
local/lower-bound class counts.  It also probes common external generator
commands in the current environment.  The intended honest outcome on current
artifacts is a blocking provenance gate: no full generator artifact is present,
so no canonical geometry or K/L_total promotion is licensed.
"""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
P2795 = GEN / "p2795_s1745_post_subclass_saturation_no_new_class_frontier_certificate.json"
OUT = GEN / "p2796_s1746_full_generator_artifact_absence_provenance_gate.json"
MD = GEN / "p2796_s1746_full_generator_artifact_absence_provenance_gate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

GENERATOR_COMMAND_CANDIDATES = ["genreg", "geng", "shortg", "labelg", "dreadnaut", "plantri"]
FULL_GENERATOR_KEYS = [
    "accepted_as_full_16node_canonical_generator_certificate",
    "canonical_16node_generator_certified",
    "canonical_generator_certified",
    "full_16node_canonical_generator_certified",
]
NEGATIVE_EXPORT_FLAGS = [
    "canonical_16node_generator_certified",
    "canonical_geometry_source_exported",
    "strict_spectral_source_law_exported",
    "global_full_spectrum_geometry_theorem_exported",
    "kernel_geometry_closure_exported",
    "kernel_fully_expresses_nadsoliton_characterISTICS".lower(),
    "role_bearing_ltotal_promoted",
    "bridge_closure_exported",
    "selector_closure_exported",
    "toe_closure_exported",
]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def walk_values(obj: Any, path: str = "") -> list[tuple[str, Any]]:
    rows: list[tuple[str, Any]] = []
    if isinstance(obj, dict):
        for key, value in obj.items():
            child = f"{path}.{key}" if path else key
            rows.append((child, value))
            rows.extend(walk_values(value, child))
    elif isinstance(obj, list):
        for index, value in enumerate(obj):
            child = f"{path}[{index}]"
            rows.extend(walk_values(value, child))
    return rows


def generated_json_files() -> list[Path]:
    return sorted(path for path in GEN.glob("*.json") if path.is_file())


def status_mentions_full_generator(payload: dict[str, Any]) -> bool:
    text = json.dumps(payload, sort_keys=True, ensure_ascii=False).lower()
    return "full connected 16-node 4-regular generator" in text or "full 16-node" in text


def artifact_audit() -> dict[str, Any]:
    json_paths = generated_json_files()
    full_generator_truth_rows = []
    false_generator_guard_rows = []
    graph6_rows = []
    status_rows = []
    for path in json_paths:
        try:
            payload = read_json(path)
        except json.JSONDecodeError:
            continue
        if isinstance(payload.get("status"), str):
            status_rows.append({"file": rel(path), "status": payload["status"]})
        for key_path, value in walk_values(payload):
            key_tail = key_path.split(".")[-1]
            if key_tail in FULL_GENERATOR_KEYS:
                row = {"file": rel(path), "key_path": key_path, "value": value}
                if value is True:
                    full_generator_truth_rows.append(row)
                else:
                    false_generator_guard_rows.append(row)
            if key_tail == "graph6" and isinstance(value, str):
                graph6_rows.append({"file": rel(path), "key_path": key_path, "graph6_sha256": hashlib.sha256(value.encode("ascii", errors="ignore")).hexdigest()})
    command_rows = [{"command": command, "path": shutil.which(command), "available": shutil.which(command) is not None} for command in GENERATOR_COMMAND_CANDIDATES]
    mention_files = [rel(path) for path in json_paths if status_mentions_full_generator(read_json(path))]
    return {
        "generated_json_file_count": len(json_paths),
        "status_row_count": len(status_rows),
        "status_rows_with_no_closure": sum(1 for row in status_rows if "NO_CLOSURE" in row["status"]),
        "graph6_payload_row_count": len(graph6_rows),
        "full_generator_truth_rows": full_generator_truth_rows,
        "full_generator_truth_row_count": len(full_generator_truth_rows),
        "false_generator_guard_row_count": len(false_generator_guard_rows),
        "full_generator_mention_files": mention_files,
        "full_generator_mention_file_count": len(mention_files),
        "external_generator_command_rows": command_rows,
        "external_generator_command_available_count": sum(1 for row in command_rows if row["available"]),
        "required_artifact_present": len(full_generator_truth_rows) > 0,
        "finite_certificate_statement": "No generated JSON artifact currently asserts a true full 16-node canonical generator certificate, while current graph6/provenance rows remain local/subclass/lower-bound evidence only.",
    }


def acceptance_matrix(witness: dict[str, Any], p2795: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2795_saturation_gate_present": p2795.get("status") == "P2795_POST_SUBCLASS_SATURATION_NO_NEW_CLASS_FRONTIER_CERTIFICATE_NO_CLOSURE",
        "generated_json_scanned": witness["generated_json_file_count"] > 0,
        "graph6_payloads_exist_but_are_not_full_generator": witness["graph6_payload_row_count"] > 0 and witness["full_generator_truth_row_count"] == 0,
        "no_true_full_generator_certificate_found": witness["full_generator_truth_row_count"] == 0,
        "external_generator_toolchain_not_available_in_environment": witness["external_generator_command_available_count"] == 0,
        "full_generator_artifact_present": False,
        "strict_nadsoliton_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_full_generator_absence_provenance_gate": all(facts[key] for key in [
            "p2795_saturation_gate_present",
            "generated_json_scanned",
            "graph6_payloads_exist_but_are_not_full_generator",
            "no_true_full_generator_certificate_found",
        ]),
        "accepted_as_full_16node_canonical_generator_certificate": False,
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "The repository contains local graph6/provenance and no-new-class certificates, but no artifact or available toolchain certifies the full connected 16-node 4-regular class.  Promotion remains blocked until such an artifact/toolchain or a strict spectral source law is supplied.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    w = payload["full_generator_artifact_absence_witness"]
    lines = [
        "# P2796/S1746 full-generator artifact absence provenance gate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Repository/toolchain audit result",
        f"- generated_json_file_count={w['generated_json_file_count']}",
        f"- graph6_payload_row_count={w['graph6_payload_row_count']}",
        f"- full_generator_truth_row_count={w['full_generator_truth_row_count']}",
        f"- false_generator_guard_row_count={w['false_generator_guard_row_count']}",
        f"- full_generator_mention_file_count={w['full_generator_mention_file_count']}",
        f"- external_generator_command_available_count={w['external_generator_command_available_count']}",
        f"- required_artifact_present={w['required_artifact_present']}",
        "",
        "## Decision",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2795 = read_json(P2795)
    witness = artifact_audit()
    acceptance = acceptance_matrix(witness, p2795)
    payload = {
        "status": "P2796_FULL_GENERATOR_ARTIFACT_ABSENCE_PROVENANCE_GATE_NO_CLOSURE",
        "input_hashes": {"P2795": sha(P2795)},
        "input_statuses": {"P2795": p2795.get("status")},
        "audited_question": "Is there currently a repository artifact or local toolchain certifying the full connected 16-node 4-regular generator class required after P2795?",
        "full_generator_artifact_absence_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not manufacture another named-subclass replay.  The next honest move is to supply/import a real certified full connected 16-node 4-regular generator artifact/toolchain with graph6/hash provenance and then run the full quotient/charpoly/complement/orbit audit; alternatively export a strict nadsoliton spectral action/source law fixing the admissible class, target spectrum, and K/L_total coupling before testing.  Otherwise preserve the P2697-P2796 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2796/S1746 full-generator artifact absence provenance gate", "## P2796/S1746 full-generator artifact absence provenance gate\n\n`P2796/S1746` audits the repository after P2795 for the required full connected 16-node 4-regular generator artifact/toolchain.  It scans generated JSON certificates for true full-generator acceptance flags, graph6/provenance rows, and available local generator commands.  The current result is a blocking provenance gate: graph6/provenance evidence exists, but no generated artifact asserts a true full-generator certificate and no local external generator command is available in the audited environment.  This is not a generator, not a strict spectral source law, and not a `K`/`L_total` variational coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2796/S1746 full-generator absence Ltotal guard", "## P2796/S1746 full-generator absence Ltotal guard\n\n`P2796/S1746` adds no variational source term.  It records that the required full connected 16-node 4-regular generator artifact/toolchain is absent on current repository artifacts, so no `K`/`L_total` source or canonical geometry theorem can be promoted from the local graph bookkeeping.\n")
    append_once(AGENTS, "Current full-generator artifact absence guardrail (P2796/S1746, 2026-06-16)", "## Current full-generator artifact absence guardrail (P2796/S1746, 2026-06-16)\n\n- P2796 audits current generated JSON/provenance artifacts and local generator commands after P2795; it finds graph6/provenance rows but zero true full connected 16-node 4-regular generator certificates and no available audited external generator command.\n- Do not continue named-subclass replay or promote local graph bookkeeping to canonical geometry; the missing object is still an actual certified full generator artifact/toolchain or a strict spectral source law.\n- Do not promote the absence/provenance gate to canonical geometry, strict spectral source law, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.\n")
    return payload


if __name__ == "__main__":
    main()
