#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from itertools import product
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
OUT = GEN / "p2517_s1467_strict_damping_dual_key_axiom_boundary_certificate.json"
MD = GEN / "p2517_s1467_strict_damping_dual_key_axiom_boundary_certificate.md"

SOURCE_FILES = {
    "P2516_DUAL_KEY_ACCEPTANCE": GEN / "p2516_s1466_strict_damping_dual_key_source_acceptance_matrix.json",
}

KEYS = ["beta_eta_numeric_source", "m2_operator_signature_source"]
STATES = ["absent", "axiom", "strict"]


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
        "new_packet": "P2517|S1467|dual-key axiom boundary|strict damping axiom boundary|non-strict source augmentation|axiom-augmented damping source",
        "precursor_packets": "P2516|S1466|dual-key source acceptance|P2515|operator-order signature acceptance|P2414|strict damping parameter identifiability",
        "axiom_language": "axiom-augmented|non-strict|proper subset obstruction|source acceptance|strict source theorem",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
        "closure_blockers": "bridge theorem|role-transfer theorem|physical-value generator|role-bearing L_total|selector closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def ternary_acceptance_table() -> dict[str, Any]:
    rows = []
    for numeric_state, signature_state in product(STATES, repeat=2):
        states = {
            "beta_eta_numeric_source": numeric_state,
            "m2_operator_signature_source": signature_state,
        }
        absent = [key for key, value in states.items() if value == "absent"]
        axiom = [key for key, value in states.items() if value == "axiom"]
        strict = [key for key, value in states.items() if value == "strict"]
        both_present = not absent
        strict_source = both_present and not axiom and len(strict) == len(KEYS)
        axiom_augmented = both_present and bool(axiom)
        blocked = bool(absent)
        if strict_source:
            classification = "strict_source_candidate_acceptance"
        elif axiom_augmented:
            classification = "non_strict_axiom_augmented_acceptance_only"
        else:
            classification = "blocked_missing_key"
        rows.append({
            "states": states,
            "absent_keys": absent,
            "axiom_keys": axiom,
            "strict_keys": strict,
            "strict_source_candidate_acceptance": strict_source,
            "non_strict_axiom_augmented_acceptance_only": axiom_augmented,
            "blocked_missing_key": blocked,
            "classification": classification,
        })
    return {
        "key_state_domain": STATES,
        "rows": rows,
        "row_count": len(rows),
        "strict_accepting_row_count": sum(1 for row in rows if row["strict_source_candidate_acceptance"]),
        "axiom_augmented_row_count": sum(1 for row in rows if row["non_strict_axiom_augmented_acceptance_only"]),
        "blocked_row_count": sum(1 for row in rows if row["blocked_missing_key"]),
        "strict_acceptance_requires_both_keys_strict": all(
            row["states"] == {"beta_eta_numeric_source": "strict", "m2_operator_signature_source": "strict"}
            for row in rows
            if row["strict_source_candidate_acceptance"]
        ),
        "axiom_rows_never_strict": all(
            not row["strict_source_candidate_acceptance"]
            for row in rows
            if row["axiom_keys"]
        ),
        "missing_key_rows_always_blocked": all(
            row["blocked_missing_key"]
            for row in rows
            if row["absent_keys"]
        ),
    }


def minimal_upgrade_paths() -> dict[str, Any]:
    rows = []
    base = {key: "absent" for key in KEYS}
    rows.append({
        "start_state": base,
        "upgrade": "add both keys as strict theorems",
        "end_state": {key: "strict" for key in KEYS},
        "result": "strict_source_candidate_acceptance",
        "strict": True,
    })
    for key in KEYS:
        end = {other: "strict" for other in KEYS}
        end[key] = "axiom"
        rows.append({
            "start_state": base,
            "upgrade": f"add {key} as axiom and the other key as strict theorem",
            "end_state": end,
            "result": "non_strict_axiom_augmented_acceptance_only",
            "strict": False,
        })
    rows.append({
        "start_state": base,
        "upgrade": "add both keys as axioms",
        "end_state": {key: "axiom" for key in KEYS},
        "result": "non_strict_axiom_augmented_acceptance_only",
        "strict": False,
    })
    return {
        "rows": rows,
        "strict_upgrade_path_count": sum(1 for row in rows if row["strict"]),
        "non_strict_upgrade_path_count": sum(1 for row in rows if not row["strict"]),
        "only_all_strict_path_is_strict": sum(1 for row in rows if row["strict"]) == 1,
    }


def build_axiom_boundary_certificate(p2516: dict[str, Any]) -> dict[str, Any]:
    table = ternary_acceptance_table()
    upgrades = minimal_upgrade_paths()
    p2516_dual_key = p2516.get("dual_key_acceptance_normal_form_exported") is True
    return {
        "frontier_atom_under_attack": "strict_damping_beta_eta_source",
        "p2516_dual_key_normal_form_inherited": p2516_dual_key,
        "certificate_type": "dual-key axiom boundary certificate for strict vs non-strict damping source acceptance",
        "ternary_acceptance_table": table,
        "minimal_upgrade_paths": upgrades,
        "strict_acceptance_requires_both_keys_as_strict_theorems": table["strict_acceptance_requires_both_keys_strict"],
        "axiom_augmented_rows_are_non_strict": table["axiom_rows_never_strict"],
        "missing_key_rows_blocked": table["missing_key_rows_always_blocked"],
        "only_all_strict_upgrade_path_is_strict": upgrades["only_all_strict_path_is_strict"],
        "axiom_boundary_exported": True,
        "strict_damping_beta_eta_source_exported": False,
        "axiom_augmented_non_strict_source_closure_exported": False,
        "damping_compression_bridge_component_ready": False,
        "full_bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_claimed": False,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2517/S1467 strict damping dual-key axiom boundary certificate

`P2517/S1467` refines the P2516 dual-key source acceptance matrix by separating strict theorem status from axiom-augmented status.  Each required key, `beta_eta_numeric_source` and `m2_operator_signature_source`, is classified as `absent`, `axiom`, or `strict`.  The ternary table has exactly one strict accepting row: both keys strict.  Rows with a missing key remain blocked, while rows where both keys are present but at least one key is supplied by axiom are explicitly non-strict/axiom-augmented only.

This prevents a common false pass: adding the numeric target or the `m=2` operator signature as an axiom may define a non-strict working closure, but it is not a strict source theorem and does not complete the damping bridge, role-transfer audit, QW-2191 selector closure, role-bearing `L_total`, physical-value generator, or ToE closure.
"""
    lag_section = """
## P2517/S1467 dual-key axiom boundary guard

`P2517/S1467` marks the strict damping source boundary: both P2516 keys must be strict theorems for strict source acceptance.  Any axiom-supplied key is non-strict by construction, and any missing key blocks the source.  Therefore no nonlinear compression-flow source theorem or role-bearing `L_total` term is licensed by axiom-augmented bookkeeping.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2517/S1467 strict damping dual-key axiom boundary certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2517/S1467 dual-key axiom boundary guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2516 = theorem(sources["P2516_DUAL_KEY_ACCEPTANCE"], "strict_damping_dual_key_source_acceptance_matrix")
    cert = build_axiom_boundary_certificate(p2516)
    table = cert["ternary_acceptance_table"]
    upgrades = cert["minimal_upgrade_paths"]
    theorem_export = {
        "theorem_name": "P2517_T1_strict_damping_dual_key_axiom_boundary_certificate",
        "audited_chain": ["P2516/S1466"],
        "strict_damping_dual_key_axiom_boundary_certificate": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2516_dual_key_normal_form_inherited": cert["p2516_dual_key_normal_form_inherited"],
        "ternary_row_count": table["row_count"],
        "strict_accepting_row_count": table["strict_accepting_row_count"],
        "axiom_augmented_row_count": table["axiom_augmented_row_count"],
        "blocked_row_count": table["blocked_row_count"],
        "strict_acceptance_requires_both_keys_as_strict_theorems": cert["strict_acceptance_requires_both_keys_as_strict_theorems"],
        "axiom_augmented_rows_are_non_strict": cert["axiom_augmented_rows_are_non_strict"],
        "missing_key_rows_blocked": cert["missing_key_rows_blocked"],
        "only_all_strict_upgrade_path_is_strict": cert["only_all_strict_upgrade_path_is_strict"],
        "strict_upgrade_path_count": upgrades["strict_upgrade_path_count"],
        "non_strict_upgrade_path_count": upgrades["non_strict_upgrade_path_count"],
        "axiom_boundary_exported": cert["axiom_boundary_exported"],
        "strict_damping_beta_eta_source_exported": False,
        "axiom_augmented_non_strict_source_closure_exported": False,
        "damping_compression_bridge_component_ready": False,
        "full_bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_claimed": False,
        "not_licensed": [
            "P2517 exports an axiom-boundary classification, not a strict damping source theorem.",
            "Axiom-supplied keys are non-strict and cannot be silently promoted to strict source closure.",
            "Missing keys remain blocked by P2516 proper-subset obstruction.",
            "No damping bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Replace both P2516 keys by strict theorems from nadsoliton dynamics; otherwise explicitly label the route axiom-augmented/non-strict.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2516_dual_key_inherited": theorem_export["p2516_dual_key_normal_form_inherited"],
        "ternary_table_shape_ok": theorem_export["ternary_row_count"] == 9,
        "strict_row_unique": theorem_export["strict_accepting_row_count"] == 1 and theorem_export["strict_acceptance_requires_both_keys_as_strict_theorems"],
        "axiom_rows_non_strict": theorem_export["axiom_augmented_rows_are_non_strict"] and theorem_export["axiom_augmented_row_count"] == 3,
        "missing_rows_blocked": theorem_export["missing_key_rows_blocked"] and theorem_export["blocked_row_count"] == 5,
        "source_not_exported": not theorem_export["strict_damping_beta_eta_source_exported"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "strict_damping_beta_eta_source_exported",
            "axiom_augmented_non_strict_source_closure_exported",
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
        "packet_id": "P2517",
        "stage_id": "S1467",
        "status": "STRICT_DAMPING_DUAL_KEY_AXIOM_BOUNDARY_CERTIFICATE_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_dual_key_axiom_boundary_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_dual_key_axiom_boundary_certificate"]["theorem_export"]
    lines = [
        "# P2517/S1467 strict damping dual-key axiom boundary certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2516 dual-key normal form inherited: `{t['p2516_dual_key_normal_form_inherited']}`.",
        f"- Ternary row count: `{t['ternary_row_count']}`.",
        f"- Strict accepting row count: `{t['strict_accepting_row_count']}`.",
        f"- Axiom-augmented row count: `{t['axiom_augmented_row_count']}`.",
        f"- Blocked row count: `{t['blocked_row_count']}`.",
        f"- Strict acceptance requires both keys as strict theorems: `{t['strict_acceptance_requires_both_keys_as_strict_theorems']}`.",
        f"- Axiom rows are non-strict: `{t['axiom_augmented_rows_are_non_strict']}`.",
        f"- Strict source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports an axiom-boundary classification only. Axiom-supplied keys are non-strict and missing keys remain blocked; no strict source atom, bridge theorem, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total term, physical-value generator, or ToE closure is exported.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_dual_key_axiom_boundary_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_dual_key_axiom_boundary_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
