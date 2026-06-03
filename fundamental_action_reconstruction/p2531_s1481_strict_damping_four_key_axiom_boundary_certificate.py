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
OUT = GEN / "p2531_s1481_strict_damping_four_key_axiom_boundary_certificate.json"
MD = GEN / "p2531_s1481_strict_damping_four_key_axiom_boundary_certificate.md"

SOURCE_FILES = {
    "P2530_FOUR_KEY_IRREDUNDANCY": GEN / "p2530_s1480_strict_damping_four_key_irredundancy_witness_certificate.json",
}

SOURCE_KEYS = [
    "M_multiplicative_character_law_source",
    "P_prime_log_proportionality_source",
    "A_slope_value_or_prime_anchor_source",
    "O_m2_operator_signature_source",
]
NUMERIC_KEYS = SOURCE_KEYS[:3]
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
        "new_packet": "P2531|S1481|four-key axiom boundary|strict damping four-key axiom|four-key ternary source status|axiom boundary four-key",
        "intended_research_nonduplication": "ternary.*source|source status lattice|axiom.*strict.*lattice|four-key.*axiom|strict theorem status|non-strict.*source",
        "precursor_packets": "P2530|S1480|four-key irredundancy witness|P2517|dual-key axiom boundary|P2529|numeric subkey rank lattice",
        "source_obligation_language": "axiom-augmented|non-strict|strict source theorem|proper subset obstruction|strict_damping_beta_eta_source",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def classify_status_assignment(status_by_key: dict[str, str]) -> dict[str, Any]:
    numeric_statuses = [status_by_key[key] for key in NUMERIC_KEYS]
    all_statuses = [status_by_key[key] for key in SOURCE_KEYS]
    beta_eta_missing = any(status == "absent" for status in numeric_statuses)
    beta_eta_strict = all(status == "strict" for status in numeric_statuses)
    beta_eta_axiom_augmented = not beta_eta_missing and not beta_eta_strict
    strict_missing = any(status == "absent" for status in all_statuses)
    strict_accept = all(status == "strict" for status in all_statuses)
    strict_axiom_augmented = not strict_missing and not strict_accept
    if strict_accept:
        verdict = "strict_accept"
    elif strict_missing:
        verdict = "blocked_missing_key"
    else:
        verdict = "non_strict_axiom_augmented_only"
    return {
        "status_by_key": status_by_key,
        "beta_eta_numeric_source_strict_accepts": beta_eta_strict,
        "beta_eta_numeric_source_axiom_augmented_only": beta_eta_axiom_augmented,
        "beta_eta_numeric_source_blocked_missing_key": beta_eta_missing,
        "strict_damping_source_strict_accepts": strict_accept,
        "strict_damping_source_axiom_augmented_only": strict_axiom_augmented,
        "strict_damping_source_blocked_missing_key": strict_missing,
        "strict_damping_source_verdict": verdict,
    }


def ternary_status_lattice() -> list[dict[str, Any]]:
    rows = []
    for statuses in product(STATUSES, repeat=len(SOURCE_KEYS)):
        rows.append(classify_status_assignment(dict(zip(SOURCE_KEYS, statuses))))
    return rows


def append_doc_sections() -> None:
    eq_section = """
`P2531/S1481` refines the P2530 four-key irredundancy contract by adding a ternary source-status boundary.  Each required key is classified as `absent`, `axiom`, or `strict`; the `3^4=81` status table has exactly one strict accepting row, where all four keys are strict theorems.  Rows with any absent key are blocked, and rows with all keys present but at least one axiom key are explicitly non-strict/axiom-augmented only.  Thus axiom-supplying the multiplicative, prime-log, slope-anchor, or `m=2` operator key cannot be silently promoted to strict damping source closure.
""".strip()
    lag_section = """
`P2531/S1481` marks the four-key strict damping source boundary: all four P2530 keys must be strict theorems for strict source acceptance.  Any axiom-supplied key remains non-strict, and any absent key blocks the source; therefore no nonlinear compression-flow source, role-bearing `L_total`, bridge completion, role-transfer theorem, QW-2191 discharge, or ToE closure is licensed.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "P2531/S1481", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "P2531/S1481", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2530 = theorem(sources["P2530_FOUR_KEY_IRREDUNDANCY"], "strict_damping_four_key_irredundancy_witness_certificate")
    lattice = ternary_status_lattice()
    strict_rows = [row for row in lattice if row["strict_damping_source_strict_accepts"]]
    non_strict_rows = [row for row in lattice if row["strict_damping_source_axiom_augmented_only"]]
    blocked_rows = [row for row in lattice if row["strict_damping_source_blocked_missing_key"]]
    beta_eta_strict_rows = [row for row in lattice if row["beta_eta_numeric_source_strict_accepts"]]
    theorem_export = {
        "frontier_atom_under_attack": "four_key_strict_vs_axiom_augmented_source_status_boundary",
        "p2530_four_key_irredundancy_inherited": bool(p2530.get("four_key_irredundancy_witness_exported", False)),
        "source_keys": SOURCE_KEYS,
        "status_values": STATUSES,
        "ternary_status_row_count": len(lattice),
        "strict_damping_strict_accepting_row_count": len(strict_rows),
        "strict_damping_axiom_augmented_only_row_count": len(non_strict_rows),
        "strict_damping_blocked_missing_key_row_count": len(blocked_rows),
        "beta_eta_numeric_strict_accepting_row_count": len(beta_eta_strict_rows),
        "exactly_one_strict_accepting_row_all_keys_strict": len(strict_rows) == 1 and all(value == "strict" for value in strict_rows[0]["status_by_key"].values()),
        "all_axiom_augmented_rows_non_strict": len(non_strict_rows) == 15 and all(not row["strict_damping_source_strict_accepts"] for row in non_strict_rows),
        "all_missing_key_rows_blocked": len(blocked_rows) == 65 and all(row["strict_damping_source_blocked_missing_key"] for row in blocked_rows),
        "four_key_axiom_boundary_exported": True,
        "axiom_augmented_strict_damping_source_exported": False,
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
        "strict_damping_four_key_axiom_boundary_certificate": {
            "ternary_status_lattice": lattice,
            "source_obligation_normal_form": "strict_damping_beta_eta_source is strict iff all four P2530 keys are strict theorems; axiom-supplied present keys are non-strict/axiom-augmented only",
        },
    }
    gatekeepers = {
        "p2530_inherited": theorem_export["p2530_four_key_irredundancy_inherited"],
        "one_strict_accept_all_strict": theorem_export["exactly_one_strict_accepting_row_all_keys_strict"],
        "axiom_rows_non_strict_and_missing_rows_blocked": theorem_export["all_axiom_augmented_rows_non_strict"] and theorem_export["all_missing_key_rows_blocked"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "axiom_augmented_strict_damping_source_exported",
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
        "packet_id": "P2531",
        "stage_id": "S1481",
        "status": "STRICT_DAMPING_FOUR_KEY_AXIOM_BOUNDARY_CERTIFICATE_NONSTRICT_AXIOM_AUGMENTED_ONLY_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_four_key_axiom_boundary_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_four_key_axiom_boundary_certificate"]["theorem_export"]
    lines = [
        "# P2531/S1481 strict damping four-key axiom boundary certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2530 four-key irredundancy inherited: `{t['p2530_four_key_irredundancy_inherited']}`.",
        f"- Ternary status rows: `{t['ternary_status_row_count']}`.",
        f"- Strict accepting rows: `{t['strict_damping_strict_accepting_row_count']}`.",
        f"- Axiom-augmented-only rows: `{t['strict_damping_axiom_augmented_only_row_count']}`.",
        f"- Blocked missing-key rows: `{t['strict_damping_blocked_missing_key_row_count']}`.",
        f"- Exactly one strict accepting row with all keys strict: `{t['exactly_one_strict_accepting_row_all_keys_strict']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports only a ternary axiom-boundary table for the conditional four-key source contract. It does not source any key, does not promote axiom-augmented rows to strict closure, and exports no bridge completion, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_four_key_axiom_boundary_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_four_key_axiom_boundary_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
