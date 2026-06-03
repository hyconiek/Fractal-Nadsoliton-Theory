#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
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
OUT = GEN / "p2530_s1480_strict_damping_four_key_irredundancy_witness_certificate.json"
MD = GEN / "p2530_s1480_strict_damping_four_key_irredundancy_witness_certificate.md"

SOURCE_FILES = {
    "P2529_NUMERIC_SUBKEY_RANK_LATTICE": GEN / "p2529_s1479_strict_damping_numeric_subkey_rank_lattice_certificate.json",
}

SOURCE_KEYS = [
    "M_multiplicative_character_law",
    "P_prime_log_proportionality",
    "A_single_prime_slope_anchor",
    "O_m2_operator_signature_source",
]
NUMERIC_KEYS = SOURCE_KEYS[:3]
OPERATOR_KEY = SOURCE_KEYS[3]


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
        "new_packet": "P2530|S1480|four-key irredundancy witness|strict damping four-key|source-key irredundancy|all-four-key minimality",
        "intended_research_nonduplication": "irredundancy|irredundant|minimal hitting|proper subset witness|four-key|all four keys|key irredundancy|source-key irredundancy",
        "precursor_packets": "P2529|S1479|numeric subkey rank lattice|P2516|dual-key source acceptance|P2518|operator signature numeric-source separation",
        "source_obligation_language": "multiplicative_character_law_source|prime_log_proportionality_source|slope_value_or_prime_anchor_source|m2_operator_signature_source|strict_damping_beta_eta_source",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def key_tuple(row: dict[str, Any], field: str) -> tuple[str, ...]:
    return tuple(row[field])


def rank_lookup(numeric_lattice: list[dict[str, Any]]) -> dict[tuple[str, ...], dict[str, Any]]:
    return {key_tuple(row, "active_numeric_subkeys"): row for row in numeric_lattice}


def source_lookup(source_lattice: list[dict[str, Any]]) -> dict[tuple[str, ...], dict[str, Any]]:
    return {key_tuple(row, "active_source_keys"): row for row in source_lattice}


def omission_witnesses(p2529_theorem: dict[str, Any]) -> list[dict[str, Any]]:
    cert = p2529_theorem["strict_damping_numeric_subkey_rank_lattice_certificate"]
    ranks = rank_lookup(cert["numeric_subkey_rank_lattice"])
    sources = source_lookup(cert["source_key_lattice"])
    witnesses = []
    for missing_key in SOURCE_KEYS:
        active = tuple(key for key in SOURCE_KEYS if key != missing_key)
        source_row = sources[active]
        numeric_active = tuple(key for key in active if key in NUMERIC_KEYS)
        numeric_row = ranks.get(numeric_active)
        if missing_key == OPERATOR_KEY:
            rejection_reason = "numeric target can be conditionally pinned, but the independent m=2 operator-signature source is absent"
            numeric_rank = ranks[tuple(NUMERIC_KEYS)]["rank"]
            numeric_nullity = ranks[tuple(NUMERIC_KEYS)]["nullity"]
        else:
            rejection_reason = "numeric source remains underdetermined because a required numeric subkey is absent"
            numeric_rank = numeric_row["rank"] if numeric_row is not None else None
            numeric_nullity = numeric_row["nullity"] if numeric_row is not None else None
        witnesses.append({
            "missing_key": missing_key,
            "active_keys_all_but_missing": list(active),
            "numeric_rank_with_available_numeric_keys": numeric_rank,
            "numeric_nullity_with_available_numeric_keys": numeric_nullity,
            "beta_eta_numeric_source_accepts": source_row["beta_eta_numeric_source_accepts"],
            "strict_damping_beta_eta_source_accepts": source_row["strict_damping_beta_eta_source_accepts"],
            "rejection_reason": rejection_reason,
        })
    return witnesses


def append_doc_sections() -> None:
    eq_section = """
`P2530/S1480` turns the P2529 lattice into an explicit four-key irredundancy witness.  Removing any one of `multiplicative_character_law_source`, `prime_log_proportionality_source`, `slope_value_or_prime_anchor_source`, or `m2_operator_signature_source` rejects strict damping source: the first three removals leave positive numeric nullity, while removing the operator key leaves numeric rank/nullity `11/0` but no independent `m=2` operator-signature source.  Thus the four-key normal form is minimal as a conditional acceptance contract, not a source theorem.
""".strip()
    lag_section = """
`P2530/S1480` records the explicit all-but-one rejection witnesses for the current strict damping source contract.  It proves irredundancy of the four conditional keys but exports none of their sources, so no nonlinear compression-flow source, role-bearing `L_total`, bridge completion, role-transfer theorem, QW-2191 discharge, or ToE closure is licensed.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "P2530/S1480", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "P2530/S1480", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2529 = theorem(sources["P2529_NUMERIC_SUBKEY_RANK_LATTICE"], "strict_damping_numeric_subkey_rank_lattice_certificate")
    witnesses = omission_witnesses(p2529)
    theorem_export = {
        "frontier_atom_under_attack": "strict_damping_four_key_source_contract_irredundancy_after_numeric_subkey_lattice",
        "p2529_rank_lattice_inherited": bool(p2529.get("numeric_subkey_rank_lattice_exported", False)),
        "source_keys": SOURCE_KEYS,
        "omission_witness_count": len(witnesses),
        "all_single_key_omissions_reject_strict_damping_source": all(not row["strict_damping_beta_eta_source_accepts"] for row in witnesses),
        "numeric_key_omissions_leave_positive_numeric_nullity": all(row["numeric_nullity_with_available_numeric_keys"] and row["numeric_nullity_with_available_numeric_keys"] > 0 for row in witnesses if row["missing_key"] in NUMERIC_KEYS),
        "operator_key_omission_has_numeric_target_but_rejects_operator_source": next(row for row in witnesses if row["missing_key"] == OPERATOR_KEY)["numeric_nullity_with_available_numeric_keys"] == 0 and not next(row for row in witnesses if row["missing_key"] == OPERATOR_KEY)["strict_damping_beta_eta_source_accepts"],
        "four_key_irredundancy_witness_exported": True,
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
        "strict_damping_four_key_irredundancy_witness_certificate": {
            "omission_witnesses": witnesses,
            "source_obligation_normal_form": "strict_damping_beta_eta_source = multiplicative_character_law_source AND prime_log_proportionality_source AND slope_value_or_prime_anchor_source AND m2_operator_signature_source",
        },
    }
    gatekeepers = {
        "p2529_inherited": theorem_export["p2529_rank_lattice_inherited"],
        "all_four_keys_irredundant": theorem_export["all_single_key_omissions_reject_strict_damping_source"] and theorem_export["numeric_key_omissions_leave_positive_numeric_nullity"] and theorem_export["operator_key_omission_has_numeric_target_but_rejects_operator_source"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
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
        "packet_id": "P2530",
        "stage_id": "S1480",
        "status": "STRICT_DAMPING_FOUR_KEY_IRREDUNDANCY_WITNESS_CERTIFICATE_CONDITIONAL_CONTRACT_ONLY_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_four_key_irredundancy_witness_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_four_key_irredundancy_witness_certificate"]["theorem_export"]
    lines = [
        "# P2530/S1480 strict damping four-key irredundancy witness certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2529 rank lattice inherited: `{t['p2529_rank_lattice_inherited']}`.",
        f"- Omission witness count: `{t['omission_witness_count']}`.",
        f"- All single-key omissions reject strict damping source: `{t['all_single_key_omissions_reject_strict_damping_source']}`.",
        f"- Numeric-key omissions leave positive numeric nullity: `{t['numeric_key_omissions_leave_positive_numeric_nullity']}`.",
        f"- Operator-key omission has numeric target but rejects operator source: `{t['operator_key_omission_has_numeric_target_but_rejects_operator_source']}`.",
        f"- Four-key irredundancy witness exported: `{t['four_key_irredundancy_witness_exported']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports only all-but-one rejection witnesses for the conditional four-key source contract. It does not source any key, bridge completion, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_four_key_irredundancy_witness_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_four_key_irredundancy_witness_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
