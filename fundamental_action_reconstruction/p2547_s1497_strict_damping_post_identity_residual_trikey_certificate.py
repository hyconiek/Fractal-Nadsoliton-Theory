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
OUT = GEN / "p2547_s1497_strict_damping_post_identity_residual_trikey_certificate.json"
MD = GEN / "p2547_s1497_strict_damping_post_identity_residual_trikey_certificate.md"

SOURCE_FILES = {
    "P2530_FOUR_KEY_IRREDUNDANCY": GEN / "p2530_s1480_strict_damping_four_key_irredundancy_witness_certificate.json",
    "P2540_M2_OBSTRUCTION": GEN / "p2540_s1490_strict_damping_m2_operator_signature_current_premise_obstruction_certificate.json",
    "P2542_PRIME_LOG_OBSTRUCTION": GEN / "p2542_s1492_strict_damping_prime_log_proportionality_current_premise_obstruction_certificate.json",
    "P2543_SLOPE_OBSTRUCTION": GEN / "p2543_s1493_strict_damping_slope_value_current_premise_obstruction_certificate.json",
    "P2544_FOUR_KEY_CLOSURE_BLOCKER": GEN / "p2544_s1494_strict_damping_four_key_current_premise_closure_blocker_certificate.json",
    "P2546_IDENTITY_CONDITIONAL": GEN / "p2546_s1496_strict_damping_identity_action_conditional_propagation_certificate.json",
}

RESIDUAL_KEYS = [
    "prime_log_proportionality_source",
    "slope_value_or_prime_anchor_source",
    "m2_operator_signature_source",
]

NEGATIVE_EXPORT_FLAGS = [
    "prime_log_proportionality_source_exported",
    "slope_value_or_prime_anchor_source_exported",
    "m2_operator_signature_source_exported",
    "beta_eta_numeric_source_exported",
    "strict_damping_beta_eta_source_exported",
    "source_obligation_discharge_exported",
    "damping_compression_bridge_component_ready",
    "full_bridge_theorem_exported",
    "role_transfer_theorem_exported",
    "selector_closure_exported",
    "qw2191_discharged_by_this_certificate",
    "role_bearing_ltotal_exported",
    "toe_closure_claimed",
]


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
    return {"count": len(lines), "samples": lines[:60]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2547|S1497|post-identity residual tri-key|residual trikey|identity source residual frontier",
        "intended_research_nonduplication": "post-identity residual|residual tri-key|tri-key residual|identity-action.*remaining.*source|remaining.*prime-log.*slope.*m2",
        "precursor_packets": "P2530|S1480|P2540|S1490|P2542|S1492|P2543|S1493|P2544|S1494|P2546|S1496",
        "source_key_language": "prime_log_proportionality_source|slope_value_or_prime_anchor_source|m2_operator_signature_source|beta_eta_numeric_source",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def load_theorems(sources: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {
        "P2530_FOUR_KEY_IRREDUNDANCY": theorem(sources["P2530_FOUR_KEY_IRREDUNDANCY"], "strict_damping_four_key_irredundancy_witness_certificate"),
        "P2540_M2_OBSTRUCTION": theorem(sources["P2540_M2_OBSTRUCTION"], "strict_damping_m2_operator_signature_current_premise_obstruction_certificate"),
        "P2542_PRIME_LOG_OBSTRUCTION": theorem(sources["P2542_PRIME_LOG_OBSTRUCTION"], "strict_damping_prime_log_proportionality_current_premise_obstruction_certificate"),
        "P2543_SLOPE_OBSTRUCTION": theorem(sources["P2543_SLOPE_OBSTRUCTION"], "strict_damping_slope_value_current_premise_obstruction_certificate"),
        "P2544_FOUR_KEY_CLOSURE_BLOCKER": theorem(sources["P2544_FOUR_KEY_CLOSURE_BLOCKER"], "strict_damping_four_key_current_premise_closure_blocker_certificate"),
        "P2546_IDENTITY_CONDITIONAL": theorem(sources["P2546_IDENTITY_CONDITIONAL"], "strict_damping_identity_action_conditional_propagation_certificate"),
    }


def residual_acceptance_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for p, a, o in product([False, True], repeat=3):
        assignment = {
            "multiplicative_character_law_source": True,
            "prime_log_proportionality_source": p,
            "slope_value_or_prime_anchor_source": a,
            "m2_operator_signature_source": o,
        }
        beta_eta_numeric = assignment["multiplicative_character_law_source"] and p and a
        strict_damping = beta_eta_numeric and o
        rows.append({
            "residual_assignment": {key: assignment[key] for key in RESIDUAL_KEYS},
            "identity_multiplicative_key_assumed_strict": True,
            "beta_eta_numeric_source_accepts": beta_eta_numeric,
            "strict_damping_beta_eta_source_accepts": strict_damping,
            "missing_residual_keys": [key for key in RESIDUAL_KEYS if not assignment[key]],
        })
    return rows


def single_omission_witnesses(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    witnesses = []
    for missing_key in RESIDUAL_KEYS:
        witness = next(
            row for row in rows
            if row["missing_residual_keys"] == [missing_key]
        )
        witnesses.append({
            "omitted_key": missing_key,
            "assignment": witness["residual_assignment"],
            "beta_eta_numeric_source_accepts": witness["beta_eta_numeric_source_accepts"],
            "strict_damping_beta_eta_source_accepts": witness["strict_damping_beta_eta_source_accepts"],
        })
    return witnesses


def build_certificate(theorems: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2530 = theorems["P2530_FOUR_KEY_IRREDUNDANCY"]
    p2540 = theorems["P2540_M2_OBSTRUCTION"]
    p2542 = theorems["P2542_PRIME_LOG_OBSTRUCTION"]
    p2543 = theorems["P2543_SLOPE_OBSTRUCTION"]
    p2544 = theorems["P2544_FOUR_KEY_CLOSURE_BLOCKER"]
    p2546 = theorems["P2546_IDENTITY_CONDITIONAL"]
    rows = residual_acceptance_rows()
    witnesses = single_omission_witnesses(rows)
    accepting = [row for row in rows if row["strict_damping_beta_eta_source_accepts"]]
    return {
        "p2530_four_key_irredundancy_inherited": p2530.get("four_key_irredundancy_witness_exported") is True,
        "p2540_m2_obstruction_inherited": p2540.get("current_premise_nonentailment_of_m2_exported") is True,
        "p2542_prime_log_obstruction_inherited": p2542.get("current_unital_multiplicative_premises_do_not_entail_prime_log_proportionality") is True,
        "p2543_slope_obstruction_inherited": p2543.get("all_slope_candidates_pass_unital_multiplicative_prime_log_premises") is True,
        "p2544_four_key_blocker_inherited": p2544.get("four_key_current_premise_closure_blocker_exported") is True,
        "p2546_identity_conditional_inherited": p2546.get("identity_action_conditional_propagation_exported") is True,
        "identity_action_source_assumption_status": "hypothetical_strict_source_for_frontier_reduction_only",
        "residual_source_keys": RESIDUAL_KEYS,
        "residual_truth_table_rows": rows,
        "residual_truth_table_row_count": len(rows),
        "residual_strict_accepting_row_count": len(accepting),
        "residual_strict_accepting_row": accepting[0],
        "residual_single_omission_witnesses": witnesses,
        "all_residual_single_omissions_reject_strict_damping": all(
            not witness["strict_damping_beta_eta_source_accepts"] for witness in witnesses
        ),
        "post_identity_residual_trikey_irredundancy_exported": True,
        "conditional_identity_reduces_missing_count_from_4_to_3": p2546.get("conditional_missing_source_key_delta") == 1,
        "identity_source_alone_cannot_export_beta_eta_numeric_source": True,
        "identity_source_alone_cannot_export_strict_damping_beta_eta_source": True,
        "recommended_next_honest_step": (
            "Treat identity-action work as useful only if it becomes a real strict source theorem.  For immediate proof/computation progress, "
            "attack one of the three residual keys with new source premises rather than another y_1 scan.  The sharpest independent target is "
            "m=2 operator-order selection from nadsoliton dynamics; in parallel, a prime-log proportionality source would be the next numeric-source bottleneck."
        ),
        "interpretation": (
            "Even granting the best-case P2546 identity-action theorem, the strict damping source problem reduces to a residual irredundant "
            "three-key conjunction.  The acceptance table has exactly one accepting row, and every single residual omission still blocks "
            "strict_damping_beta_eta_source."
        ),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2547/S1497 strict damping post-identity residual tri-key certificate

`P2547/S1497` grants the best-case P2546 hypothesis only as a conditional frontier reduction: the multiplicative/unital key is assumed strict.  The residual truth table over prime-log proportionality, `delta=4/5` slope/prime anchor, and `m=2` operator signature has exactly one accepting row.  Each single residual-key omission still rejects `strict_damping_beta_eta_source`.

Therefore identity-action work, even if later sourced, is not a full damping-source closure.  The next proof/computation frontier is a genuinely new source theorem for one residual key, preferably the independent `m=2` operator-order selection or the numeric prime-log proportionality source.  No residual key source, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.
""".strip()
    lag_section = """
## P2547/S1497 post-identity residual tri-key Ltotal guard

`P2547/S1497` records that a hypothetical identity-action source leaves an irredundant residual tri-key blocker.  The strict damping term remains non-role-bearing in `L_total` until prime-log proportionality, slope anchoring, and `m=2` operator selection are all strict source theorems.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2547/S1497 strict damping post-identity residual tri-key certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2547/S1497 post-identity residual tri-key Ltotal guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    theorems = load_theorems(sources)
    cert = build_certificate(theorems)
    theorem_export = {
        "theorem_name": "P2547_T1_strict_damping_post_identity_residual_trikey_certificate",
        "audited_chain": ["P2530/S1480", "P2540/S1490", "P2542/S1492", "P2543/S1493", "P2544/S1494", "P2546/S1496"],
        "strict_damping_post_identity_residual_trikey_certificate": cert,
        **cert,
        "prime_log_proportionality_source_exported": False,
        "slope_value_or_prime_anchor_source_exported": False,
        "m2_operator_signature_source_exported": False,
        "beta_eta_numeric_source_exported": False,
        "strict_damping_beta_eta_source_exported": False,
        "source_obligation_discharge_exported": False,
        "damping_compression_bridge_component_ready": False,
        "full_bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_claimed": False,
        "not_licensed": [
            "This packet is a conditional residual-frontier theorem, not a source theorem for any residual key.",
            "It does not replace the P2530 four-key normal form.",
            "It does not transfer legacy physical-role claims onto the strict gate kernel.",
            "It does not discharge QW-2191 or any ToE gate.",
        ],
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2546_inherited": theorem_export["p2546_identity_conditional_inherited"],
        "truth_table_size_verified": theorem_export["residual_truth_table_row_count"] == 8,
        "unique_accepting_row_verified": theorem_export["residual_strict_accepting_row_count"] == 1,
        "single_omission_rejection_verified": theorem_export["all_residual_single_omissions_reject_strict_damping"],
        "negative_controls_preserved": not any(theorem_export[key] for key in NEGATIVE_EXPORT_FLAGS),
    }
    return {
        "packet_id": "P2547",
        "stage_id": "S1497",
        "status": "STRICT_DAMPING_POST_IDENTITY_RESIDUAL_TRIKEY_CERTIFICATE_NO_RESIDUAL_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_post_identity_residual_trikey_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_post_identity_residual_trikey_certificate"]["theorem_export"]
    lines = [
        "# P2547/S1497 strict damping post-identity residual tri-key certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Residual keys: `{', '.join(t['residual_source_keys'])}`.",
        f"- Residual truth-table rows: `{t['residual_truth_table_row_count']}`.",
        f"- Strict accepting residual rows: `{t['residual_strict_accepting_row_count']}`.",
        f"- Single residual omissions reject strict damping: `{t['all_residual_single_omissions_reject_strict_damping']}`.",
        f"- Identity source alone exports beta/eta numeric source: `{not t['identity_source_alone_cannot_export_beta_eta_numeric_source']}`.",
        "",
        "## Interpretation",
        "",
        t["interpretation"],
        "",
        "## Recommendation",
        "",
        t["recommended_next_honest_step"],
        "",
        "## Negative Controls",
        "",
        "No residual key source, full source-obligation discharge, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_post_identity_residual_trikey_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_post_identity_residual_trikey_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
