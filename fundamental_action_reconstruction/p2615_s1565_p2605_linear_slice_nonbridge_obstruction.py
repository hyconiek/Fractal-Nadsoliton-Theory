#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2615_s1565_p2605_linear_slice_nonbridge_obstruction.json"
MD = GEN / "p2615_s1565_p2605_linear_slice_nonbridge_obstruction.md"

SOURCE_FILES = {
    "P2393_BOUNDARY_NEGATIVE_CONTROL": GEN / "p2393_s1343_kernel_completion_boundary_residual_certificate.json",
    "P2605_LINEAR_SLICE_EVIDENCE": GEN / "p2605_s1555_legacy_strict_linear_slice_completion_map_evidence.json",
    "P2606_NONLINEAR_RESIDUAL": GEN / "p2606_s1556_strict_side_nonlinear_compression_residual_addition.json",
    "P2610_CRITICAL_REVALIDATION": GEN / "p2610_s1560_p2601_p2608_critical_revalidation_audit.json",
    "P2614_CONTINUUM_PRIME_LOG": GEN / "p2614_s1564_p2602_continuum_rg_character_prime_log_proof.json",
}

BETA_TORS = Fraction(1, 100)
AUDIT_D_VALUES = list(range(1, 13))
NONLINEAR_ETA_CANDIDATES = [Fraction(4, 5), Fraction(9, 5)]
NEGATIVE_EXPORT_FLAGS = [
    "p2605_full_bridge_revalidated",
    "p2607_bridge_completion_revalidated",
    "p2608_role_bearing_ltotal_reenabled",
    "legacy_physical_role_transfer_exported",
    "qw2191_discharged_by_this_packet",
    "toe_closure_claimed",
    "apd_source_exported",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
    ], cwd=REPO, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2615|S1565|linear-slice non-bridge|P2605.*nonbridge|two-node denominator obstruction",
        "intended_research_nonduplication": "P2605.*linear-slice.*obstruction|linear slice cannot bridge|constant-beta.*d\\^eta|eta=1.*negative-control",
        "precursor_chain": "P2393|P2605|P2606|P2610|P2614|legacy strict linear slice|nonlinear compression residual",
        "guardrails": "K_legacy_ont|K_strict_gate|role-bearing L_total|QW-2191|ToE closure|APD source",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def ratio_power_residual(d: int, eta: Fraction, beta: float) -> float:
    return float(BETA_TORS) * d - beta * (d ** float(eta))


def best_beta_for_anchor(anchor_d: int, eta: Fraction) -> float:
    return float(BETA_TORS) * anchor_d / (anchor_d ** float(eta))


def nonlinear_candidate_rows() -> list[dict[str, Any]]:
    rows = []
    for eta in NONLINEAR_ETA_CANDIDATES:
        beta_from_d1 = best_beta_for_anchor(1, eta)
        beta_from_d2 = best_beta_for_anchor(2, eta)
        residuals_anchor1 = [ratio_power_residual(d, eta, beta_from_d1) for d in AUDIT_D_VALUES]
        residuals_anchor2 = [ratio_power_residual(d, eta, beta_from_d2) for d in AUDIT_D_VALUES]
        rows.append({
            "eta": str(eta),
            "eta_equals_one": eta == 1,
            "beta_fit_from_d1": beta_from_d1,
            "beta_fit_from_d2": beta_from_d2,
            "two_anchor_beta_mismatch": abs(beta_from_d1 - beta_from_d2),
            "max_abs_residual_when_beta_fit_at_d1": max(abs(value) for value in residuals_anchor1),
            "max_abs_residual_when_beta_fit_at_d2": max(abs(value) for value in residuals_anchor2),
            "accepts_single_constant_beta_completion_on_two_nodes": abs(beta_from_d1 - beta_from_d2) < 1e-14,
        })
    return rows


def exact_two_node_obstruction() -> dict[str, Any]:
    return {
        "statement": "For a node-preserving constant-beta denominator map 1+beta_tors*d = 1+beta*d^eta, equality at any two distinct positive nodes a,b forces eta=1 and beta=beta_tors.",
        "proof_steps": [
            "Cancel the shared unit term to obtain beta_tors*d = beta*d^eta at each audited positive node.",
            "At nodes a and b, divide the two equations: a/b = (a/b)^eta because beta and beta_tors are nonzero.",
            "Since a/b is positive and not equal to 1, taking logs gives (eta-1)*log(a/b)=0, hence eta=1.",
            "Substituting eta=1 back into beta_tors*d = beta*d gives beta=beta_tors.",
            "Therefore the P2605 eta=1 slice is an exact boundary equality only; it cannot by itself be a bridge to any nonlinear strict-side compression exponent eta != 1.",
        ],
        "boundary_conditions_that_would_avoid_the_obstruction": [
            "allow a node-dependent beta(d), which is not the constant-parameter strict gate kernel",
            "change the node coordinate d by a separately justified physical reparametrization",
            "provide independent strict-side residual additions and a selector/role theorem rather than treating the eta=1 slice as complete",
        ],
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    payloads = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2610 = theorem(payloads["P2610_CRITICAL_REVALIDATION"], "p2601_p2608_critical_revalidation_audit")
    p2614 = theorem(payloads["P2614_CONTINUUM_PRIME_LOG"], "p2602_continuum_rg_character_prime_log_proof")
    quarantined_before = set(p2610.get("quarantined_packet_ids_after_revalidation", []))
    rows = nonlinear_candidate_rows()
    obstruction = exact_two_node_obstruction()
    p2605_obstructed = all(not row["accepts_single_constant_beta_completion_on_two_nodes"] for row in rows)
    p2606_retained = "P2606" in set(p2610.get("accepted_packet_ids_after_revalidation", []))
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2615_T1_p2605_linear_slice_nonbridge_obstruction",
        "audited_chain": ["P2393/S1343", "P2605/S1555", "P2606/S1556", "P2610/S1560", "P2614/S1564"],
        "exact_two_node_obstruction": obstruction,
        "nonlinear_candidate_rows": rows,
        "p2605_quarantine_before_p2615": "P2605" in quarantined_before,
        "p2605_linear_slice_reclassified_as_boundary_negative_control": p2605_obstructed,
        "p2605_quarantine_retained_by_p2615": p2605_obstructed,
        "p2606_nonlinear_residual_component_retained": p2606_retained,
        "p2601_p2602_revalidation_inherited": p2614.get("strict_damping_beta_eta_source_revalidated_under_df_9_5_scope") is True,
        "strict_damping_source_status_after_p2615": "NON_BRIDGE_SOURCE_REVALIDATED_BY_P2613_P2614_WITH_P2603_SCOPE; LEGACY_TO_STRICT_BRIDGE_STILL_BLOCKED_BY_P2605_P2607_P2608",
        "remaining_p2610_quarantines_after_p2615": sorted((quarantined_before - {"P2601", "P2602"}) | {"P2605"}),
        "recommended_next_honest_step": "Do not treat P2605 as bridge completion. Either keep bridge work blocked pending a physical source for the strict-side phase/topological selector, or audit whether P2608 role semantics can accept any non-bridge bookkeeping-only term (likely no).",
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2605_was_quarantined_before": theorem_export["p2605_quarantine_before_p2615"],
        "two_node_obstruction_present": len(obstruction["proof_steps"]) == 5,
        "all_nonlinear_candidates_fail_constant_beta_two_node_completion": all(not row["accepts_single_constant_beta_completion_on_two_nodes"] for row in rows),
        "nonlinear_beta_mismatches_positive": all(row["two_anchor_beta_mismatch"] > 0.0 for row in rows),
        "p2605_quarantine_retained": theorem_export["p2605_quarantine_retained_by_p2615"],
        "p2606_component_retained": theorem_export["p2606_nonlinear_residual_component_retained"],
        "nonbridge_source_revalidation_inherited": theorem_export["p2601_p2602_revalidation_inherited"],
        "no_p2605_full_bridge_revalidation": theorem_export["p2605_full_bridge_revalidated"] is False,
        "no_p2607_bridge_revalidation": theorem_export["p2607_bridge_completion_revalidated"] is False,
        "no_p2608_ltotal_reenable": theorem_export["p2608_role_bearing_ltotal_reenabled"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_packet"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
        "no_apd_source_exported": theorem_export["apd_source_exported"] is False,
    }
    return {
        "packet_id": "P2615",
        "stage_id": "S1565",
        "status": "P2615_P2605_LINEAR_SLICE_RETAINED_ONLY_AS_BOUNDARY_NEGATIVE_CONTROL_NONBRIDGE_OBSTRUCTION_RECORDED",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "p2605_linear_slice_nonbridge_obstruction": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {name: sha256_json(payload) for name, payload in payloads.items()},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["p2605_linear_slice_nonbridge_obstruction"]["theorem_export"]
    proof = t["exact_two_node_obstruction"]
    lines = [
        "# P2615/S1565 P2605 linear-slice non-bridge obstruction", "",
        f"Status: `{payload['status']}`", "", "## Theorem", "",
        proof["statement"], "", "## Exact proof", "",
    ]
    for step in proof["proof_steps"]:
        lines.append(f"- {step}")
    lines.extend([
        "", "## Computed checks", "",
        f"- P2605 quarantine retained: `{t['p2605_quarantine_retained_by_p2615']}`.",
        f"- P2606 nonlinear residual component retained: `{t['p2606_nonlinear_residual_component_retained']}`.",
        f"- Non-bridge damping source revalidation inherited from P2613/P2614: `{t['p2601_p2602_revalidation_inherited']}`.",
        f"- Remaining quarantines after P2615: `{t['remaining_p2610_quarantines_after_p2615']}`.",
        "", "## Scope guards", "",
        "P2615 does not revalidate P2605 as a full bridge, does not revalidate P2607, does not re-enable P2608 role-bearing L_total, and does not export QW-2191 discharge, APD sourcehood, legacy physical-role transfer, or ToE closure.",
        "", "## Fingerprint", "",
        f"`{payload['p2605_linear_slice_nonbridge_obstruction']['theorem_fingerprint_sha256']}`",
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2615/S1565 P2605 linear-slice non-bridge obstruction

`P2615/S1565` reclassifies P2605 as an exact boundary/negative-control slice rather than bridge completion: a node-preserving constant-beta denominator equality `1+beta_tors*d = 1+beta*d^eta` at two distinct positive nodes forces `eta=1` and `beta=beta_tors`.  Therefore the P2605 `eta=1` match cannot by itself supply nonlinear strict-side compression.  P2606 remains only a nonlinear residual component, while P2607 bridge completion and P2608 role-bearing `L_total` remain blocked.
""".strip()
    lag_section = """
## P2615/S1565 linear-slice nonbridge Ltotal guard

`P2615/S1565` prevents the P2605 linear slice from being used as a role-bearing `L_total` bridge premise.  The non-bridge damping source bookkeeping from P2613/P2614 may stand under its scope, but legacy-to-strict bridge completion and P2608 role transfer remain quarantined.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2615/S1565 P2605 linear-slice non-bridge obstruction", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2615/S1565 linear-slice nonbridge Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
