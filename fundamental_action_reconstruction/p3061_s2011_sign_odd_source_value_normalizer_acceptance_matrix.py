#!/usr/bin/env python3
"""P3061/S2011: sign-odd source-value normalizer acceptance matrix.

P3060 showed that sign-even magnitude/support normalizers cannot fix the
P3059 coefficient sign.  P3061 constructs the missing kind of object instead:
a sign-odd source-value normalizer that would multiply a polarity-odd source
value sigma into the coefficient module.  The matrix separates formal
acceptance (what would work if a strict non-premise signed value existed) from
current export (what the repo actually supplies now).
"""
from __future__ import annotations

import hashlib, json, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3060_s2010_coefficient_sign_normalization_impossibility_verifier import OUT as P3060

OUT = GEN / "p3061_s2011_sign_odd_source_value_normalizer_acceptance_matrix.json"
MD = GEN / "p3061_s2011_sign_odd_source_value_normalizer_acceptance_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "sign_odd_strict_sourced_normalizer": r"sign-odd.*normalizer|strict-sourced coefficient-sign|coefficient-sign normalizer|source-value normalizer",
    "nonpremise_signed_value": r"non-premise signed value|nonpremise signed value|strict.*signed value|signed value.*source law",
    "polarity_selection_without_closure": r"plus_polarity|minus_polarity|polarity.*selection|no selector.*closure",
}

SOURCE_VALUE_CASES = [
    {"case": "no_sigma_export", "sigma": None, "source_status": "absent", "strict_nonpremise": False, "current_exported": False},
    {"case": "zero_sigma", "sigma": 0, "source_status": "zero_value", "strict_nonpremise": True, "current_exported": False},
    {"case": "premise_sigma_plus", "sigma": 1, "source_status": "premise_or_convention", "strict_nonpremise": False, "current_exported": False},
    {"case": "premise_sigma_minus", "sigma": -1, "source_status": "premise_or_convention", "strict_nonpremise": False, "current_exported": False},
    {"case": "strict_nonpremise_sigma_plus", "sigma": 1, "source_status": "formal_missing_object", "strict_nonpremise": True, "current_exported": False},
    {"case": "strict_nonpremise_sigma_minus", "sigma": -1, "source_status": "formal_missing_object", "strict_nonpremise": True, "current_exported": False},
]


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def evaluate_case(row: dict[str, Any]) -> dict[str, Any]:
    sigma = row["sigma"]
    nonzero = sigma in (-1, 1)
    formal_acceptance_if_exported = bool(nonzero and row["strict_nonpremise"])
    selected_polarity = "plus_polarity" if sigma == 1 else "minus_polarity" if sigma == -1 else None
    current_acceptance = bool(formal_acceptance_if_exported and row["current_exported"])
    blocker = "none" if current_acceptance else (
        "missing exported strict non-premise signed source value" if formal_acceptance_if_exported else
        "absent/zero/premise sigma cannot normalize the global coefficient sign"
    )
    return {
        **row,
        "sign_odd_under_sigma_flip": nonzero,
        "selected_polarity": selected_polarity,
        "formal_acceptance_if_exported": formal_acceptance_if_exported,
        "current_acceptance": current_acceptance,
        "blocker": blocker,
    }


def build_payload() -> dict[str, Any]:
    p3060 = read_json(P3060)
    grep_rows = content_grep()
    matrix = [evaluate_case(row) for row in SOURCE_VALUE_CASES]
    obligations = [
        {"obligation": "content_first_grep_before_source_value_matrix", "satisfied": True, "detail": "searched by sign-odd normalizer and non-premise signed-value content"},
        {"obligation": "construct_sign_odd_source_value_normalizer", "satisfied": True, "detail": "sigma-normalizer separates formal plus/minus polarity when sigma is strict non-premise and nonzero"},
        {"obligation": "separate_formal_acceptance_from_current_export", "satisfied": True, "detail": "formal missing-object cases are not counted as current exports"},
        {"obligation": "export_nonpremise_signed_value", "satisfied": False, "detail": "current artifacts provide no exported strict non-premise sigma value"},
        {"obligation": "export_selector_or_ltotal", "satisfied": False, "detail": "no QW-2191 discharge, selector closure, L_total, bridge, role transfer, or ToE closure follows"},
    ]
    return {
        "status": "P3061_SIGN_ODD_SOURCE_VALUE_NORMALIZER_ACCEPTANCE_MATRIX_NO_CURRENT_EXPORT",
        "input_hashes": {"P3060": hashlib.sha256(P3060.read_bytes()).hexdigest() if P3060.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": grep_rows,
            "normalizer_object": {
                "object": "SignOddSourceValueCoefficientSignNormalizer",
                "source_symbol": "sigma_selector",
                "target_missing_object": "strict_sourced_nonpremise_signed_value_for_coefficient_sign_normalization",
                "rule": "choose the coefficient-sign polarity by sign(sigma_selector) only when sigma_selector is strict, non-premise, nonzero, and coupled to G_selector",
                "acceptance_boundary": "formal construction is allowed, but current export requires an actual strict source law/value artifact",
            },
            "acceptance_matrix": matrix,
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(row["hit_count"] for row in grep_rows),
            "source_value_cases": len(matrix),
            "formal_acceptance_if_exported_cases": sum(1 for row in matrix if row["formal_acceptance_if_exported"]),
            "current_accepted_cases": sum(1 for row in matrix if row["current_acceptance"]),
            "premise_or_convention_cases_rejected": sum(1 for row in matrix if row["source_status"] == "premise_or_convention" and not row["current_acceptance"]),
            "p3060_status_seen": p3060.get("status"),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
        },
        "proof_obligations": obligations,
        "decision": {
            "bounded_no_go": "P3061 constructs the missing sign-odd source-value normalizer as an acceptance matrix.  Two formal rows would select plus_polarity or minus_polarity if a strict non-premise nonzero sigma_selector value were exported, but zero rows are currently accepted because current artifacts do not export that signed source value.  Premise/convention sigma rows are rejected.  No G_selector, QW-2191 discharge, selector closure, L_total, bridge/role transfer, or ToE closure is exported.",
            "negative_export_flags": {k: False for k in ["strict_nonpremise_sigma_selector_exported", "global_coefficient_sign_normalization_exported", "actual_g_selector_formula_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "The next proof-grade move should not build more formal sigma rows.  It must either exhibit one concrete strict source law computing a nonzero non-premise sigma_selector value coupled to G_selector, or pivot to a different P3057 atom with content-first grep and no selector/readout/L_total/bridge/ToE export.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3061/S2011 sign-odd source-value normalizer acceptance matrix", "", f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- source-value cases: `{c['source_value_cases']}`", f"- formal acceptance-if-exported cases: `{c['formal_acceptance_if_exported_cases']}`", f"- current accepted cases: `{c['current_accepted_cases']}`", f"- premise/convention cases rejected: `{c['premise_or_convention_cases_rejected']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3061/S2011 sign-odd source-value normalizer acceptance matrix", "## P3061/S2011 sign-odd source-value normalizer acceptance matrix\n\n`P3061/S2011` constructs the missing sign-odd normalizer shape after `P3060`: a `SignOddSourceValueCoefficientSignNormalizer` using a putative `sigma_selector` value.  The finite acceptance matrix has `6` source-value cases.  `2` rows would formally select `plus_polarity` or `minus_polarity` if an exported strict non-premise nonzero `sigma_selector` existed, but `0` rows are currently accepted because no such source law/value artifact is exported; premise/convention sigma rows are rejected.  No actual `G_selector`, `QW-2191` discharge, selector closure, `L_total`, bridge/role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3061/S2011 sign-odd source-value normalizer `L_total` guard", "## P3061/S2011 sign-odd source-value normalizer `L_total` guard\n\n`P3061/S2011` adds no physical `L_total` term.  It is a formal acceptance matrix for a missing signed source value and exports no unit-bearing action/EOM carrier or nonproxy variational chain rule.\n")
    append_once(AGENTS, "Current sign-odd source-value normalizer guardrail (P3061/S2011, 2026-06-23)", "## Current sign-odd source-value normalizer guardrail (P3061/S2011, 2026-06-23)\n\n- P3061 constructs the formal sign-odd normalizer shape requested after P3060: a `sigma_selector` source-value rule coupled to coefficient-sign normalization for `G_selector`.\n- The finite acceptance matrix checks `6` source-value cases; `2` rows would work only if an exported strict non-premise nonzero `sigma_selector` existed, while `0` rows are currently accepted and premise/convention sigma rows are rejected.\n- Do not promote formal sigma rows, premise/convention signs, or source-value placeholders to `QW-2191` discharge, selector closure, observed physics, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move must exhibit one concrete strict source law computing a nonzero non-premise `sigma_selector` value coupled to `G_selector`, or pivot to a different P3057 atom while preserving the P3048-P3061 bounded no-export boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
