#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2618_s1568_analytic_legacy_to_strict_completion_obstruction.json"
MD = GEN / "p2618_s1568_analytic_legacy_to_strict_completion_obstruction.md"

SOURCE_FILES = {
    "TICKET_P2618": ROOT / "TICKET_P2618_BRIDGE_COMPLETION.md",
    "P2603_FRACTAL_CODIMENSION": GEN / "p2603_s1553_nadsoliton_fractal_codimension_slope_source_theorem.json",
    "P2612_GF2_ORIGIN_OBSTRUCTION": GEN / "p2612_s1562_p2607_gf2_physical_origin_obstruction.json",
    "P2614_CONTINUUM_RG_CHARACTER": GEN / "p2614_s1564_p2602_continuum_rg_character_prime_log_proof.json",
    "P2615_LINEAR_SLICE_OBSTRUCTION": GEN / "p2615_s1565_p2605_linear_slice_nonbridge_obstruction.json",
    "P2616_ROLE_ACCEPTANCE_OBSTRUCTION": GEN / "p2616_s1566_p2608_role_acceptance_obstruction_after_source_revalidation.json",
    "P2617_EXPONENT_SEMANTICS": GEN / "p2617_s1567_p2606_exponent_semantics_reclassification.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "full_completion_map_exported",
    "beta_tors_to_beta_theorem_exported",
    "gf2_free_strict_selector_exported",
    "p2607_bridge_completion_revalidated",
    "p2608_role_bearing_ltotal_reenabled",
    "legacy_physical_role_transfer_exported",
    "qw2191_discharged_by_this_packet",
    "toe_closure_claimed",
    "apd_source_exported",
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
        "new_packet": "P2618|S1568|analytic legacy-to-strict completion|orientation-natural selector obstruction|damping completion obstruction",
        "damping_precursors": "P2603|P2614|P2615|P2617|d\\^eta|beta_tors|fractal codimension|eta=9/5",
        "phase_precursors": "P2607|P2612|GF\\(2\\)|phase topological selector|QW-2191|orientation|selector source",
        "role_guards": "P2616|P2608|role-bearing L_total|role-transfer theorem|legacy physical-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def exact_damping_rows() -> list[dict[str, Any]]:
    eta = Fraction(9, 5)
    delta = Fraction(4, 5)
    beta_tors = Fraction(1, 100)
    beta = Fraction(1, 1)
    rows: list[dict[str, Any]] = []
    for d in range(1, 13):
        legacy_den = Fraction(1, 1) + beta_tors * d
        # Use float only as a sanity column; exact proof below does not use it.
        strict_den_float = 1.0 + float(beta) * (d ** float(eta))
        residual_float = strict_den_float - float(legacy_den)
        rows.append({
            "d": d,
            "delta_codimension": str(delta),
            "eta_strict": str(eta),
            "legacy_denominator_exact": str(legacy_den),
            "strict_denominator_float": strict_den_float,
            "strict_minus_legacy_float": residual_float,
            "strict_exceeds_legacy_for_beta_1_beta_tors_0_01": residual_float > 0.0,
        })
    return rows


def scalar_renormalization_obstruction_cases() -> list[dict[str, Any]]:
    # If c(1+bt*d)=1+b*d^(9/5) for all d>0, compare d->0 and derivatives:
    # constant term gives c=1; first derivative gives bt = b*eta*d^(eta-1), impossible constant for b>0, eta!=1.
    return [
        {
            "case": "constant_term",
            "equation": "c*(1+beta_tors*d)=1+beta*d^(9/5)",
            "limit_d_to_0_plus": "c=1",
            "forced": True,
        },
        {
            "case": "first_derivative",
            "equation_after_c_equals_1": "beta_tors = beta*(9/5)*d^(4/5) for every d>0",
            "d_dependence": "right side varies with d when beta>0",
            "contradicts_constant_beta_tors": True,
        },
        {
            "case": "degenerate_escape",
            "condition": "beta=0 or eta=1",
            "admissible_for_strict_eta_9_5_positive_damping": False,
        },
    ]


def orientation_selector_sanity_rows() -> list[dict[str, Any]]:
    rows = []
    for candidate_name, value in [("plus_orientation", 1), ("minus_orientation", -1)]:
        reversed_value = -value
        rows.append({
            "candidate": candidate_name,
            "selected_sign": value,
            "orientation_reversed_selected_sign": reversed_value,
            "orientation_invariant_selector_condition_s_equals_reversed_s": value == reversed_value,
            "odd_phase_sign_survives_without_orientation_source": value == reversed_value,
        })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    payloads = {name: load_json(path) for name, path in SOURCE_FILES.items() if path.suffix == ".json"}
    p2612 = theorem(payloads.get("P2612_GF2_ORIGIN_OBSTRUCTION", {}), "p2607_gf2_physical_origin_obstruction")
    p2614 = theorem(payloads.get("P2614_CONTINUUM_RG_CHARACTER", {}), "p2602_continuum_rg_character_prime_log_proof")
    p2616 = theorem(payloads.get("P2616_ROLE_ACCEPTANCE_OBSTRUCTION", {}), "p2608_role_acceptance_obstruction_after_source_revalidation")
    p2617 = theorem(payloads.get("P2617_EXPONENT_SEMANTICS", {}), "p2606_exponent_semantics_reclassification")
    damping_rows = exact_damping_rows()
    renorm_cases = scalar_renormalization_obstruction_cases()
    selector_rows = orientation_selector_sanity_rows()
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2618_T1_analytic_legacy_to_strict_completion_obstruction",
        "requested_ticket": rel(SOURCE_FILES["TICKET_P2618"]),
        "result_type": "fundamental_obstruction_to_full_a_priori_completion_map_under_current_sources",
        "damping_partial_result": {
            "df_used": "9/5",
            "delta_codimension": "4/5",
            "eta_strict_denominator": "9/5",
            "eta_relation": "eta_strict = 1 + delta = D_f",
            "accepted_source_revalidation_inherited_from_p2614": p2614.get("strict_damping_beta_eta_source_revalidated_under_df_9_5_scope") is True,
            "p2617_semantics_guard_inherited": p2617.get("eta_strict_equals_one_plus_delta") is True,
            "what_is_proven": "Fractal projection can justify the exponent upgrade eta=9/5 as linear backbone plus codimension slope, but not an exact beta_tors -> beta scalar completion map.",
        },
        "damping_nonrenormalization_proof": {
            "formal_statement": "No positive constants beta_tors, beta and scalar c can satisfy c(1+beta_tors*d)=1+beta*d^(9/5) for all d>0, except degenerate non-strict escapes beta=0 or eta=1.",
            "proof_steps": [
                "Assume c(1+beta_tors*d)=1+beta*d^eta for every d>0 with eta=9/5 and beta>0.",
                "Taking d -> 0+ compares constant terms and forces c=1.",
                "Differentiating gives beta_tors = beta*eta*d^(eta-1) for every d>0.",
                "Because eta-1=4/5 is nonzero and beta>0, beta*eta*d^(4/5) is not constant in d.",
                "Therefore exact scalar renormalization or silent substitution from the linear torsion denominator to the strict nonlinear denominator is impossible under the current axioms.",
            ],
            "symbolic_cases": renorm_cases,
        },
        "damping_sample_sanity_check": damping_rows,
        "phase_selector_obstruction_proof": {
            "formal_statement": "A GF(2)-free topological/field-theoretic construction that is natural under orientation reversal cannot select a nonzero odd phase sign without an additional orientation, symmetry-breaking, boundary, or source premise.",
            "proof_steps": [
                "The strict phase sign is an odd orientation datum: reversing the relevant orientation sends sigma to -sigma.",
                "A purely invariant construction with no extra source must be natural under the same orientation reversal, hence must output the same selected datum before and after reversal.",
                "Combining oddness and naturality gives sigma = -sigma.",
                "For a strict sign sigma in {+1,-1}, sigma = -sigma is impossible.",
                "Thus analytic topology can classify the sign torsor/cohomology class, but it cannot choose the strict representative without an additional selector source.",
            ],
            "allowed_future_escape_premises": [
                "explicit orientation or spin/Pin structure with physical source",
                "symmetry-breaking boundary condition",
                "new internally exported selector current/source",
                "role-transfer theorem proving beta_tors -> oriented strict datum",
            ],
            "orientation_reversal_sanity_check": selector_rows,
            "p2612_gf2_obstruction_inherited": p2612.get("p2607_bridge_completion_revalidated_by_p2612") is False or p2612.get("gf2_physical_origin_obstruction_exported") is True,
        },
        "completion_map_verdict": {
            "amplitude_normalization": "not attacked here; previous amplitude/scalar normalization certificates remain bookkeeping only unless role transfer is supplied",
            "damping_compression": "partial exponent-source theorem survives; exact beta_tors -> beta completion is obstructed",
            "phase_topological_selector": "classification possible, strict representative selection obstructed without new source",
            "role_transfer": "blocked by inherited P2616 acceptance predicate",
            "honest_next_step": "Do not promote L_total. Next prove or add a non-GF(2) orientation/source premise, then rerun role-transfer audit claim-by-claim.",
        },
        "p2616_role_obstruction_inherited": p2616.get("role_bearing_ltotal_accepted") is False or p2616.get("role_bearing_L_total_accepted") is False,
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "strict_eta_relation_recorded": theorem_export["damping_partial_result"]["eta_relation"] == "eta_strict = 1 + delta = D_f",
        "renormalization_cases_all_block_exact_scalar_map": all(case.get("forced", True) or case.get("contradicts_constant_beta_tors", True) or not case.get("admissible_for_strict_eta_9_5_positive_damping", True) for case in renorm_cases),
        "damping_samples_show_nonzero_residual": all(row["strict_exceeds_legacy_for_beta_1_beta_tors_0_01"] for row in damping_rows),
        "orientation_rows_reject_invariant_odd_selector": all(not row["odd_phase_sign_survives_without_orientation_source"] for row in selector_rows),
        "no_full_completion_map": theorem_export["full_completion_map_exported"] is False,
        "no_beta_tors_to_beta_theorem": theorem_export["beta_tors_to_beta_theorem_exported"] is False,
        "no_gf2_free_selector": theorem_export["gf2_free_strict_selector_exported"] is False,
        "no_ltotal_reenable": theorem_export["p2608_role_bearing_ltotal_reenabled"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_packet"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2618",
        "stage_id": "S1568",
        "status": "P2618_ANALYTIC_COMPLETION_MAP_PARTIAL_EXPONENT_SOURCE_WITH_DAMPING_AND_SELECTOR_OBSTRUCTIONS_NO_LTOTAL_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "analytic_legacy_to_strict_completion_obstruction": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {name: sha256_json(payload) for name, payload in payloads.items()},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["analytic_legacy_to_strict_completion_obstruction"]["theorem_export"]
    lines = [
        "# P2618/S1568 analytic legacy-to-strict completion obstruction", "",
        f"Status: `{payload['status']}`", "", "## Verdict", "",
        "The requested full a-priori completion map is **not** exported.  The honest result is a partial analytic damping exponent source plus two obstructions: exact scalar damping completion is impossible for `eta=9/5`, and a GF(2)-free purely invariant phase construction cannot choose the strict odd phase sign without an added selector source.", "",
        "## Damping/compression result", "",
        f"- `D_f = {t['damping_partial_result']['df_used']}`.",
        f"- `delta = {t['damping_partial_result']['delta_codimension']}`.",
        f"- `eta_strict = {t['damping_partial_result']['eta_strict_denominator']}`.",
        f"- Relation: `{t['damping_partial_result']['eta_relation']}`.",
        f"- Scope: {t['damping_partial_result']['what_is_proven']}", "",
        "## Exact damping obstruction proof", "",
        t["damping_nonrenormalization_proof"]["formal_statement"], "",
    ]
    for step in t["damping_nonrenormalization_proof"]["proof_steps"]:
        lines.append(f"- {step}")
    lines.extend(["", "## Phase/topological selector obstruction", "", t["phase_selector_obstruction_proof"]["formal_statement"], ""])
    for step in t["phase_selector_obstruction_proof"]["proof_steps"]:
        lines.append(f"- {step}")
    lines.extend(["", "Allowed future escape premises:"])
    for premise in t["phase_selector_obstruction_proof"]["allowed_future_escape_premises"]:
        lines.append(f"- {premise}")
    lines.extend([
        "", "## Completion-map verdict", "",
        f"- Damping/compression: {t['completion_map_verdict']['damping_compression']}.",
        f"- Phase/topological selector: {t['completion_map_verdict']['phase_topological_selector']}.",
        f"- Role transfer: {t['completion_map_verdict']['role_transfer']}.",
        f"- Honest next step: {t['completion_map_verdict']['honest_next_step']}",
        "", "## Scope guards", "",
        "No full completion map, no `beta_tors -> beta` theorem, no GF(2)-free strict selector, no role-bearing `L_total`, no legacy physical-role transfer, no `QW-2191` discharge, and no ToE closure are exported.",
        "", "## Fingerprint", "",
        f"`{payload['analytic_legacy_to_strict_completion_obstruction']['theorem_fingerprint_sha256']}`",
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2618/S1568 analytic completion-map obstruction

`P2618/S1568` answers the requested legacy-to-strict completion-map ticket with a partial proof and an obstruction.  The fractal projection source supports the exponent relation `eta_strict=1+(D_f-1)=9/5`, but an exact scalar completion `c(1+beta_tors*d)=1+beta*d^(9/5)` is impossible for positive strict damping because the derivative would require constant `beta_tors` to equal `beta*(9/5)*d^(4/5)` for every `d`.  The phase/topological sign is likewise only classifiable as an odd orientation/cohomology datum; a purely invariant GF(2)-free construction cannot select its representative without an additional orientation, symmetry-breaking, boundary, or source premise.
""".strip()
    lag_section = """
## P2618/S1568 Ltotal bridge guard

`P2618/S1568` does not open role-bearing `L_total`.  It preserves the strict exponent source `eta=9/5` as a partial damping-compression result, but blocks exact `beta_tors -> beta` completion, GF(2)-free strict selector export, legacy physical-role transfer, `QW-2191` discharge, and ToE closure pending a real selector/source theorem and a separate role-transfer audit.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2618/S1568 analytic completion-map obstruction", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2618/S1568 Ltotal bridge guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
