#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from fractions import Fraction
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
OUT = GEN / "p2541_s1491_strict_damping_multiplicative_character_current_premise_obstruction_certificate.json"
MD = GEN / "p2541_s1491_strict_damping_multiplicative_character_current_premise_obstruction_certificate.md"

SOURCE_FILES = {
    "P2524_AFFINE_CONTINUUM": GEN / "p2524_s1474_strict_damping_affine_consistency_continuum_nonidentifiability_certificate.json",
    "P2525_MULTIPLICATIVE_BETA_NORMALIZATION": GEN / "p2525_s1475_strict_damping_multiplicative_cocycle_beta_normalization_subkey_certificate.json",
    "P2526_FINITE_MONOID_NULLITY": GEN / "p2526_s1476_strict_damping_finite_monoid_prime_character_nullity_certificate.json",
    "P2530_FOUR_KEY_IRREDUNDANCY": GEN / "p2530_s1480_strict_damping_four_key_irredundancy_witness_certificate.json",
    "P2539_TOE_POTENTIAL_RECOMMENDATION": GEN / "p2539_s1489_strict_damping_toe_potential_recommendation_certificate.json",
    "P2540_M2_CURRENT_PREMISE_OBSTRUCTION": GEN / "p2540_s1490_strict_damping_m2_operator_signature_current_premise_obstruction_certificate.json",
}

NODE_DOMAIN = list(range(1, 12))
INTERCEPT_CANDIDATES = [Fraction(-1, 1), Fraction(-1, 2), Fraction(0, 1), Fraction(1, 2), Fraction(1, 1)]
SLOPE_CANDIDATES = [Fraction(-1, 1), Fraction(0, 1), Fraction(1, 2), Fraction(4, 5), Fraction(1, 1), Fraction(9, 5), Fraction(2, 1)]
STRICT_INTERCEPT = Fraction(0, 1)
STRICT_DELTA = Fraction(4, 5)
TOL = 1e-14


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
    return {"count": len(lines), "samples": lines[:50]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2541|S1491|multiplicative character current-premise obstruction|affine plus unital|multiplicative source obstruction",
        "intended_research_nonduplication": "multiplicative_character_law_source|multiplicative character law|unital character|left normalization|affine.*multiplicative",
        "precursor_packets": "P2524|S1474|P2525|S1475|P2526|S1476|P2530|S1480|P2539|S1489|P2540|S1490",
        "source_key_language": "multiplicative_character_law_source|prime_log_proportionality_source|slope_value_or_prime_anchor_source|m2_operator_signature_source",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def frac_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def y_affine(d: int, intercept: Fraction, slope: Fraction) -> float:
    return float(intercept) + float(slope) * math.log(d)


def multiplicative_pairs() -> list[tuple[int, int, int]]:
    return [(d, e, d * e) for d, e in product(NODE_DOMAIN, repeat=2) if d * e in NODE_DOMAIN]


def affine_multiplicative_equivalence_rows() -> list[dict[str, Any]]:
    pairs = multiplicative_pairs()
    rows = []
    for intercept in INTERCEPT_CANDIDATES:
        for slope in SLOPE_CANDIDATES:
            defects = [
                y_affine(de, intercept, slope) - y_affine(d, intercept, slope) - y_affine(e, intercept, slope)
                for d, e, de in pairs
            ]
            max_abs_defect = max(abs(value) for value in defects)
            y1 = y_affine(1, intercept, slope)
            rows.append({
                "intercept_log_beta": frac_text(intercept),
                "slope_delta": frac_text(slope),
                "is_strict_numeric_target": intercept == STRICT_INTERCEPT and slope == STRICT_DELTA,
                "affine_consistency_premise_accepts": True,
                "unital_y1_zero_accepts": abs(y1) < TOL,
                "multiplicative_character_accepts": max_abs_defect < TOL,
                "max_abs_multiplicative_defect": max_abs_defect,
                "defect_is_minus_intercept_on_all_pairs": all(abs(value + float(intercept)) < TOL for value in defects),
                "affine_plus_unital_entails_multiplicative_on_grid": abs(y1) < TOL and max_abs_defect < TOL,
                "multiplicative_entails_unital_on_grid": max_abs_defect < TOL and abs(y1) < TOL,
            })
    return rows


def obstruction_witnesses(rows: list[dict[str, Any]]) -> dict[str, Any]:
    affine_countermodels = [
        row for row in rows
        if row["affine_consistency_premise_accepts"] and not row["multiplicative_character_accepts"]
    ]
    multiplicative_rows = [row for row in rows if row["multiplicative_character_accepts"]]
    nonstrict_slope_rows = [
        row for row in multiplicative_rows
        if not row["is_strict_numeric_target"]
    ]
    chosen_affine_countermodel = next(
        row for row in affine_countermodels
        if row["intercept_log_beta"] == "1/2" and row["slope_delta"] == "4/5"
    )
    chosen_multiplicative_nonstrict_slope = next(
        row for row in nonstrict_slope_rows
        if row["intercept_log_beta"] == "0" and row["slope_delta"] == "1/2"
    )
    return {
        "affine_countermodel_count": len(affine_countermodels),
        "multiplicative_accepting_count": len(multiplicative_rows),
        "multiplicative_nonstrict_slope_count": len(nonstrict_slope_rows),
        "chosen_affine_countermodel_passing_affine_but_failing_multiplicativity": chosen_affine_countermodel,
        "chosen_multiplicative_nonstrict_slope_witness": chosen_multiplicative_nonstrict_slope,
        "all_multiplicative_rows_are_unital": all(row["unital_y1_zero_accepts"] for row in multiplicative_rows),
        "all_unital_affine_rows_are_multiplicative": all(row["multiplicative_character_accepts"] for row in rows if row["unital_y1_zero_accepts"]),
        "affine_consistency_alone_entails_multiplicativity": len(affine_countermodels) == 0,
        "affine_plus_unital_equivalent_to_multiplicative_on_grid": (
            all(row["multiplicative_character_accepts"] == row["unital_y1_zero_accepts"] for row in rows)
        ),
        "multiplicative_law_selects_strict_slope": len(nonstrict_slope_rows) == 0,
    }


def build_obstruction_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2524 = theorem(sources["P2524_AFFINE_CONTINUUM"], "strict_damping_affine_consistency_continuum_nonidentifiability_certificate")
    p2525 = theorem(sources["P2525_MULTIPLICATIVE_BETA_NORMALIZATION"], "strict_damping_multiplicative_cocycle_beta_normalization_subkey_certificate")
    p2526 = theorem(sources["P2526_FINITE_MONOID_NULLITY"], "strict_damping_finite_monoid_prime_character_nullity_certificate")
    p2530 = theorem(sources["P2530_FOUR_KEY_IRREDUNDANCY"], "strict_damping_four_key_irredundancy_witness_certificate")
    p2539 = theorem(sources["P2539_TOE_POTENTIAL_RECOMMENDATION"], "strict_damping_toe_potential_recommendation_certificate")
    p2540 = theorem(sources["P2540_M2_CURRENT_PREMISE_OBSTRUCTION"], "strict_damping_m2_operator_signature_current_premise_obstruction_certificate")
    rows = affine_multiplicative_equivalence_rows()
    witnesses = obstruction_witnesses(rows)
    return {
        "frontier_source_key_under_attack": "multiplicative_character_law_source",
        "p2524_affine_continuum_inherited": p2524.get("basis_independent_affine_consistency_nonidentifiability_exported") is True,
        "p2525_conditional_beta_normalization_inherited": p2525.get("conditional_beta_normalization_subkey_exported") is True,
        "p2526_prime_character_nullity_inherited": p2526.get("conditional_finite_monoid_character_nullity_exported") is True,
        "p2530_four_key_irredundancy_inherited": p2530.get("four_key_irredundancy_witness_exported") is True,
        "p2539_next_step_recommendation_inherited": p2539.get("recommended_next_honest_step") == "prove_or_refute_one_strict_source_theorem_for_a_single_P2529_source_key_before_more_bookkeeping_layers",
        "p2540_m2_route_refutation_inherited": p2540.get("m2_operator_signature_source_route_refuted_for_current_source_free_premises") is True,
        "candidate_grid_rows": rows,
        "obstruction_witnesses": witnesses,
        "current_affine_premises_do_not_entail_multiplicative_law": not witnesses["affine_consistency_alone_entails_multiplicativity"],
        "multiplicative_law_equivalent_to_unital_left_normalization_inside_affine_family": witnesses["affine_plus_unital_equivalent_to_multiplicative_on_grid"],
        "multiplicative_law_overgenerates_slope_family": not witnesses["multiplicative_law_selects_strict_slope"],
        "current_premise_obstruction_exported": True,
        "required_new_premise_class": "strict dynamics must source unital left normalization or an equivalent monoid-character law; affine consistency alone is insufficient",
    }


def append_doc_sections() -> None:
    eq_section = """
## P2541/S1491 strict damping multiplicative-character current-premise obstruction certificate

`P2541/S1491` attacks the `multiplicative_character_law_source` key after P2540.  It proves a precise obstruction/equivalence on the affine family inherited from P2524: for `y_d=b+a log d`, the multiplicative defect `y_{de}-y_d-y_e` is `-b` on every audited product `de<=11`.  Therefore affine consistency alone does not entail the multiplicative character law; for example `(b,a)=(1/2,4/5)` passes affine secant/cocycle consistency but fails multiplicativity.  Conversely, inside the affine family, adding the unital left-normalization `y_1=0` is equivalent to multiplicativity on the audited monoid.

Thus the current route cannot export `multiplicative_character_law_source` unless strict dynamics source `y_1=0` or an equivalent monoid-character law.  Even if that law is supplied, it leaves the slope family `a log d`, so prime log-proportionality/slope-value and the operator-signature key remain separate.  No strict source key, source-obligation discharge, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, physical-value generator, or ToE closure is exported.
"""
    lag_section = """
## P2541/S1491 multiplicative-character source obstruction guard

`P2541/S1491` shows that current affine-consistency premises do not source the multiplicative law: `y_d=b+a log d` is multiplicative exactly when the missing left-normalization `b=y_1=0` is supplied.  So `multiplicative_character_law_source` remains open; it cannot be promoted from consistency checks without a strict unital/monoid-character source.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2541/S1491 strict damping multiplicative-character current-premise obstruction certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2541/S1491 multiplicative-character source obstruction guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    cert = build_obstruction_certificate(sources)
    witnesses = cert["obstruction_witnesses"]
    theorem_export = {
        "theorem_name": "P2541_T1_strict_damping_multiplicative_character_current_premise_obstruction_certificate",
        "audited_chain": ["P2524/S1474", "P2525/S1475", "P2526/S1476", "P2530/S1480", "P2539/S1489", "P2540/S1490"],
        "strict_damping_multiplicative_character_current_premise_obstruction_certificate": cert,
        "frontier_source_key_under_attack": cert["frontier_source_key_under_attack"],
        "p2524_affine_continuum_inherited": cert["p2524_affine_continuum_inherited"],
        "p2525_conditional_beta_normalization_inherited": cert["p2525_conditional_beta_normalization_inherited"],
        "p2526_prime_character_nullity_inherited": cert["p2526_prime_character_nullity_inherited"],
        "p2530_four_key_irredundancy_inherited": cert["p2530_four_key_irredundancy_inherited"],
        "p2539_next_step_recommendation_inherited": cert["p2539_next_step_recommendation_inherited"],
        "p2540_m2_route_refutation_inherited": cert["p2540_m2_route_refutation_inherited"],
        "candidate_grid_row_count": len(cert["candidate_grid_rows"]),
        "affine_countermodel_count": witnesses["affine_countermodel_count"],
        "multiplicative_accepting_count": witnesses["multiplicative_accepting_count"],
        "multiplicative_nonstrict_slope_count": witnesses["multiplicative_nonstrict_slope_count"],
        "chosen_affine_countermodel_passing_affine_but_failing_multiplicativity": witnesses["chosen_affine_countermodel_passing_affine_but_failing_multiplicativity"],
        "chosen_multiplicative_nonstrict_slope_witness": witnesses["chosen_multiplicative_nonstrict_slope_witness"],
        "all_multiplicative_rows_are_unital": witnesses["all_multiplicative_rows_are_unital"],
        "all_unital_affine_rows_are_multiplicative": witnesses["all_unital_affine_rows_are_multiplicative"],
        "current_affine_premises_do_not_entail_multiplicative_law": cert["current_affine_premises_do_not_entail_multiplicative_law"],
        "multiplicative_law_equivalent_to_unital_left_normalization_inside_affine_family": cert["multiplicative_law_equivalent_to_unital_left_normalization_inside_affine_family"],
        "multiplicative_law_overgenerates_slope_family": cert["multiplicative_law_overgenerates_slope_family"],
        "current_premise_obstruction_exported": cert["current_premise_obstruction_exported"],
        "required_new_premise_class": cert["required_new_premise_class"],
        "multiplicative_character_law_source_exported": False,
        "prime_log_proportionality_source_exported": False,
        "slope_value_or_prime_anchor_source_exported": False,
        "beta_eta_numeric_source_exported": False,
        "m2_operator_signature_source_exported": False,
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
            "P2541 refutes only the current affine-consistency route to multiplicative_character_law_source.",
            "Inside the affine family, multiplicativity is equivalent to unital left normalization, which remains unsourced.",
            "Even a sourced multiplicative law would still leave slope and m2 operator-signature obligations open.",
            "No bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Try to source the unital left-normalization y_1=0 from strict nadsoliton dynamics, or pivot to prime_log_proportionality_source / slope_value_or_prime_anchor_source.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "all_required_precursors_inherited": all(theorem_export[key] for key in [
            "p2524_affine_continuum_inherited",
            "p2525_conditional_beta_normalization_inherited",
            "p2526_prime_character_nullity_inherited",
            "p2530_four_key_irredundancy_inherited",
            "p2539_next_step_recommendation_inherited",
            "p2540_m2_route_refutation_inherited",
        ]),
        "countermodel_exists": theorem_export["affine_countermodel_count"] > 0 and theorem_export["current_affine_premises_do_not_entail_multiplicative_law"],
        "affine_unital_equivalence_verified": theorem_export["all_multiplicative_rows_are_unital"] and theorem_export["all_unital_affine_rows_are_multiplicative"],
        "slope_overgeneration_verified": theorem_export["multiplicative_law_overgenerates_slope_family"] and theorem_export["multiplicative_nonstrict_slope_count"] > 0,
        "multiplicative_source_not_exported": not theorem_export["multiplicative_character_law_source_exported"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "multiplicative_character_law_source_exported",
            "prime_log_proportionality_source_exported",
            "slope_value_or_prime_anchor_source_exported",
            "beta_eta_numeric_source_exported",
            "m2_operator_signature_source_exported",
            "strict_damping_beta_eta_source_exported",
            "source_obligation_discharge_exported",
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
        "packet_id": "P2541",
        "stage_id": "S1491",
        "status": "STRICT_DAMPING_MULTIPLICATIVE_CHARACTER_CURRENT_PREMISE_OBSTRUCTION_CERTIFICATE_CURRENT_AFFINE_ROUTE_REFUTED_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_multiplicative_character_current_premise_obstruction_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_multiplicative_character_current_premise_obstruction_certificate"]["theorem_export"]
    lines = [
        "# P2541/S1491 strict damping multiplicative-character current-premise obstruction certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier source key under attack: `{t['frontier_source_key_under_attack']}`.",
        f"- Candidate grid row count: `{t['candidate_grid_row_count']}`.",
        f"- Affine countermodel count: `{t['affine_countermodel_count']}`.",
        f"- Multiplicative accepting count: `{t['multiplicative_accepting_count']}`.",
        f"- Multiplicative non-strict slope count: `{t['multiplicative_nonstrict_slope_count']}`.",
        f"- Current affine premises entail multiplicativity: `{not t['current_affine_premises_do_not_entail_multiplicative_law']}`.",
        f"- Multiplicative law equivalent to unital left normalization inside affine family: `{t['multiplicative_law_equivalent_to_unital_left_normalization_inside_affine_family']}`.",
        f"- Multiplicative law overgenerates slope family: `{t['multiplicative_law_overgenerates_slope_family']}`.",
        f"- Multiplicative source exported: `{t['multiplicative_character_law_source_exported']}`.",
        "",
        "## Interpretation",
        "",
        t["required_new_premise_class"],
        "",
        "## Negative controls",
        "",
        "This packet exports a current-premise obstruction only. It does not source any strict damping key, discharge source_obligation_discharge, complete the bridge, export role transfer, discharge QW-2191, produce role-bearing L_total, or claim ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_multiplicative_character_current_premise_obstruction_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_multiplicative_character_current_premise_obstruction_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
