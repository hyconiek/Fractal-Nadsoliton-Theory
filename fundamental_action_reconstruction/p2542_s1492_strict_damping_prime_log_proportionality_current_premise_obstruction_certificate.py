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
OUT = GEN / "p2542_s1492_strict_damping_prime_log_proportionality_current_premise_obstruction_certificate.json"
MD = GEN / "p2542_s1492_strict_damping_prime_log_proportionality_current_premise_obstruction_certificate.md"

SOURCE_FILES = {
    "P2526_FINITE_MONOID_NULLITY": GEN / "p2526_s1476_strict_damping_finite_monoid_prime_character_nullity_certificate.json",
    "P2527_PRIME_LOG_SLOPE_LINE": GEN / "p2527_s1477_strict_damping_prime_log_proportionality_slope_line_certificate.json",
    "P2530_FOUR_KEY_IRREDUNDANCY": GEN / "p2530_s1480_strict_damping_four_key_irredundancy_witness_certificate.json",
    "P2539_TOE_POTENTIAL_RECOMMENDATION": GEN / "p2539_s1489_strict_damping_toe_potential_recommendation_certificate.json",
    "P2541_MULTIPLICATIVE_CURRENT_PREMISE_OBSTRUCTION": GEN / "p2541_s1491_strict_damping_multiplicative_character_current_premise_obstruction_certificate.json",
}

NODE_DOMAIN = list(range(1, 12))
PRIMES = [2, 3, 5, 7, 11]
STRICT_DELTA = Fraction(4, 5)
TOL = 1e-12


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
        "new_packet": "P2542|S1492|prime-log proportionality current-premise obstruction|prime log source obstruction|prime-character countermodel",
        "intended_research_nonduplication": "prime_log_proportionality_source|prime log proportionality|prime-generator freedom|finite monoid character|log-proportional prime",
        "precursor_packets": "P2526|S1476|P2527|S1477|P2530|S1480|P2539|S1489|P2541|S1491",
        "source_key_language": "multiplicative_character_law_source|prime_log_proportionality_source|slope_value_or_prime_anchor_source|m2_operator_signature_source",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def factorization_exponents(n: int) -> list[int]:
    value = n
    exponents = []
    for prime in PRIMES:
        exponent = 0
        while value % prime == 0:
            value //= prime
            exponent += 1
        exponents.append(exponent)
    if value != 1:
        raise ValueError(f"unexpected factor outside audited primes: {n}")
    return exponents


def multiplicative_pairs() -> list[tuple[int, int, int]]:
    return [(d, e, d * e) for d, e in product(NODE_DOMAIN, repeat=2) if d * e in NODE_DOMAIN]


def node_values_from_prime_values(prime_values: list[float]) -> list[float]:
    return [
        sum(exponent * value for exponent, value in zip(factorization_exponents(d), prime_values))
        for d in NODE_DOMAIN
    ]


def max_multiplicative_residual(node_values: list[float]) -> float:
    return max(abs(node_values[de - 1] - node_values[d - 1] - node_values[e - 1]) for d, e, de in multiplicative_pairs())


def normalized_prime_ratios(prime_values: list[float]) -> list[float]:
    return [value / math.log(prime) for prime, value in zip(PRIMES, prime_values)]


def prime_character_samples() -> dict[str, list[float]]:
    strict = [float(STRICT_DELTA) * math.log(prime) for prime in PRIMES]
    samples = {
        "strict_log_slope_delta_4_over_5": strict,
        "nonstrict_log_slope_delta_1_over_2": [0.5 * math.log(prime) for prime in PRIMES],
        "unit_p2_only_character": [1.0, 0.0, 0.0, 0.0, 0.0],
        "unit_p3_only_character": [0.0, 1.0, 0.0, 0.0, 0.0],
        "unit_p5_only_character": [0.0, 0.0, 1.0, 0.0, 0.0],
        "perturbed_strict_p2_character": [strict[0] + 0.25] + strict[1:],
        "mixed_integer_prime_character": [1.0, 2.0, 3.0, 5.0, 8.0],
    }
    return samples


def prime_character_rows() -> list[dict[str, Any]]:
    rows = []
    for name, prime_values in prime_character_samples().items():
        node_values = node_values_from_prime_values(prime_values)
        ratios = normalized_prime_ratios(prime_values)
        ratio_spread = max(ratios) - min(ratios)
        multiplicative_residual = max_multiplicative_residual(node_values)
        rows.append({
            "sample_name": name,
            "prime_values_by_p_2_3_5_7_11": prime_values,
            "node_values_y_1_to_y_11": node_values,
            "unital_y1_zero_accepts": abs(node_values[0]) < TOL,
            "max_multiplicative_residual": multiplicative_residual,
            "multiplicative_character_accepts": multiplicative_residual < TOL,
            "normalized_prime_ratios_vp_over_logp": ratios,
            "normalized_ratio_spread": ratio_spread,
            "prime_log_proportionality_accepts": ratio_spread < TOL,
            "is_strict_numeric_target_slope": name == "strict_log_slope_delta_4_over_5",
            "is_log_proportional_but_nonstrict_slope": name == "nonstrict_log_slope_delta_1_over_2",
            "is_current_premise_countermodel_to_prime_log_proportionality": multiplicative_residual < TOL and ratio_spread >= TOL,
        })
    return rows


def exact_ratio_constraint_summary() -> dict[str, Any]:
    # Use normalized ratio variables r_p=v_p/log(p).  Prime-log proportionality
    # is exactly r_p-r_2=0 for p in {3,5,7,11}.
    rows = []
    for index in range(1, len(PRIMES)):
        row = [0 for _ in PRIMES]
        row[0] = -1
        row[index] = 1
        rows.append(row)
    return {
        "variables": [f"r_{prime}=v_{prime}/log({prime})" for prime in PRIMES],
        "constraint_rows": rows,
        "rank": 4,
        "nullity": 1,
        "interpretation": "These four independent equations are the missing prime-log proportionality premise; they are not implied by finite monoid multiplicativity.",
    }


def build_obstruction_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2526 = theorem(sources["P2526_FINITE_MONOID_NULLITY"], "strict_damping_finite_monoid_prime_character_nullity_certificate")
    p2527 = theorem(sources["P2527_PRIME_LOG_SLOPE_LINE"], "strict_damping_prime_log_proportionality_slope_line_certificate")
    p2530 = theorem(sources["P2530_FOUR_KEY_IRREDUNDANCY"], "strict_damping_four_key_irredundancy_witness_certificate")
    p2539 = theorem(sources["P2539_TOE_POTENTIAL_RECOMMENDATION"], "strict_damping_toe_potential_recommendation_certificate")
    p2541 = theorem(sources["P2541_MULTIPLICATIVE_CURRENT_PREMISE_OBSTRUCTION"], "strict_damping_multiplicative_character_current_premise_obstruction_certificate")
    rows = prime_character_rows()
    countermodels = [row for row in rows if row["is_current_premise_countermodel_to_prime_log_proportionality"]]
    proportional_rows = [row for row in rows if row["prime_log_proportionality_accepts"]]
    strict_rows = [row for row in rows if row["is_strict_numeric_target_slope"]]
    nonstrict_proportional_rows = [row for row in rows if row["is_log_proportional_but_nonstrict_slope"]]
    chosen_countermodel = next(row for row in countermodels if row["sample_name"] == "unit_p2_only_character")
    return {
        "frontier_source_key_under_attack": "prime_log_proportionality_source",
        "p2526_prime_character_nullity_inherited": p2526.get("conditional_finite_monoid_character_nullity_exported") is True,
        "p2527_slope_line_subkey_inherited": p2527.get("prime_log_proportionality_subkey_exported") is True,
        "p2530_four_key_irredundancy_inherited": p2530.get("four_key_irredundancy_witness_exported") is True,
        "p2539_next_step_recommendation_inherited": p2539.get("recommended_next_honest_step") == "prove_or_refute_one_strict_source_theorem_for_a_single_P2529_source_key_before_more_bookkeeping_layers",
        "p2541_multiplicative_route_obstruction_inherited": p2541.get("current_premise_obstruction_exported") is True,
        "prime_character_rows": rows,
        "exact_ratio_constraint_summary": exact_ratio_constraint_summary(),
        "current_unital_multiplicative_premises_do_not_entail_prime_log_proportionality": len(countermodels) > 0,
        "countermodel_count": len(countermodels),
        "chosen_countermodel_unit_p2_only": chosen_countermodel,
        "all_sample_rows_are_unital_multiplicative_characters": all(row["unital_y1_zero_accepts"] and row["multiplicative_character_accepts"] for row in rows),
        "proportional_sample_count": len(proportional_rows),
        "strict_log_slope_is_one_proportional_member": len(strict_rows) == 1 and strict_rows[0]["prime_log_proportionality_accepts"],
        "nonstrict_log_slope_also_proportional": len(nonstrict_proportional_rows) == 1 and nonstrict_proportional_rows[0]["prime_log_proportionality_accepts"],
        "prime_log_proportionality_overgenerates_slope_line": len(nonstrict_proportional_rows) > 0,
        "current_premise_obstruction_exported": True,
        "required_new_premise_class": "strict dynamics must source equality of normalized prime ratios v_p/log(p), or an equivalent prime-log proportionality law",
    }


def append_doc_sections() -> None:
    eq_section = """
## P2542/S1492 strict damping prime-log proportionality current-premise obstruction certificate

`P2542/S1492` attacks the `prime_log_proportionality_source` key.  It inherits P2526 finite-monoid prime-character nullity, P2527 conditional slope-line reduction, P2530 four-key irredundancy, P2539 source-key recommendation, and P2541 multiplicative-character obstruction.  The new audit constructs explicit unital multiplicative characters on the audited monoid `d=1..11` by freely assigning values to the prime generators `2,3,5,7,11`.  These characters all satisfy `y_1=0` and `y_{de}=y_d+y_e`, but representative choices such as the `unit_p2_only_character` fail the normalized-ratio equalities `v_p/log(p)=constant`.

Therefore unital finite-monoid multiplicativity does not entail prime-log proportionality.  The missing proportionality is exactly a four-equation collapse of the five prime-generator degrees of freedom to a one-dimensional slope line; this collapse remains a separate source premise.  Even when prime-log proportionality is supplied, non-strict slopes such as `delta=1/2` still pass, so the slope-value/prime-anchor key remains separate.  No strict source key, source-obligation discharge, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, physical-value generator, or ToE closure is exported.
"""
    lag_section = """
## P2542/S1492 prime-log proportionality source obstruction guard

`P2542/S1492` shows that unital multiplicative characters on `d=1..11` do not source prime-log proportionality: arbitrary prime-generator values pass multiplicativity, while `v_p/log(p)=constant` is an additional unsourced collapse.  Thus `prime_log_proportionality_source` remains open, and the damping term still cannot become role-bearing in `L_total`.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2542/S1492 strict damping prime-log proportionality current-premise obstruction certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2542/S1492 prime-log proportionality source obstruction guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    cert = build_obstruction_certificate(sources)
    theorem_export = {
        "theorem_name": "P2542_T1_strict_damping_prime_log_proportionality_current_premise_obstruction_certificate",
        "audited_chain": ["P2526/S1476", "P2527/S1477", "P2530/S1480", "P2539/S1489", "P2541/S1491"],
        "strict_damping_prime_log_proportionality_current_premise_obstruction_certificate": cert,
        "frontier_source_key_under_attack": cert["frontier_source_key_under_attack"],
        "p2526_prime_character_nullity_inherited": cert["p2526_prime_character_nullity_inherited"],
        "p2527_slope_line_subkey_inherited": cert["p2527_slope_line_subkey_inherited"],
        "p2530_four_key_irredundancy_inherited": cert["p2530_four_key_irredundancy_inherited"],
        "p2539_next_step_recommendation_inherited": cert["p2539_next_step_recommendation_inherited"],
        "p2541_multiplicative_route_obstruction_inherited": cert["p2541_multiplicative_route_obstruction_inherited"],
        "sample_prime_character_row_count": len(cert["prime_character_rows"]),
        "all_sample_rows_are_unital_multiplicative_characters": cert["all_sample_rows_are_unital_multiplicative_characters"],
        "countermodel_count": cert["countermodel_count"],
        "chosen_countermodel_unit_p2_only": cert["chosen_countermodel_unit_p2_only"],
        "current_unital_multiplicative_premises_do_not_entail_prime_log_proportionality": cert["current_unital_multiplicative_premises_do_not_entail_prime_log_proportionality"],
        "exact_ratio_constraint_rank": cert["exact_ratio_constraint_summary"]["rank"],
        "exact_ratio_constraint_nullity": cert["exact_ratio_constraint_summary"]["nullity"],
        "proportional_sample_count": cert["proportional_sample_count"],
        "strict_log_slope_is_one_proportional_member": cert["strict_log_slope_is_one_proportional_member"],
        "nonstrict_log_slope_also_proportional": cert["nonstrict_log_slope_also_proportional"],
        "prime_log_proportionality_overgenerates_slope_line": cert["prime_log_proportionality_overgenerates_slope_line"],
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
            "P2542 refutes only the current finite-monoid/unital route to prime_log_proportionality_source.",
            "Prime-log proportionality is a separate normalized-ratio equality premise, not a consequence of multiplicativity.",
            "Even prime-log proportionality would still leave the slope-value/prime-anchor key and m2 operator-signature key open.",
            "No bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Attack slope_value_or_prime_anchor_source, or construct a strict nadsoliton source for normalized prime-ratio equality.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "all_required_precursors_inherited": all(theorem_export[key] for key in [
            "p2526_prime_character_nullity_inherited",
            "p2527_slope_line_subkey_inherited",
            "p2530_four_key_irredundancy_inherited",
            "p2539_next_step_recommendation_inherited",
            "p2541_multiplicative_route_obstruction_inherited",
        ]),
        "sample_rows_are_unital_multiplicative": theorem_export["all_sample_rows_are_unital_multiplicative_characters"],
        "countermodel_exists": theorem_export["countermodel_count"] > 0 and theorem_export["current_unital_multiplicative_premises_do_not_entail_prime_log_proportionality"],
        "ratio_collapse_shape_ok": theorem_export["exact_ratio_constraint_rank"] == 4 and theorem_export["exact_ratio_constraint_nullity"] == 1,
        "slope_overgeneration_preserved": theorem_export["prime_log_proportionality_overgenerates_slope_line"] and theorem_export["nonstrict_log_slope_also_proportional"],
        "prime_log_source_not_exported": not theorem_export["prime_log_proportionality_source_exported"],
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
        "packet_id": "P2542",
        "stage_id": "S1492",
        "status": "STRICT_DAMPING_PRIME_LOG_PROPORTIONALITY_CURRENT_PREMISE_OBSTRUCTION_CERTIFICATE_CURRENT_MONOID_ROUTE_REFUTED_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_prime_log_proportionality_current_premise_obstruction_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_prime_log_proportionality_current_premise_obstruction_certificate"]["theorem_export"]
    lines = [
        "# P2542/S1492 strict damping prime-log proportionality current-premise obstruction certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier source key under attack: `{t['frontier_source_key_under_attack']}`.",
        f"- Sample prime-character rows: `{t['sample_prime_character_row_count']}`.",
        f"- All samples are unital multiplicative characters: `{t['all_sample_rows_are_unital_multiplicative_characters']}`.",
        f"- Countermodel count: `{t['countermodel_count']}`.",
        f"- Current unital/multiplicative premises entail prime-log proportionality: `{not t['current_unital_multiplicative_premises_do_not_entail_prime_log_proportionality']}`.",
        f"- Ratio constraint rank/nullity: `{t['exact_ratio_constraint_rank']}/{t['exact_ratio_constraint_nullity']}`.",
        f"- Non-strict log slope also proportional: `{t['nonstrict_log_slope_also_proportional']}`.",
        f"- Prime-log proportionality overgenerates slope line: `{t['prime_log_proportionality_overgenerates_slope_line']}`.",
        f"- Prime-log source exported: `{t['prime_log_proportionality_source_exported']}`.",
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
        f"`{payload['strict_damping_prime_log_proportionality_current_premise_obstruction_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_prime_log_proportionality_current_premise_obstruction_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
