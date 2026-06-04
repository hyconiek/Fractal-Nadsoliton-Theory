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
OUT = GEN / "p2543_s1493_strict_damping_slope_value_current_premise_obstruction_certificate.json"
MD = GEN / "p2543_s1493_strict_damping_slope_value_current_premise_obstruction_certificate.md"

SOURCE_FILES = {
    "P2527_PRIME_LOG_SLOPE_LINE": GEN / "p2527_s1477_strict_damping_prime_log_proportionality_slope_line_certificate.json",
    "P2528_PRIME_SLOPE_ANCHOR_EQUIVALENCE": GEN / "p2528_s1478_strict_damping_prime_slope_anchor_equivalence_certificate.json",
    "P2530_FOUR_KEY_IRREDUNDANCY": GEN / "p2530_s1480_strict_damping_four_key_irredundancy_witness_certificate.json",
    "P2539_TOE_POTENTIAL_RECOMMENDATION": GEN / "p2539_s1489_strict_damping_toe_potential_recommendation_certificate.json",
    "P2542_PRIME_LOG_CURRENT_PREMISE_OBSTRUCTION": GEN / "p2542_s1492_strict_damping_prime_log_proportionality_current_premise_obstruction_certificate.json",
}

NODE_DOMAIN = list(range(1, 12))
PRIMES = [2, 3, 5, 7, 11]
SLOPE_CANDIDATES = [Fraction(-1, 1), Fraction(0, 1), Fraction(1, 2), Fraction(4, 5), Fraction(1, 1), Fraction(9, 5), Fraction(2, 1)]
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
        "new_packet": "P2543|S1493|slope value current-premise obstruction|prime-anchor source obstruction|slope source obstruction",
        "intended_research_nonduplication": "slope_value_or_prime_anchor_source|slope value source|prime slope anchor|delta=4/5|slope-line countermodel",
        "precursor_packets": "P2527|S1477|P2528|S1478|P2530|S1480|P2539|S1489|P2542|S1492",
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


def node_values_for_slope(slope: Fraction) -> list[float]:
    prime_values = [float(slope) * math.log(prime) for prime in PRIMES]
    return [
        sum(exponent * value for exponent, value in zip(factorization_exponents(d), prime_values))
        for d in NODE_DOMAIN
    ]


def max_multiplicative_residual(node_values: list[float]) -> float:
    return max(abs(node_values[de - 1] - node_values[d - 1] - node_values[e - 1]) for d, e, de in multiplicative_pairs())


def prime_values_for_slope(slope: Fraction) -> list[float]:
    return [float(slope) * math.log(prime) for prime in PRIMES]


def slope_line_rows() -> list[dict[str, Any]]:
    rows = []
    for slope in SLOPE_CANDIDATES:
        prime_values = prime_values_for_slope(slope)
        node_values = node_values_for_slope(slope)
        ratios = [value / math.log(prime) for prime, value in zip(PRIMES, prime_values)]
        ratio_spread = max(ratios) - min(ratios)
        residual = max_multiplicative_residual(node_values)
        rows.append({
            "slope_delta": frac_text(slope),
            "eta_if_slope_delta": frac_text(Fraction(1, 1) + slope),
            "is_strict_delta_4_over_5": slope == STRICT_DELTA,
            "prime_values_by_p_2_3_5_7_11": prime_values,
            "node_values_y_1_to_y_11": node_values,
            "unital_y1_zero_accepts": abs(node_values[0]) < TOL,
            "multiplicative_character_accepts": residual < TOL,
            "max_multiplicative_residual": residual,
            "prime_log_proportionality_accepts": ratio_spread < TOL,
            "normalized_ratio_spread": ratio_spread,
            "slope_value_or_prime_anchor_accepts": slope == STRICT_DELTA,
            "is_current_premise_countermodel_to_slope_value": slope != STRICT_DELTA,
        })
    return rows


def prime_anchor_residual_rows() -> list[dict[str, Any]]:
    rows = []
    for slope in SLOPE_CANDIDATES:
        for prime in PRIMES:
            candidate_value = float(slope) * math.log(prime)
            strict_value = float(STRICT_DELTA) * math.log(prime)
            rows.append({
                "slope_delta": frac_text(slope),
                "prime_anchor": prime,
                "candidate_value": candidate_value,
                "strict_anchor_value": strict_value,
                "anchor_residual": candidate_value - strict_value,
                "anchor_accepts": abs(candidate_value - strict_value) < TOL,
            })
    return rows


def build_obstruction_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2527 = theorem(sources["P2527_PRIME_LOG_SLOPE_LINE"], "strict_damping_prime_log_proportionality_slope_line_certificate")
    p2528 = theorem(sources["P2528_PRIME_SLOPE_ANCHOR_EQUIVALENCE"], "strict_damping_prime_slope_anchor_equivalence_certificate")
    p2530 = theorem(sources["P2530_FOUR_KEY_IRREDUNDANCY"], "strict_damping_four_key_irredundancy_witness_certificate")
    p2539 = theorem(sources["P2539_TOE_POTENTIAL_RECOMMENDATION"], "strict_damping_toe_potential_recommendation_certificate")
    p2542 = theorem(sources["P2542_PRIME_LOG_CURRENT_PREMISE_OBSTRUCTION"], "strict_damping_prime_log_proportionality_current_premise_obstruction_certificate")
    rows = slope_line_rows()
    countermodels = [row for row in rows if row["is_current_premise_countermodel_to_slope_value"]]
    strict_rows = [row for row in rows if row["is_strict_delta_4_over_5"]]
    anchor_rows = prime_anchor_residual_rows()
    chosen_countermodel = next(row for row in countermodels if row["slope_delta"] == "1/2")
    return {
        "frontier_source_key_under_attack": "slope_value_or_prime_anchor_source",
        "p2527_slope_line_inherited": p2527.get("prime_log_proportionality_subkey_exported") is True,
        "p2528_prime_anchor_equivalence_inherited": p2528.get("conditional_prime_slope_anchor_equivalence_exported") is True,
        "p2530_four_key_irredundancy_inherited": p2530.get("four_key_irredundancy_witness_exported") is True,
        "p2539_next_step_recommendation_inherited": p2539.get("recommended_next_honest_step") == "prove_or_refute_one_strict_source_theorem_for_a_single_P2529_source_key_before_more_bookkeeping_layers",
        "p2542_prime_log_route_obstruction_inherited": p2542.get("current_premise_obstruction_exported") is True,
        "slope_line_rows": rows,
        "prime_anchor_residual_rows": anchor_rows,
        "current_slope_line_premises_do_not_entail_delta_4_over_5": len(countermodels) > 0,
        "countermodel_count": len(countermodels),
        "chosen_countermodel_delta_1_over_2": chosen_countermodel,
        "strict_delta_row_count": len(strict_rows),
        "all_slope_candidates_pass_unital_multiplicative_prime_log_premises": all(
            row["unital_y1_zero_accepts"] and row["multiplicative_character_accepts"] and row["prime_log_proportionality_accepts"]
            for row in rows
        ),
        "only_strict_delta_passes_slope_value_or_prime_anchor_filter": (
            len(strict_rows) == 1
            and all(row["slope_value_or_prime_anchor_accepts"] == row["is_strict_delta_4_over_5"] for row in rows)
        ),
        "prime_anchor_equivalence_available_but_unsourced": p2528.get("conditional_prime_slope_anchor_equivalence_exported") is True and not p2528.get("slope_value_source_exported", False),
        "current_premise_obstruction_exported": True,
        "required_new_premise_class": "strict dynamics must source the numeric slope delta=4/5 or an equivalent prime-value anchor v_p=(4/5)log(p)",
    }


def append_doc_sections() -> None:
    eq_section = """
## P2543/S1493 strict damping slope-value current-premise obstruction certificate

`P2543/S1493` attacks the final numeric source key `slope_value_or_prime_anchor_source`.  It inherits P2527 prime-log proportionality slope-line reduction, P2528 conditional prime-anchor equivalence, P2530 four-key irredundancy, P2539 source-key recommendation, and P2542 prime-log current-premise obstruction.  The new audit enumerates the same slope line `v_p=a log p` for candidate slopes including `a=1/2` and `a=4/5`.  Every audited slope is unital, multiplicative, and prime-log proportional, but only `a=4/5` satisfies the slope-value/prime-anchor target.

Therefore the current slope-line premises do not entail the strict slope `delta=4/5`; `delta=1/2` is an explicit current-premise countermodel.  P2528 shows that any single prime anchor would conditionally select the strict slope, but the anchor value/placement remains unsourced.  No strict source key, source-obligation discharge, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, physical-value generator, or ToE closure is exported.
"""
    lag_section = """
## P2543/S1493 slope-value source obstruction guard

`P2543/S1493` shows that the prime-log slope line still does not source `delta=4/5`: slopes such as `delta=1/2` satisfy unital multiplicativity and prime-log proportionality.  A strict numeric slope or prime-value anchor remains an external source obligation, so the damping term is not yet role-bearing in `L_total`.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2543/S1493 strict damping slope-value current-premise obstruction certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2543/S1493 slope-value source obstruction guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    cert = build_obstruction_certificate(sources)
    theorem_export = {
        "theorem_name": "P2543_T1_strict_damping_slope_value_current_premise_obstruction_certificate",
        "audited_chain": ["P2527/S1477", "P2528/S1478", "P2530/S1480", "P2539/S1489", "P2542/S1492"],
        "strict_damping_slope_value_current_premise_obstruction_certificate": cert,
        "frontier_source_key_under_attack": cert["frontier_source_key_under_attack"],
        "p2527_slope_line_inherited": cert["p2527_slope_line_inherited"],
        "p2528_prime_anchor_equivalence_inherited": cert["p2528_prime_anchor_equivalence_inherited"],
        "p2530_four_key_irredundancy_inherited": cert["p2530_four_key_irredundancy_inherited"],
        "p2539_next_step_recommendation_inherited": cert["p2539_next_step_recommendation_inherited"],
        "p2542_prime_log_route_obstruction_inherited": cert["p2542_prime_log_route_obstruction_inherited"],
        "slope_line_row_count": len(cert["slope_line_rows"]),
        "all_slope_candidates_pass_unital_multiplicative_prime_log_premises": cert["all_slope_candidates_pass_unital_multiplicative_prime_log_premises"],
        "countermodel_count": cert["countermodel_count"],
        "chosen_countermodel_delta_1_over_2": cert["chosen_countermodel_delta_1_over_2"],
        "strict_delta_row_count": cert["strict_delta_row_count"],
        "current_slope_line_premises_do_not_entail_delta_4_over_5": cert["current_slope_line_premises_do_not_entail_delta_4_over_5"],
        "only_strict_delta_passes_slope_value_or_prime_anchor_filter": cert["only_strict_delta_passes_slope_value_or_prime_anchor_filter"],
        "prime_anchor_equivalence_available_but_unsourced": cert["prime_anchor_equivalence_available_but_unsourced"],
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
            "P2543 refutes only the current prime-log slope-line route to slope_value_or_prime_anchor_source.",
            "The strict slope delta=4/5 remains a target value, not a source theorem.",
            "P2528 prime anchors are conditionally sufficient but remain unsourced.",
            "No bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Stop treating strict damping beta/eta as sourced until a new nadsoliton theorem supplies one of the four blocked source keys; prioritize a new source principle rather than more repair bookkeeping.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "all_required_precursors_inherited": all(theorem_export[key] for key in [
            "p2527_slope_line_inherited",
            "p2528_prime_anchor_equivalence_inherited",
            "p2530_four_key_irredundancy_inherited",
            "p2539_next_step_recommendation_inherited",
            "p2542_prime_log_route_obstruction_inherited",
        ]),
        "slope_line_countermodel_exists": theorem_export["countermodel_count"] > 0 and theorem_export["current_slope_line_premises_do_not_entail_delta_4_over_5"],
        "strict_filter_shape_ok": theorem_export["strict_delta_row_count"] == 1 and theorem_export["only_strict_delta_passes_slope_value_or_prime_anchor_filter"],
        "all_slope_rows_pass_current_premises": theorem_export["all_slope_candidates_pass_unital_multiplicative_prime_log_premises"],
        "anchor_equivalence_not_source": theorem_export["prime_anchor_equivalence_available_but_unsourced"] and not theorem_export["slope_value_or_prime_anchor_source_exported"],
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
        "packet_id": "P2543",
        "stage_id": "S1493",
        "status": "STRICT_DAMPING_SLOPE_VALUE_CURRENT_PREMISE_OBSTRUCTION_CERTIFICATE_CURRENT_SLOPE_LINE_ROUTE_REFUTED_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_slope_value_current_premise_obstruction_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_slope_value_current_premise_obstruction_certificate"]["theorem_export"]
    lines = [
        "# P2543/S1493 strict damping slope-value current-premise obstruction certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier source key under attack: `{t['frontier_source_key_under_attack']}`.",
        f"- Slope-line row count: `{t['slope_line_row_count']}`.",
        f"- All slope candidates pass current unital/multiplicative/prime-log premises: `{t['all_slope_candidates_pass_unital_multiplicative_prime_log_premises']}`.",
        f"- Countermodel count: `{t['countermodel_count']}`.",
        f"- Current slope-line premises entail delta=4/5: `{not t['current_slope_line_premises_do_not_entail_delta_4_over_5']}`.",
        f"- Only strict delta passes slope/anchor filter: `{t['only_strict_delta_passes_slope_value_or_prime_anchor_filter']}`.",
        f"- Prime-anchor equivalence available but unsourced: `{t['prime_anchor_equivalence_available_but_unsourced']}`.",
        f"- Slope source exported: `{t['slope_value_or_prime_anchor_source_exported']}`.",
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
        f"`{payload['strict_damping_slope_value_current_premise_obstruction_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_slope_value_current_premise_obstruction_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
