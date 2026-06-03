#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from fractions import Fraction
from statistics import mean, pstdev
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
OUT = GEN / "p2527_s1477_strict_damping_prime_log_proportionality_slope_line_certificate.json"
MD = GEN / "p2527_s1477_strict_damping_prime_log_proportionality_slope_line_certificate.md"

SOURCE_FILES = {
    "P2526_FINITE_MONOID_PRIME_CHARACTER_NULLITY": GEN / "p2526_s1476_strict_damping_finite_monoid_prime_character_nullity_certificate.json",
}

NODE_DOMAIN = list(range(1, 12))
PRIMES = [2, 3, 5, 7, 11]
SLOPE_CANDIDATES = [Fraction(-1, 1), Fraction(0, 1), Fraction(1, 2), Fraction(4, 5), Fraction(1, 1), Fraction(9, 5), Fraction(2, 1)]
STRICT_DELTA = Fraction(4, 5)
STRICT_ETA = Fraction(9, 5)
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
    return {"count": len(lines), "samples": lines[:40]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2527|S1477|prime log proportionality slope line|prime-log proportionality|log-proportional prime character|slope-line certificate",
        "intended_research_nonduplication": "prime log proportionality|log proportionality|prime proportionality|ratio.*log prime|prime ratio|prime-generator proportionality",
        "precursor_packets": "P2526|S1476|finite monoid prime character nullity|P2525|multiplicative-cocycle beta-normalization|P2524|affine consistency continuum",
        "source_obligation_language": "slope value source|prime proportionality source|beta_eta_numeric_source|m2_operator_signature_source|proper subset obstruction",
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
    for p in PRIMES:
        exponent = 0
        while value % p == 0:
            value //= p
            exponent += 1
        exponents.append(exponent)
    if value != 1:
        raise ValueError(f"unexpected factor outside audited primes: {n}")
    return exponents


def factorization_matrix() -> list[list[int]]:
    return [factorization_exponents(d) for d in NODE_DOMAIN]


def exact_ratio_constraint_matrix() -> list[list[Fraction]]:
    # Variables are normalized prime ratios r_p = v_p / log(p).  The constraints
    # r_p - r_2 = 0 for p in {3,5,7,11} have exact rank 4 and leave one slope.
    rows = []
    for index in range(1, len(PRIMES)):
        row = [Fraction(0, 1) for _ in PRIMES]
        row[0] = Fraction(-1, 1)
        row[index] = Fraction(1, 1)
        rows.append(row)
    return rows


def rref(matrix: list[list[Fraction]]) -> tuple[list[list[Fraction]], list[int]]:
    work = [row[:] for row in matrix]
    if not work:
        return work, []
    row_count = len(work)
    col_count = len(work[0])
    pivots: list[int] = []
    pivot_row = 0
    for col in range(col_count):
        pivot = next((r for r in range(pivot_row, row_count) if work[r][col] != 0), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][col]
        work[pivot_row] = [value / scale for value in work[pivot_row]]
        for r in range(row_count):
            if r == pivot_row or work[r][col] == 0:
                continue
            factor = work[r][col]
            work[r] = [a - factor * b for a, b in zip(work[r], work[pivot_row])]
        pivots.append(col)
        pivot_row += 1
        if pivot_row == row_count:
            break
    return work, pivots


def prime_values_for_slope(slope: Fraction) -> list[float]:
    return [float(slope) * math.log(p) for p in PRIMES]


def node_values_from_prime_values(prime_values: list[float]) -> list[float]:
    return [sum(exponent * value for exponent, value in zip(row, prime_values)) for row in factorization_matrix()]


def least_squares_slope(prime_values: list[float]) -> dict[str, Any]:
    logs = [math.log(p) for p in PRIMES]
    normal = sum(value * value for value in logs)
    slope = sum(log_p * v_p for log_p, v_p in zip(logs, prime_values)) / normal
    residuals = [v_p - slope * log_p for log_p, v_p in zip(logs, prime_values)]
    return {
        "normal_scalar_sum_log_prime_squared": normal,
        "rank_one_design_positive_normal": normal > 0,
        "least_squares_slope": slope,
        "residuals_by_prime_2_3_5_7_11": residuals,
        "max_abs_residual": max(abs(value) for value in residuals),
        "prime_log_proportionality_accepts": max(abs(value) for value in residuals) < TOL,
    }


def secant_stats(values: list[float]) -> dict[str, Any]:
    secants = []
    for i, d_i in enumerate(NODE_DOMAIN):
        for j in range(i + 1, len(NODE_DOMAIN)):
            d_j = NODE_DOMAIN[j]
            secants.append((values[j] - values[i]) / (math.log(d_j) - math.log(d_i)))
    return {
        "pair_count": len(secants),
        "secant_min": min(secants),
        "secant_max": max(secants),
        "secant_mean": mean(secants),
        "secant_population_std": pstdev(secants),
        "secant_spread": max(secants) - min(secants),
        "affine_secants_constant": max(secants) - min(secants) < TOL,
    }


def slope_line_rows() -> list[dict[str, Any]]:
    rows = []
    for slope in SLOPE_CANDIDATES:
        prime_values = prime_values_for_slope(slope)
        node_values = node_values_from_prime_values(prime_values)
        lsq = least_squares_slope(prime_values)
        stats = secant_stats(node_values)
        rows.append({
            "slope_delta": frac_text(slope),
            "eta_if_slope_delta": frac_text(Fraction(1, 1) + slope),
            "prime_values_by_p_2_3_5_7_11": prime_values,
            "node_values_y_1_to_y_11": node_values,
            "least_squares_audit": lsq,
            "secant_stats": stats,
            "is_strict_numeric_target": slope == STRICT_DELTA,
            "prime_log_proportionality_accepts": lsq["prime_log_proportionality_accepts"],
            "affine_node_line_accepts": stats["affine_secants_constant"],
        })
    return rows


def negative_sample_rows() -> list[dict[str, Any]]:
    strict = prime_values_for_slope(STRICT_DELTA)
    samples = {
        "unit_p2_only_character": [1.0, 0.0, 0.0, 0.0, 0.0],
        "unit_p3_only_character": [0.0, 1.0, 0.0, 0.0, 0.0],
        "perturbed_strict_p2_character": [strict[0] + 0.25] + strict[1:],
    }
    rows = []
    for name, prime_values in samples.items():
        node_values = node_values_from_prime_values(prime_values)
        lsq = least_squares_slope(prime_values)
        stats = secant_stats(node_values)
        rows.append({
            "sample_name": name,
            "prime_values_by_p_2_3_5_7_11": prime_values,
            "least_squares_audit": lsq,
            "secant_stats": stats,
            "prime_log_proportionality_rejects": not lsq["prime_log_proportionality_accepts"],
            "affine_node_line_rejects": not stats["affine_secants_constant"],
        })
    return rows


def append_doc_sections() -> None:
    eq_section = """
`P2527/S1477` audits the next subkey after P2526.  If the five prime-generator values are constrained to be log-proportional, `v_p = a log p`, the exact normalized-ratio constraint matrix has rank `4` and nullity `1`; combined with multiplicativity this collapses the finite prime-character family to the one-parameter affine node line `y_d=a log d`.  The finite audit accepts all tested slopes on that line and rejects representative non-proportional prime characters, so this is still not a numeric slope theorem: `a=4/5` remains an independent slope-value source obligation before `beta_eta_numeric_source` can be claimed.
""".strip()
    lag_section = """
`P2527/S1477` shows that a sourced prime log-proportionality law would reduce the P2526 five-dimensional prime-character freedom to a one-dimensional slope line, but it would not choose the strict slope `delta=4/5`.  The damping term still needs a numeric slope source and the `m=2` operator-signature source; no nonlinear compression-flow source or role-bearing `L_total` term is licensed.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "P2527/S1477", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "P2527/S1477", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2526 = theorem(sources["P2526_FINITE_MONOID_PRIME_CHARACTER_NULLITY"], "strict_damping_finite_monoid_prime_character_nullity_certificate")
    constraints = exact_ratio_constraint_matrix()
    reduced, pivots = rref(constraints)
    rank = len(pivots)
    nullity = len(PRIMES) - rank
    line_rows = slope_line_rows()
    rejected_rows = negative_sample_rows()
    strict_rows = [row for row in line_rows if row["is_strict_numeric_target"]]
    theorem_export = {
        "frontier_atom_under_attack": "prime_log_proportionality_subkey_between_multiplicative_character_family_and_numeric_slope_source",
        "p2526_prime_character_nullity_inherited": bool(p2526.get("conditional_finite_monoid_character_nullity_exported", False)),
        "prime_generators": PRIMES,
        "slope_candidate_count": len(line_rows),
        "exact_normalized_ratio_constraint_rank": rank,
        "exact_normalized_ratio_constraint_nullity": nullity,
        "one_dimensional_slope_line_after_prime_log_proportionality": rank == 4 and nullity == 1,
        "all_slope_candidates_pass_prime_log_proportionality": all(row["prime_log_proportionality_accepts"] for row in line_rows),
        "all_slope_candidates_give_affine_node_line": all(row["affine_node_line_accepts"] for row in line_rows),
        "strict_numeric_target_is_one_slope_line_member": len(strict_rows) == 1 and strict_rows[0]["prime_log_proportionality_accepts"],
        "representative_nonproportional_prime_characters_rejected": all(row["prime_log_proportionality_rejects"] and row["affine_node_line_rejects"] for row in rejected_rows),
        "prime_log_proportionality_subkey_exported": True,
        "prime_log_proportionality_source_exported": False,
        "slope_value_source_exported": False,
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
        "strict_damping_prime_log_proportionality_slope_line_certificate": {
            "exact_normalized_ratio_constraint_matrix_rows": [[str(value) for value in row] for row in constraints],
            "exact_normalized_ratio_rref_rows": [[str(value) for value in row] for row in reduced],
            "slope_line_rows": line_rows,
            "representative_nonproportional_rejection_rows": rejected_rows,
            "source_obligation_normal_form": "beta_eta_numeric_source = multiplicative_character_law_source AND prime_log_proportionality_source AND slope_value_source; strict_damping_beta_eta_source additionally requires m2_operator_signature_source",
        },
    }
    gatekeepers = {
        "p2526_inherited": theorem_export["p2526_prime_character_nullity_inherited"],
        "prime_log_proportionality_collapses_to_slope_line_only": theorem_export["one_dimensional_slope_line_after_prime_log_proportionality"] and theorem_export["all_slope_candidates_pass_prime_log_proportionality"],
        "strict_target_not_unique_without_slope_source": theorem_export["strict_numeric_target_is_one_slope_line_member"] and theorem_export["slope_candidate_count"] > 1,
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "prime_log_proportionality_source_exported",
            "slope_value_source_exported",
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
        "packet_id": "P2527",
        "stage_id": "S1477",
        "status": "STRICT_DAMPING_PRIME_LOG_PROPORTIONALITY_SLOPE_LINE_CERTIFICATE_CONDITIONAL_SUBKEY_ONLY_NO_SLOPE_VALUE_SOURCE_NO_OPERATOR_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_prime_log_proportionality_slope_line_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_prime_log_proportionality_slope_line_certificate"]["theorem_export"]
    lines = [
        "# P2527/S1477 strict damping prime-log-proportionality slope-line certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2526 prime-character nullity inherited: `{t['p2526_prime_character_nullity_inherited']}`.",
        f"- Exact normalized-ratio rank/nullity: `{t['exact_normalized_ratio_constraint_rank']}/{t['exact_normalized_ratio_constraint_nullity']}`.",
        f"- One-dimensional slope line after prime-log proportionality: `{t['one_dimensional_slope_line_after_prime_log_proportionality']}`.",
        f"- Slope candidates on line: `{t['slope_candidate_count']}`.",
        f"- Strict numeric target is one slope-line member: `{t['strict_numeric_target_is_one_slope_line_member']}`.",
        f"- Non-proportional prime characters rejected: `{t['representative_nonproportional_prime_characters_rejected']}`.",
        f"- Prime-log proportionality source exported: `{t['prime_log_proportionality_source_exported']}`.",
        f"- Slope value source exported: `{t['slope_value_source_exported']}`.",
        f"- Beta/eta numeric source exported: `{t['beta_eta_numeric_source_exported']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports only a conditional prime-log-proportionality subkey that collapses finite prime-character freedom to a one-dimensional slope line. It does not source that law, does not choose the numeric slope, does not source the m=2 operator signature, and exports no bridge completion, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_prime_log_proportionality_slope_line_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_prime_log_proportionality_slope_line_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
