#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from fractions import Fraction
from itertools import product
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
OUT = GEN / "p2526_s1476_strict_damping_finite_monoid_prime_character_nullity_certificate.json"
MD = GEN / "p2526_s1476_strict_damping_finite_monoid_prime_character_nullity_certificate.md"

SOURCE_FILES = {
    "P2525_MULTIPLICATIVE_BETA_NORMALIZATION": GEN / "p2525_s1475_strict_damping_multiplicative_cocycle_beta_normalization_subkey_certificate.json",
}

NODE_DOMAIN = list(range(1, 12))
PRIMES = [2, 3, 5, 7, 11]
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
        "new_packet": "P2526|S1476|finite monoid prime character nullity|prime-character nullity|multiplicative character rank|nullity certificate|prime generator freedom",
        "intended_research_nonduplication": "prime generator|prime-value|finite monoid|monoid character|factorization matrix|prime basis|multiplicative.*rank|nullity",
        "precursor_packets": "P2525|S1475|multiplicative-cocycle beta-normalization|P2524|affine consistency continuum|P2523|pairwise secant consistency",
        "source_obligation_language": "prime proportionality source|slope value source|beta_eta_numeric_source|m2_operator_signature_source|proper subset obstruction",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def multiplicative_pairs() -> list[tuple[int, int, int]]:
    return [(d, e, d * e) for d, e in product(NODE_DOMAIN, repeat=2) if d * e in NODE_DOMAIN]


def constraint_matrix() -> list[list[Fraction]]:
    matrix = []
    for d, e, de in multiplicative_pairs():
        row = [Fraction(0, 1) for _ in NODE_DOMAIN]
        row[de - 1] += Fraction(1, 1)
        row[d - 1] -= Fraction(1, 1)
        row[e - 1] -= Fraction(1, 1)
        matrix.append(row)
    return matrix


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


def mat_vec(matrix: list[list[int]], vector: list[float]) -> list[float]:
    return [sum(coef * value for coef, value in zip(row, vector)) for row in matrix]


def max_constraint_residual(values: list[float]) -> float:
    return max(abs(values[de - 1] - values[d - 1] - values[e - 1]) for d, e, de in multiplicative_pairs())


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


def sample_prime_character_rows() -> list[dict[str, Any]]:
    samples = {
        "strict_prime_log_character": [float(STRICT_DELTA) * math.log(p) for p in PRIMES],
        "unit_p2_only_character": [1.0, 0.0, 0.0, 0.0, 0.0],
        "unit_p3_only_character": [0.0, 1.0, 0.0, 0.0, 0.0],
        "perturbed_strict_p2_character": [float(STRICT_DELTA) * math.log(2) + 0.25] + [float(STRICT_DELTA) * math.log(p) for p in PRIMES[1:]],
    }
    rows = []
    matrix = factorization_matrix()
    for name, prime_values in samples.items():
        node_values = mat_vec(matrix, prime_values)
        stats = secant_stats(node_values)
        rows.append({
            "sample_name": name,
            "prime_values_by_p_2_3_5_7_11": prime_values,
            "node_values_y_1_to_y_11": node_values,
            "max_multiplicative_constraint_residual": max_constraint_residual(node_values),
            "multiplicative_constraints_accept": max_constraint_residual(node_values) < TOL,
            "affine_secants_constant": stats["affine_secants_constant"],
            "secant_spread": stats["secant_spread"],
            "is_strict_numeric_target": name == "strict_prime_log_character",
        })
    return rows


def prime_ratio_rows() -> list[dict[str, Any]]:
    return [{
        "prime": p,
        "strict_prime_value": float(STRICT_DELTA) * math.log(p),
        "strict_value_over_log_prime": float(STRICT_DELTA),
    } for p in PRIMES]


def append_doc_sections() -> None:
    eq_section = """
`P2526/S1476` strengthens the P2525 guard by removing the affine assumption and auditing raw node variables `y_1,...,y_11` under the finite multiplicative law `y_{de}=y_d+y_e` for all audited products `de<=11`.  The exact rational constraint matrix has rank `6` and nullity `5`: the solution space is the finite monoid-character family freely determined by the prime-generator values at `2,3,5,7,11`.  Thus multiplicativity alone does not derive the strict log slope; affine/log-proportionality across prime generators and then a slope-value source are still independent obligations before `beta_eta_numeric_source` can be claimed.
""".strip()
    lag_section = """
`P2526/S1476` records that a multiplicative node law, even if sourced, leaves a five-dimensional prime-character freedom on the audited `d=1..11` monoid.  The strict damping term still needs a source for prime log-proportionality and for the numeric slope `delta=4/5`, in addition to the `m=2` operator signature; no nonlinear compression-flow source or role-bearing `L_total` term is licensed.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "P2526/S1476", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "P2526/S1476", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2525 = theorem(sources["P2525_MULTIPLICATIVE_BETA_NORMALIZATION"], "strict_damping_multiplicative_cocycle_beta_normalization_subkey_certificate")
    constraints = constraint_matrix()
    reduced, pivots = rref(constraints)
    rank = len(pivots)
    nullity = len(NODE_DOMAIN) - rank
    free_columns = [col for col in range(len(NODE_DOMAIN)) if col not in pivots]
    samples = sample_prime_character_rows()
    theorem_export = {
        "frontier_atom_under_attack": "multiplicative_character_law_strength_for_strict_damping_numeric_source",
        "p2525_beta_normalization_subkey_inherited": bool(p2525.get("conditional_beta_normalization_subkey_exported", False)),
        "node_domain": NODE_DOMAIN,
        "prime_generators": PRIMES,
        "multiplicative_pair_count_on_domain_1_to_11": len(multiplicative_pairs()),
        "constraint_row_count": len(constraints),
        "constraint_column_count": len(NODE_DOMAIN),
        "exact_constraint_rank": rank,
        "exact_constraint_nullity": nullity,
        "rref_pivot_variables_y_indices_1_based": [col + 1 for col in pivots],
        "rref_free_variables_y_indices_1_based": [col + 1 for col in free_columns],
        "canonical_prime_parameter_variables_y_indices_1_based": PRIMES,
        "factorization_matrix_rows_y_1_to_y_11_columns_p_2_3_5_7_11": factorization_matrix(),
        "finite_monoid_character_family_dimension": len(PRIMES),
        "factorization_parameterization_dimension_matches_nullity": nullity == len(PRIMES),
        "multiplicative_law_alone_leaves_prime_generator_freedom": nullity == len(PRIMES),
        "strict_log_character_is_one_member_of_prime_character_family": any(row["is_strict_numeric_target"] and row["multiplicative_constraints_accept"] for row in samples),
        "non_affine_prime_characters_pass_multiplicativity": any(row["multiplicative_constraints_accept"] and not row["affine_secants_constant"] for row in samples),
        "prime_log_proportionality_needed_to_recover_affine_slope_line": True,
        "slope_value_source_still_needed_for_delta_4_over_5": True,
        "conditional_finite_monoid_character_nullity_exported": True,
        "multiplicative_character_law_source_exported": False,
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
        "strict_damping_finite_monoid_prime_character_nullity_certificate": {
            "multiplicative_pairs": [{"d": d, "e": e, "de": de} for d, e, de in multiplicative_pairs()],
            "constraint_matrix_rows": [[str(value) for value in row] for row in constraints],
            "rref_rows": [[str(value) for value in row] for row in reduced],
            "sample_prime_character_rows": samples,
            "prime_ratio_rows_for_strict_log_character": prime_ratio_rows(),
            "source_obligation_normal_form": "beta_eta_numeric_source = multiplicative_character_law_source AND prime_log_proportionality_source AND slope_value_source; strict_damping_beta_eta_source additionally requires m2_operator_signature_source",
        },
    }
    gatekeepers = {
        "p2525_inherited": theorem_export["p2525_beta_normalization_subkey_inherited"],
        "rank_nullity_exact_prime_dimension": rank == 6 and nullity == 5 and theorem_export["factorization_parameterization_dimension_matches_nullity"],
        "multiplicativity_not_numeric_slope_source": theorem_export["non_affine_prime_characters_pass_multiplicativity"] and theorem_export["slope_value_source_still_needed_for_delta_4_over_5"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "multiplicative_character_law_source_exported",
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
        "packet_id": "P2526",
        "stage_id": "S1476",
        "status": "STRICT_DAMPING_FINITE_MONOID_PRIME_CHARACTER_NULLITY_CERTIFICATE_CONDITIONAL_CHARACTER_FAMILY_ONLY_NO_PRIME_PROPORTIONALITY_SOURCE_NO_SLOPE_SOURCE_NO_OPERATOR_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_finite_monoid_prime_character_nullity_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_finite_monoid_prime_character_nullity_certificate"]["theorem_export"]
    lines = [
        "# P2526/S1476 strict damping finite-monoid prime-character nullity certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2525 beta-normalization subkey inherited: `{t['p2525_beta_normalization_subkey_inherited']}`.",
        f"- Multiplicative product pairs on `d=1..11`: `{t['multiplicative_pair_count_on_domain_1_to_11']}`.",
        f"- Exact constraint rank/nullity: `{t['exact_constraint_rank']}/{t['exact_constraint_nullity']}`.",
        f"- Canonical prime parameter variables: `{t['canonical_prime_parameter_variables_y_indices_1_based']}`.",
        f"- Factorization parameterization dimension matches nullity: `{t['factorization_parameterization_dimension_matches_nullity']}`.",
        f"- Multiplicative law alone leaves prime-generator freedom: `{t['multiplicative_law_alone_leaves_prime_generator_freedom']}`.",
        f"- Non-affine prime characters pass multiplicativity: `{t['non_affine_prime_characters_pass_multiplicativity']}`.",
        f"- Prime log-proportionality source exported: `{t['prime_log_proportionality_source_exported']}`.",
        f"- Beta/eta numeric source exported: `{t['beta_eta_numeric_source_exported']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports only a conditional finite-monoid character-family/nullity theorem. It does not source the multiplicative law, prime log-proportionality, slope value, m=2 operator signature, bridge completion, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_finite_monoid_prime_character_nullity_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_finite_monoid_prime_character_nullity_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
