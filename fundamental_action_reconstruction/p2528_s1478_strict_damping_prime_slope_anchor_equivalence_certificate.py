#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from fractions import Fraction
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
OUT = GEN / "p2528_s1478_strict_damping_prime_slope_anchor_equivalence_certificate.json"
MD = GEN / "p2528_s1478_strict_damping_prime_slope_anchor_equivalence_certificate.md"

SOURCE_FILES = {
    "P2527_PRIME_LOG_PROPORTIONALITY_SLOPE_LINE": GEN / "p2527_s1477_strict_damping_prime_log_proportionality_slope_line_certificate.json",
}

PRIMES = [2, 3, 5, 7, 11]
NODE_DOMAIN = list(range(1, 12))
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
        "new_packet": "P2528|S1478|prime slope anchor equivalence|prime-value slope anchor|single-prime slope anchor|prime-generator anchor equivalence",
        "intended_research_nonduplication": "prime value anchor|prime slope anchor|prime-generator anchor|generator value anchor|single-prime|one-prime|prime.*anchor",
        "precursor_packets": "P2527|S1477|prime log proportionality slope line|P2526|finite monoid prime character nullity|P2521|single-node anchor equivalence",
        "source_obligation_language": "anchor placement source|prime value source|slope value source|beta_eta_numeric_source|proper subset obstruction",
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


def node_values_for_slope(slope: Fraction) -> list[float]:
    prime_values = [float(slope) * math.log(p) for p in PRIMES]
    return [sum(exponent * value for exponent, value in zip(factorization_exponents(d), prime_values)) for d in NODE_DOMAIN]


def secant_spread(values: list[float]) -> float:
    secants = []
    for i, d_i in enumerate(NODE_DOMAIN):
        for j in range(i + 1, len(NODE_DOMAIN)):
            d_j = NODE_DOMAIN[j]
            secants.append((values[j] - values[i]) / (math.log(d_j) - math.log(d_i)))
    return max(secants) - min(secants)


def prime_anchor_rows() -> list[dict[str, Any]]:
    rows = []
    for prime in PRIMES:
        log_p = math.log(prime)
        target_value = float(STRICT_DELTA) * log_p
        solved_slope = target_value / log_p
        node_values = node_values_for_slope(STRICT_DELTA)
        rows.append({
            "prime_anchor": prime,
            "log_prime": log_p,
            "strict_prime_value_anchor": target_value,
            "determinant_log_prime_positive": log_p > 0,
            "solved_slope_delta": solved_slope,
            "solved_eta": 1.0 + solved_slope,
            "recovers_strict_delta": abs(solved_slope - float(STRICT_DELTA)) < TOL,
            "recovers_strict_eta": abs(1.0 + solved_slope - float(STRICT_ETA)) < TOL,
            "reconstructed_node_values_y_1_to_y_11": node_values,
            "reconstructed_node_secant_spread": secant_spread(node_values),
            "reconstructs_affine_strict_node_line": secant_spread(node_values) < TOL,
        })
    return rows


def candidate_grid_audit() -> list[dict[str, Any]]:
    rows = []
    for prime in PRIMES:
        log_p = math.log(prime)
        strict_value = float(STRICT_DELTA) * log_p
        accepted = []
        candidates = []
        for slope in SLOPE_CANDIDATES:
            candidate_value = float(slope) * log_p
            residual = candidate_value - strict_value
            accepted_flag = abs(residual) < TOL
            if accepted_flag:
                accepted.append(frac_text(slope))
            candidates.append({
                "slope_delta_candidate": frac_text(slope),
                "eta_if_slope_delta": frac_text(Fraction(1, 1) + slope),
                "candidate_prime_value": candidate_value,
                "anchor_residual": residual,
                "accepted_by_prime_anchor": accepted_flag,
            })
        rows.append({
            "prime_anchor": prime,
            "accepted_slope_candidates": accepted,
            "only_strict_slope_accepted_on_grid": accepted == ["4/5"],
            "candidate_rows": candidates,
        })
    return rows


def append_doc_sections() -> None:
    eq_section = """
`P2528/S1478` audits the remaining slope-value subkey after P2527.  Once prime log-proportionality has reduced the P2526 character family to `v_p=a log p`, any single prime-generator value anchor `v_p=(4/5)log p` with `p in {2,3,5,7,11}` has positive determinant `log p` and conditionally recovers `a=4/5`, hence `eta=9/5`.  The all-prime audit and finite slope-candidate grid show equivalence of the five single-prime anchors, not a source theorem: the strict dynamics still must source which prime/value anchor, or otherwise source the numeric slope directly.
""".strip()
    lag_section = """
`P2528/S1478` shows that single prime-value anchors are sufficient but placement-nonunique slope keys on the P2527 slope line.  Since the prime anchor placement/value remains unsourced, this still does not export `beta_eta_numeric_source`, the `m=2` operator-signature source, nonlinear compression-flow source, or a role-bearing `L_total` term.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "P2528/S1478", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "P2528/S1478", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2527 = theorem(sources["P2527_PRIME_LOG_PROPORTIONALITY_SLOPE_LINE"], "strict_damping_prime_log_proportionality_slope_line_certificate")
    anchor_rows = prime_anchor_rows()
    grid_rows = candidate_grid_audit()
    theorem_export = {
        "frontier_atom_under_attack": "slope_value_source_subkey_after_prime_log_proportionality_slope_line",
        "p2527_slope_line_inherited": bool(p2527.get("prime_log_proportionality_subkey_exported", False)),
        "prime_anchor_count": len(anchor_rows),
        "slope_candidate_count_per_anchor": len(SLOPE_CANDIDATES),
        "all_prime_anchors_have_positive_log_determinant": all(row["determinant_log_prime_positive"] for row in anchor_rows),
        "all_single_prime_anchors_recover_strict_delta": all(row["recovers_strict_delta"] for row in anchor_rows),
        "all_single_prime_anchors_recover_strict_eta": all(row["recovers_strict_eta"] for row in anchor_rows),
        "all_single_prime_anchors_reconstruct_strict_node_line": all(row["reconstructs_affine_strict_node_line"] for row in anchor_rows),
        "all_candidate_grid_audits_accept_only_strict_slope": all(row["only_strict_slope_accepted_on_grid"] for row in grid_rows),
        "single_prime_anchor_placement_nonunique": len(anchor_rows) == len(PRIMES) and all(row["recovers_strict_delta"] for row in anchor_rows),
        "conditional_prime_slope_anchor_equivalence_exported": True,
        "prime_anchor_placement_source_exported": False,
        "prime_anchor_value_source_exported": False,
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
        "strict_damping_prime_slope_anchor_equivalence_certificate": {
            "prime_anchor_rows": anchor_rows,
            "candidate_grid_audit_rows": grid_rows,
            "proof_sketch": "On the P2527 line v_p=a log(p), the single equation v_p=(4/5)log(p) has determinant log(p)>0 for each audited prime p, hence a=4/5. This proves conditional sufficiency and placement equivalence, not a source for the anchor.",
            "source_obligation_normal_form": "beta_eta_numeric_source = multiplicative_character_law_source AND prime_log_proportionality_source AND slope_value_or_prime_anchor_source; strict_damping_beta_eta_source additionally requires m2_operator_signature_source",
        },
    }
    gatekeepers = {
        "p2527_inherited": theorem_export["p2527_slope_line_inherited"],
        "single_prime_anchors_conditionally_pin_slope": theorem_export["all_single_prime_anchors_recover_strict_delta"] and theorem_export["all_candidate_grid_audits_accept_only_strict_slope"],
        "anchor_placement_nonunique_not_source": theorem_export["single_prime_anchor_placement_nonunique"] and not theorem_export["slope_value_source_exported"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "prime_anchor_placement_source_exported",
            "prime_anchor_value_source_exported",
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
        "packet_id": "P2528",
        "stage_id": "S1478",
        "status": "STRICT_DAMPING_PRIME_SLOPE_ANCHOR_EQUIVALENCE_CERTIFICATE_CONDITIONAL_ANCHOR_ONLY_NO_ANCHOR_SOURCE_NO_SLOPE_SOURCE_NO_OPERATOR_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_prime_slope_anchor_equivalence_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_prime_slope_anchor_equivalence_certificate"]["theorem_export"]
    lines = [
        "# P2528/S1478 strict damping prime-slope anchor equivalence certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2527 slope line inherited: `{t['p2527_slope_line_inherited']}`.",
        f"- Prime anchor count: `{t['prime_anchor_count']}`.",
        f"- All prime anchors have positive determinant: `{t['all_prime_anchors_have_positive_log_determinant']}`.",
        f"- All single-prime anchors recover strict delta: `{t['all_single_prime_anchors_recover_strict_delta']}`.",
        f"- All grid audits accept only strict slope: `{t['all_candidate_grid_audits_accept_only_strict_slope']}`.",
        f"- Single-prime anchor placement nonunique: `{t['single_prime_anchor_placement_nonunique']}`.",
        f"- Slope value source exported: `{t['slope_value_source_exported']}`.",
        f"- Beta/eta numeric source exported: `{t['beta_eta_numeric_source_exported']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports only conditional equivalence/sufficiency of single-prime slope anchors on the P2527 slope line. It does not source anchor placement, anchor value, slope value, m=2 operator signature, bridge completion, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_prime_slope_anchor_equivalence_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_prime_slope_anchor_equivalence_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
