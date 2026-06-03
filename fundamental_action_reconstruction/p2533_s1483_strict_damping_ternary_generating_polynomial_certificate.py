#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from math import factorial
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
OUT = GEN / "p2533_s1483_strict_damping_ternary_generating_polynomial_certificate.json"
MD = GEN / "p2533_s1483_strict_damping_ternary_generating_polynomial_certificate.md"

SOURCE_FILES = {
    "P2532_STRICTIZATION_DISTANCE": GEN / "p2532_s1482_strict_damping_four_key_strictization_distance_certificate.json",
}

SOURCE_KEYS = [
    "M_multiplicative_character_law_source",
    "P_prime_log_proportionality_source",
    "A_slope_value_or_prime_anchor_source",
    "O_m2_operator_signature_source",
]
STATUS_SYMBOLS = {
    "absent": "u",
    "axiom": "a",
    "strict": "s",
}
POLY_VARIABLE_ORDER = ["u_absent", "a_axiom", "s_strict"]


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
        "new_packet": "P2533|S1483|ternary generating polynomial|strictization polynomial|source boundary polynomial|theorem deficit polynomial|upgrade-edge polynomial",
        "intended_research_nonduplication": "generating polynomial|deficit polynomial|upgrade-edge polynomial|source boundary polynomial|ternary polynomial|coefficient certificate",
        "precursor_packets": "P2532|S1482|strictization distance|one-step strictization graph|P2531|four-key axiom boundary",
        "source_obligation_language": "axiom-augmented|strict theorem|absent source|strict_damping_beta_eta_source|theorem-upgrade",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def multinomial4(absent_count: int, axiom_count: int, strict_count: int) -> int:
    if absent_count + axiom_count + strict_count != len(SOURCE_KEYS):
        return 0
    return factorial(len(SOURCE_KEYS)) // (factorial(absent_count) * factorial(axiom_count) * factorial(strict_count))


def monomial_label(absent_count: int, axiom_count: int, strict_count: int) -> str:
    parts = []
    for symbol, power in [("u", absent_count), ("a", axiom_count), ("s", strict_count)]:
        if power == 0:
            continue
        parts.append(symbol if power == 1 else f"{symbol}^{power}")
    return "1" if not parts else "*".join(parts)


def coefficient_row(absent_count: int, axiom_count: int, strict_count: int, coefficient: int) -> dict[str, Any]:
    return {
        "exponents_u_a_s": [absent_count, axiom_count, strict_count],
        "monomial": monomial_label(absent_count, axiom_count, strict_count),
        "coefficient": coefficient,
        "theorem_deficit": absent_count + axiom_count,
        "absent_source_gap": absent_count,
        "axiom_upgrade_distance_if_present": axiom_count if absent_count == 0 and axiom_count > 0 else None,
    }


def full_ternary_polynomial() -> list[dict[str, Any]]:
    rows = []
    n = len(SOURCE_KEYS)
    for absent_count in range(n + 1):
        for axiom_count in range(n - absent_count + 1):
            strict_count = n - absent_count - axiom_count
            rows.append(coefficient_row(absent_count, axiom_count, strict_count, multinomial4(absent_count, axiom_count, strict_count)))
    return rows


def filtered_polynomial(rows: list[dict[str, Any]], predicate: str) -> list[dict[str, Any]]:
    out = []
    for row in rows:
        absent_count, axiom_count, strict_count = row["exponents_u_a_s"]
        include = (
            predicate == "strict_accept" and absent_count == 0 and axiom_count == 0 and strict_count == len(SOURCE_KEYS)
        ) or (
            predicate == "present_axiom_augmented" and absent_count == 0 and axiom_count > 0
        ) or (
            predicate == "blocked_missing_key" and absent_count > 0
        )
        if include:
            out.append(row)
    return out


def evaluate_at_ones(rows: list[dict[str, Any]]) -> int:
    return sum(row["coefficient"] for row in rows)


def polynomial_by_stat(rows: list[dict[str, Any]], stat_key: str) -> dict[str, int]:
    totals: dict[str, int] = {}
    for row in rows:
        stat_value = row[stat_key]
        if stat_value is None:
            continue
        totals[str(stat_value)] = totals.get(str(stat_value), 0) + row["coefficient"]
    return dict(sorted(totals.items(), key=lambda item: int(item[0])))


def derivative_edge_counts(rows: list[dict[str, Any]]) -> dict[str, Any]:
    absent_edge_count = sum(row["coefficient"] * row["exponents_u_a_s"][0] for row in rows)
    axiom_edge_count = sum(row["coefficient"] * row["exponents_u_a_s"][1] for row in rows)
    key_edge_count = 2 * (3 ** (len(SOURCE_KEYS) - 1))
    return {
        "absent_source_theorem_introduction_edges_from_du_du": absent_edge_count,
        "axiom_to_strict_theorem_upgrade_edges_from_da_da": axiom_edge_count,
        "total_one_step_edges_from_degree_derivative": absent_edge_count + axiom_edge_count,
        "uniform_per_key_upgrade_edge_count": key_edge_count,
        "all_key_upgrade_edge_counts": {key: key_edge_count for key in SOURCE_KEYS},
    }


def build_polynomial_certificate(p2532_theorem: dict[str, Any]) -> dict[str, Any]:
    full = full_ternary_polynomial()
    strict_accept = filtered_polynomial(full, "strict_accept")
    present_axiom = filtered_polynomial(full, "present_axiom_augmented")
    blocked = filtered_polynomial(full, "blocked_missing_key")
    return {
        "variable_order": POLY_VARIABLE_ORDER,
        "compressed_universe_polynomial": "(u_absent + a_axiom + s_strict)^4",
        "strict_accept_polynomial": "s_strict^4",
        "present_axiom_augmented_polynomial": "(a_axiom + s_strict)^4 - s_strict^4",
        "blocked_missing_key_polynomial": "(u_absent + a_axiom + s_strict)^4 - (a_axiom + s_strict)^4",
        "full_coefficient_rows": full,
        "strict_accept_coefficient_rows": strict_accept,
        "present_axiom_augmented_coefficient_rows": present_axiom,
        "blocked_missing_key_coefficient_rows": blocked,
        "class_counts_from_polynomial": {
            "strict_accept": evaluate_at_ones(strict_accept),
            "present_non_strict_axiom_augmented": evaluate_at_ones(present_axiom),
            "blocked_missing_key": evaluate_at_ones(blocked),
        },
        "theorem_deficit_counts_from_polynomial": polynomial_by_stat(full, "theorem_deficit"),
        "present_axiom_upgrade_distance_counts_from_polynomial": polynomial_by_stat(present_axiom, "axiom_upgrade_distance_if_present"),
        "blocked_absent_source_gap_counts_from_polynomial": polynomial_by_stat(blocked, "absent_source_gap"),
        "derivative_edge_counts_from_polynomial": derivative_edge_counts(full),
        "p2532_reference_counts": {
            "class_counts": p2532_theorem["strict_damping_four_key_strictization_distance_certificate"]["class_counts"],
            "theorem_deficit_counts": p2532_theorem["theorem_deficit_counts"],
            "present_axiom_upgrade_distance_counts": p2532_theorem["present_axiom_upgrade_distance_counts"],
            "blocked_absent_source_gap_counts": p2532_theorem["blocked_absent_source_gap_counts"],
            "one_step_edge_type_counts": p2532_theorem["one_step_edge_type_counts"],
            "upgraded_key_counts": p2532_theorem["upgraded_key_counts"],
        },
    }


def append_doc_sections() -> None:
    eq_section = """
`P2533/S1483` compresses the P2531/P2532 four-key ternary table into exact generating polynomials.  With `u=absent`, `a=axiom`, and `s=strict`, the universe is `(u+a+s)^4`; strict acceptance is exactly `s^4`; non-strict axiom-present rows are `(a+s)^4-s^4`; and missing-key blocked rows are `(u+a+s)^4-(a+s)^4`.  Coefficient extraction reproduces class counts `1/15/65`, theorem-deficit counts `1,8,24,32,16`, axiom-upgrade counts `4,6,4,1`, absent-source-gap counts `32,24,8,1`, and derivative edge counts `108/108`.  This is a symbolic/combinatorial proof certificate for the source-boundary table, not a source theorem.
""".strip()
    lag_section = """
`P2533/S1483` records the polynomial normal form of the four-key strict damping source boundary.  The symbolic certificate proves that strict acceptance occupies only the `s^4` monomial, while all other monomials remain either non-strict axiom bookkeeping or missing-source blockers.  It exports no key source, beta-eta source, `m=2` operator source, nonlinear compression-flow source, role-bearing `L_total`, bridge completion, role-transfer theorem, QW-2191 discharge, or ToE closure.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "P2533/S1483", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "P2533/S1483", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2532 = theorem(sources["P2532_STRICTIZATION_DISTANCE"], "strict_damping_four_key_strictization_distance_certificate")
    cert = build_polynomial_certificate(p2532)
    edge_counts = cert["derivative_edge_counts_from_polynomial"]
    reference = cert["p2532_reference_counts"]
    theorem_export = {
        "frontier_atom_under_attack": "four_key_ternary_source_boundary_generating_polynomial",
        "p2532_strictization_distance_inherited": bool(p2532.get("strictization_distance_certificate_exported", False)),
        "source_keys": SOURCE_KEYS,
        "status_symbols": STATUS_SYMBOLS,
        "compressed_universe_polynomial": cert["compressed_universe_polynomial"],
        "strict_accept_polynomial": cert["strict_accept_polynomial"],
        "present_axiom_augmented_polynomial": cert["present_axiom_augmented_polynomial"],
        "blocked_missing_key_polynomial": cert["blocked_missing_key_polynomial"],
        "class_counts_from_polynomial": cert["class_counts_from_polynomial"],
        "class_counts_match_p2532": cert["class_counts_from_polynomial"] == reference["class_counts"],
        "theorem_deficit_counts_from_polynomial": cert["theorem_deficit_counts_from_polynomial"],
        "theorem_deficit_counts_match_p2532": cert["theorem_deficit_counts_from_polynomial"] == reference["theorem_deficit_counts"],
        "present_axiom_upgrade_distance_counts_from_polynomial": cert["present_axiom_upgrade_distance_counts_from_polynomial"],
        "present_axiom_upgrade_distance_counts_match_p2532": cert["present_axiom_upgrade_distance_counts_from_polynomial"] == reference["present_axiom_upgrade_distance_counts"],
        "blocked_absent_source_gap_counts_from_polynomial": cert["blocked_absent_source_gap_counts_from_polynomial"],
        "blocked_absent_source_gap_counts_match_p2532": cert["blocked_absent_source_gap_counts_from_polynomial"] == reference["blocked_absent_source_gap_counts"],
        "derivative_absent_axiom_edge_counts_from_polynomial": {
            "absent_source_theorem_introduction": edge_counts["absent_source_theorem_introduction_edges_from_du_du"],
            "axiom_to_strict_theorem_upgrade": edge_counts["axiom_to_strict_theorem_upgrade_edges_from_da_da"],
        },
        "derivative_edge_counts_match_p2532": {
            "absent_source_theorem_introduction": edge_counts["absent_source_theorem_introduction_edges_from_du_du"],
            "axiom_to_strict_theorem_upgrade": edge_counts["axiom_to_strict_theorem_upgrade_edges_from_da_da"],
        } == reference["one_step_edge_type_counts"],
        "uniform_key_upgrade_counts_from_polynomial": edge_counts["all_key_upgrade_edge_counts"],
        "uniform_key_upgrade_counts_match_p2532": edge_counts["all_key_upgrade_edge_counts"] == reference["upgraded_key_counts"],
        "strict_accept_support_is_single_s4_monomial": cert["strict_accept_coefficient_rows"] == [coefficient_row(0, 0, len(SOURCE_KEYS), 1)],
        "present_axiom_and_blocked_support_disjoint": not {tuple(row["exponents_u_a_s"]) for row in cert["present_axiom_augmented_coefficient_rows"]}.intersection({tuple(row["exponents_u_a_s"]) for row in cert["blocked_missing_key_coefficient_rows"]}),
        "polynomial_certificate_exported": True,
        "axiom_promotion_to_strict_exported": False,
        "multiplicative_character_law_source_exported": False,
        "prime_log_proportionality_source_exported": False,
        "slope_value_or_prime_anchor_source_exported": False,
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
        "strict_damping_ternary_generating_polynomial_certificate": cert,
    }
    gatekeepers = {
        "p2532_inherited": theorem_export["p2532_strictization_distance_inherited"],
        "polynomial_counts_match_p2532": theorem_export["class_counts_match_p2532"] and theorem_export["theorem_deficit_counts_match_p2532"] and theorem_export["present_axiom_upgrade_distance_counts_match_p2532"] and theorem_export["blocked_absent_source_gap_counts_match_p2532"],
        "derivative_edges_match_p2532": theorem_export["derivative_edge_counts_match_p2532"] and theorem_export["uniform_key_upgrade_counts_match_p2532"],
        "support_partition_verified": theorem_export["strict_accept_support_is_single_s4_monomial"] and theorem_export["present_axiom_and_blocked_support_disjoint"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "axiom_promotion_to_strict_exported",
            "multiplicative_character_law_source_exported",
            "prime_log_proportionality_source_exported",
            "slope_value_or_prime_anchor_source_exported",
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
        "packet_id": "P2533",
        "stage_id": "S1483",
        "status": "STRICT_DAMPING_TERNARY_GENERATING_POLYNOMIAL_CERTIFICATE_SYMBOLIC_SOURCE_BOUNDARY_ONLY_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_ternary_generating_polynomial_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_ternary_generating_polynomial_certificate"]["theorem_export"]
    lines = [
        "# P2533/S1483 strict damping ternary generating polynomial certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2532 strictization distance inherited: `{t['p2532_strictization_distance_inherited']}`.",
        f"- Universe polynomial: `{t['compressed_universe_polynomial']}`.",
        f"- Strict accept polynomial: `{t['strict_accept_polynomial']}`.",
        f"- Present axiom polynomial: `{t['present_axiom_augmented_polynomial']}`.",
        f"- Blocked missing-key polynomial: `{t['blocked_missing_key_polynomial']}`.",
        f"- Class counts from polynomial: `{t['class_counts_from_polynomial']}`.",
        f"- Deficit counts from polynomial: `{t['theorem_deficit_counts_from_polynomial']}`.",
        f"- Derivative edge counts from polynomial: `{t['derivative_absent_axiom_edge_counts_from_polynomial']}`.",
        f"- Strict accept support is single s^4 monomial: `{t['strict_accept_support_is_single_s4_monomial']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports only a symbolic generating-polynomial certificate for the conditional four-key source boundary. It does not source any key, promote axioms to strict theorems, export bridge completion, export a role-transfer theorem, discharge QW-2191, produce role-bearing L_total, or claim ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_ternary_generating_polynomial_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_ternary_generating_polynomial_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
