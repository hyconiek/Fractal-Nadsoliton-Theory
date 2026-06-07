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
OUT = GEN / "p2551_s1501_strict_damping_post_prime_log_slope_anchor_obstruction_certificate.json"
MD = GEN / "p2551_s1501_strict_damping_post_prime_log_slope_anchor_obstruction_certificate.md"

SOURCE_FILES = {
    "P2528_PRIME_SLOPE_ANCHOR_EQUIVALENCE": GEN / "p2528_s1478_strict_damping_prime_slope_anchor_equivalence_certificate.json",
    "P2543_SLOPE_VALUE_OBSTRUCTION": GEN / "p2543_s1493_strict_damping_slope_value_current_premise_obstruction_certificate.json",
    "P2547_POST_IDENTITY_RESIDUAL_TRIKEY": GEN / "p2547_s1497_strict_damping_post_identity_residual_trikey_certificate.json",
    "P2550_PRIME_LOG_ADJACENT_RATIO_BASIS": GEN / "p2550_s1500_strict_damping_prime_log_adjacent_ratio_basis_certificate.json",
}

PRIMES = [2, 3, 5, 7, 11]
STRICT_DELTA = Fraction(4, 5)
SLOPE_CANDIDATES = [Fraction(-1, 1), Fraction(0, 1), Fraction(1, 2), Fraction(4, 5), Fraction(1, 1), Fraction(9, 5), Fraction(2, 1)]
NEGATIVE_EXPORT_FLAGS = [
    "slope_value_or_prime_anchor_source_exported",
    "prime_log_proportionality_source_exported",
    "m2_operator_signature_source_exported",
    "strict_quadruple_trace_source_exported",
    "beta_eta_numeric_source_exported",
    "strict_damping_beta_eta_source_exported",
    "source_obligation_discharge_exported",
    "damping_compression_bridge_component_ready",
    "full_bridge_theorem_exported",
    "role_transfer_theorem_exported",
    "selector_closure_exported",
    "qw2191_discharged_by_this_certificate",
    "role_bearing_ltotal_exported",
    "toe_closure_claimed",
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
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2551|S1501|post-prime-log slope anchor|slope anchor obstruction|post prime log slope",
        "intended_research_nonduplication": "post-prime-log.*slope|prime-log.*slope anchor.*obstruction|adjacent ratio.*slope|slope anchor.*after.*prime-log|one-dimensional slope line",
        "slope_precursors": "P2528|S1478|P2543|S1493|slope_value_or_prime_anchor_source|prime slope anchor|delta=4/5",
        "prime_log_frontier": "P2550|S1500|adjacent ratio basis|prime_log_proportionality_source|normalized prime ratios",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def load_theorems(sources: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {
        "P2528_PRIME_SLOPE_ANCHOR_EQUIVALENCE": theorem(sources["P2528_PRIME_SLOPE_ANCHOR_EQUIVALENCE"], "strict_damping_prime_slope_anchor_equivalence_certificate"),
        "P2543_SLOPE_VALUE_OBSTRUCTION": theorem(sources["P2543_SLOPE_VALUE_OBSTRUCTION"], "strict_damping_slope_value_current_premise_obstruction_certificate"),
        "P2547_POST_IDENTITY_RESIDUAL_TRIKEY": theorem(sources["P2547_POST_IDENTITY_RESIDUAL_TRIKEY"], "strict_damping_post_identity_residual_trikey_certificate"),
        "P2550_PRIME_LOG_ADJACENT_RATIO_BASIS": theorem(sources["P2550_PRIME_LOG_ADJACENT_RATIO_BASIS"], "strict_damping_prime_log_adjacent_ratio_basis_certificate"),
    }


def frac_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def slope_line_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for slope in SLOPE_CANDIDATES:
        ratios = [slope for _ in PRIMES]
        prime_values_symbolic = [f"({frac_text(slope)})*log({prime})" for prime in PRIMES]
        prime_values_numeric = [float(slope) * math.log(prime) for prime in PRIMES]
        rows.append({
            "slope_delta": frac_text(slope),
            "eta_if_slope_delta": frac_text(Fraction(1, 1) + slope),
            "normalized_ratios_all_equal": True,
            "satisfies_p2550_adjacent_ratio_basis": True,
            "prime_values_v_p_symbolic": dict(zip([f"v_{p}" for p in PRIMES], prime_values_symbolic)),
            "prime_values_v_p_numeric": dict(zip([f"v_{p}" for p in PRIMES], prime_values_numeric)),
            "strict_delta_4_over_5_anchor_accepts": slope == STRICT_DELTA,
            "slope_value_or_prime_anchor_source_satisfied_as_target": slope == STRICT_DELTA,
            "post_prime_log_slope_countermodel": slope != STRICT_DELTA,
            "distance_from_strict_delta": frac_text(slope - STRICT_DELTA),
        })
    return rows


def single_prime_anchor_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for prime in PRIMES:
        strict_anchor = float(STRICT_DELTA) * math.log(prime)
        rows.append({
            "anchor_prime": prime,
            "anchor_equation": f"v_{prime}=(4/5)*log({prime})",
            "anchor_value_numeric": strict_anchor,
            "under_prime_log_line_vp_delta_logp_solves_delta": "4/5",
            "determinant_log_prime_nonzero": math.log(prime) != 0.0,
            "single_prime_anchor_suffices_given_prime_log_line": True,
        })
    return rows


def obstruction_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    accepting = [row for row in rows if row["strict_delta_4_over_5_anchor_accepts"]]
    countermodels = [row for row in rows if row["post_prime_log_slope_countermodel"]]
    return {
        "post_prime_log_slope_line_parameter": "delta in Q/R with v_p=delta*log(p)",
        "audited_slope_count": len(rows),
        "strict_accepting_slopes": [row["slope_delta"] for row in accepting],
        "countermodel_slopes": [row["slope_delta"] for row in countermodels],
        "countermodel_count": len(countermodels),
        "explicit_delta_one_half_countermodel": next(row for row in rows if row["slope_delta"] == "1/2"),
        "prime_log_adjacent_basis_entails_slope_value": len(countermodels) == 0,
        "single_prime_anchor_rows": single_prime_anchor_rows(),
        "single_prime_anchor_count": len(PRIMES),
        "all_single_prime_anchors_suffice_given_prime_log_line": True,
        "model_separation_statement": (
            "After prime-log proportionality, the remaining numeric freedom is a one-dimensional slope line. "
            "The exact value delta=4/5 is selected by any one strict prime anchor, but is not entailed by the prime-log basis itself."
        ),
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    source_theorems = load_theorems(sources)
    rows = slope_line_rows()
    audit = obstruction_audit(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2551_T1_strict_damping_post_prime_log_slope_anchor_obstruction_certificate",
        "audited_chain": ["P2528/S1478", "P2543/S1493", "P2547/S1497", "P2550/S1500"],
        "frontier_source_under_attack": "slope_value_or_prime_anchor_source",
        "p2528_prime_slope_anchor_equivalence_inherited": bool(source_theorems["P2528_PRIME_SLOPE_ANCHOR_EQUIVALENCE"]),
        "p2543_slope_value_current_premise_obstruction_inherited": source_theorems["P2543_SLOPE_VALUE_OBSTRUCTION"].get("current_slope_line_premises_do_not_entail_delta_4_over_5") is True,
        "p2547_post_identity_residual_trikey_inherited": source_theorems["P2547_POST_IDENTITY_RESIDUAL_TRIKEY"].get("post_identity_residual_trikey_irredundancy_exported") is True,
        "p2550_prime_log_adjacent_basis_inherited": source_theorems["P2550_PRIME_LOG_ADJACENT_RATIO_BASIS"].get("full_adjacent_ratio_basis_suffices_for_prime_log_proportionality") is True,
        "post_prime_log_slope_line_rows": rows,
        "post_prime_log_slope_anchor_obstruction_audit": audit,
        "post_prime_log_slope_value_nonentailment_exported": not audit["prime_log_adjacent_basis_entails_slope_value"],
        "single_prime_anchor_equivalence_reconfirmed_given_prime_log_line": audit["all_single_prime_anchors_suffice_given_prime_log_line"],
        "slope_value_or_prime_anchor_source_exported": False,
        "residual_after_hypothetical_identity_m2_prime_log_and_slope": [],
        "conditional_full_strict_damping_if_all_other_sources_and_slope_anchor_supplied": True,
        "recommended_next_honest_step": (
            "Do not infer delta=4/5 from prime-log proportionality. The next honest step is to seek a strict nadsoliton prime-value/slope-anchor theorem "
            "that exports one equation v_p=(4/5)log(p) or an equivalent slope selector; only after that can the strict damping numeric source be combined with the still-separate m2/trace and bridge/role-transfer gates."
        ),
        "not_licensed": [
            "No strict mechanism selecting delta=4/5 is supplied.",
            "No slope-value or prime-anchor source theorem is exported.",
            "No strict damping beta/eta source, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
    }
    for flag in NEGATIVE_EXPORT_FLAGS:
        theorem_export[flag] = False
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2543_slope_obstruction_inherited": theorem_export["p2543_slope_value_current_premise_obstruction_inherited"],
        "p2550_prime_log_basis_inherited": theorem_export["p2550_prime_log_adjacent_basis_inherited"],
        "post_prime_log_countermodels_exist": audit["countermodel_count"] > 0,
        "delta_one_half_is_explicit_countermodel": audit["explicit_delta_one_half_countermodel"]["post_prime_log_slope_countermodel"] is True,
        "single_prime_anchors_suffice_conditionally": theorem_export["single_prime_anchor_equivalence_reconfirmed_given_prime_log_line"],
        "no_false_slope_source_export": theorem_export["slope_value_or_prime_anchor_source_exported"] is False,
        "no_false_bridge_or_role_transfer": theorem_export["full_bridge_theorem_exported"] is False and theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_or_toe_claim": theorem_export["qw2191_discharged_by_this_certificate"] is False and theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2551",
        "stage_id": "S1501",
        "status": "STRICT_DAMPING_POST_PRIME_LOG_SLOPE_ANCHOR_OBSTRUCTION_CERTIFICATE_NO_SLOPE_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_post_prime_log_slope_anchor_obstruction_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_post_prime_log_slope_anchor_obstruction_certificate"]["theorem_export"]
    audit = t["post_prime_log_slope_anchor_obstruction_audit"]
    lines = [
        "# P2551/S1501 strict damping post-prime-log slope-anchor obstruction certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier source under attack: `{t['frontier_source_under_attack']}`.",
        f"- P2543 slope obstruction inherited: `{t['p2543_slope_value_current_premise_obstruction_inherited']}`.",
        f"- P2550 prime-log adjacent basis inherited: `{t['p2550_prime_log_adjacent_basis_inherited']}`.",
        f"- Audited slope count: `{audit['audited_slope_count']}`.",
        f"- Strict accepting slopes: `{audit['strict_accepting_slopes']}`.",
        f"- Countermodel slopes: `{audit['countermodel_slopes']}`.",
        f"- Explicit delta=1/2 countermodel passes prime-log basis: `{audit['explicit_delta_one_half_countermodel']['satisfies_p2550_adjacent_ratio_basis']}`.",
        f"- Prime-log adjacent basis entails slope value: `{audit['prime_log_adjacent_basis_entails_slope_value']}`.",
        f"- Single prime anchors suffice given prime-log line: `{audit['all_single_prime_anchors_suffice_given_prime_log_line']}`.",
        f"- Slope/prime-anchor source exported: `{t['slope_value_or_prime_anchor_source_exported']}`.",
        "",
        "## Interpretation",
        "",
        audit["model_separation_statement"],
        "",
        "The packet reconfirms that any one strict prime-value anchor would fix the slope after prime-log proportionality, but the existence of such an anchor is still the missing source theorem.",
        "",
        "## Recommended next honest step",
        "",
        t["recommended_next_honest_step"],
        "",
        "## Negative controls",
        "",
        "No slope-value source, prime-anchor source, strict damping beta/eta source, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, physical-value generator, or ToE closure is exported.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_post_prime_log_slope_anchor_obstruction_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2551/S1501` audits the post-prime-log residual slope key.  Even after the P2550 adjacent-ratio basis collapses the five prime ratios to a single line `v_p=delta log(p)`, the value of `delta` remains free: audited slopes such as `delta=1/2` satisfy the prime-log basis but are not the strict `delta=4/5` target.  Any one strict prime anchor `v_p=(4/5)log(p)` would fix the line, but P2551 exports no theorem that strict nadsoliton dynamics supplies such an anchor.
""".strip()
    lag_section = """
`P2551/S1501` blocks promotion of the post-prime-log slope value into a role-bearing `L_total` damping term.  The prime-log line plus a single prime anchor would determine `delta=4/5`, but the anchor remains a source obligation rather than a strict theorem.  No strict damping beta/eta source, bridge completion, or role-transfer license is exported.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2551/S1501 post-prime-log slope-anchor obstruction guard", "## P2551/S1501 post-prime-log slope-anchor obstruction guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2551/S1501 post-prime-log slope-anchor Ltotal guard", "## P2551/S1501 post-prime-log slope-anchor Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
