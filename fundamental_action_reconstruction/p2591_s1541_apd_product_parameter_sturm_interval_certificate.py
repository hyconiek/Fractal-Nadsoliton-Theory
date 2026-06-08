#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from typing import Any

import sympy as sp

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2590_s1540_apd_finite_even_moment_shell_interval_nonuniqueness_audit import (
    EXPECTED_CARDINALITY,
    INTERNAL_FOURTH_SHELL,
    INTERNAL_SECOND_SHELL,
    INTERNAL_SIXTH_SHELL,
    OUT as P2590_OUT,
    PRODUCT_PARAMETER_GRID,
    shell_support_witness,
    support_from_product,
    theorem,
)

GEN = ROOT / "generated"
OUT = GEN / "p2591_s1541_apd_product_parameter_sturm_interval_certificate.json"
MD = GEN / "p2591_s1541_apd_product_parameter_sturm_interval_certificate.md"

SOURCE_FILES = {
    "P2590_FINITE_EVEN_MOMENT_SHELL_INTERVAL": P2590_OUT,
}
INTERVAL = [300, 576]
ROOT_BOUND = 821
NEGATIVE_EXPORT_FLAGS = [
    "apd_product_parameter_interval_source_exported",
    "apd_sturm_interval_selector_source_exported",
    "apd_finite_even_moment_shell_source_exported",
    "apd_support_selection_source_exported",
    "apd_finite_support_source_exported",
    "apd_positive_measure_source_exported",
    "strict_dynamical_source_for_A_P_D_exported",
    "strict_phase_frequency_source_exported",
    "strict_damping_beta_eta_source_exported",
    "bridge_theorem_exported",
    "legacy_to_strict_completion_bridge_exported",
    "role_transfer_theorem_exported",
    "role_bearing_ltotal_exported",
    "qw2191_discharged_by_this_certificate",
    "toe_closure_claimed",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
    ], cwd=REPO, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2591|S1541|APD product parameter interval Sturm|product parameter interval Sturm.*APD",
        "intended_research_nonduplication": "APD.*Sturm interval|Sturm interval.*APD|APD.*discriminant interval|discriminant interval.*APD|APD.*even moment shell interval proof|even moment shell interval proof.*APD|APD.*quartic support interval|quartic support interval.*APD",
        "apd_precursors": "P2590|S1540|APD.*finite even moment shell interval|strict_dynamical_source_for_A_P_D",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def sign_at(poly: sp.Expr, symbol: sp.Symbol, value: int) -> int:
    evaluated = sp.expand(poly).subs(symbol, sp.Integer(value))
    return int(sp.sign(evaluated))


def variations(signs: list[int]) -> int:
    filtered = [sign for sign in signs if sign != 0]
    return sum(1 for left, right in zip(filtered, filtered[1:]) if left * right < 0)


def sturm_count(poly: sp.Expr, symbol: sp.Symbol, left: int, right: int) -> dict[str, Any]:
    sequence = [sp.expand(item) for item in sp.sturm(poly, symbol)]
    left_signs = [sign_at(item, symbol, left) for item in sequence]
    right_signs = [sign_at(item, symbol, right) for item in sequence]
    left_variations = variations(left_signs)
    right_variations = variations(right_signs)
    return {
        "interval": [left, right],
        "sturm_sequence_degrees": [int(sp.degree(item, gen=symbol)) for item in sequence],
        "left_signs": left_signs,
        "right_signs": right_signs,
        "left_variations": left_variations,
        "right_variations": right_variations,
        "root_count_open_interval": left_variations - right_variations,
    }


def support_snapshot(product_parameter: int) -> dict[str, Any]:
    witness = shell_support_witness(support_from_product(float(product_parameter)))
    return {
        "product_parameter": product_parameter,
        "squared_offsets": witness["squared_offsets"],
        "points": witness["points"],
        "has_expected_internal_second_fourth_sixth_shell": witness["has_expected_internal_second_fourth_sixth_shell"],
        "vandermonde_rank": witness["vandermonde_rank"],
        "max_abs_recovered_weight_error": witness["max_abs_recovered_weight_error"],
    }


def sturm_interval_certificate() -> dict[str, Any]:
    x, e4 = sp.symbols("x e4")
    offset_quartic = x ** 4 - 30 * x ** 3 + 273 * x ** 2 - 820 * x + e4
    discriminant = sp.factor(sp.discriminant(offset_quartic, x))
    discriminant_core = sp.factor(discriminant / 16)
    no_discriminant_roots = sturm_count(discriminant_core, e4, INTERVAL[0], INTERVAL[1])
    anchor_poly = sp.expand(offset_quartic.subs(e4, INTERVAL[0]))
    positive_roots_anchor = sturm_count(anchor_poly, x, 0, ROOT_BOUND)
    negative_roots_anchor = sturm_count(anchor_poly, x, -ROOT_BOUND, 0)
    endpoint_snapshots = [support_snapshot(value) for value in INTERVAL]
    return {
        "sympy_version": sp.__version__,
        "offset_quartic_family": "x^4 - 30*x^3 + 273*x^2 - 820*x + e4",
        "discriminant": str(discriminant),
        "discriminant_core": str(discriminant_core),
        "certified_product_parameter_interval": INTERVAL,
        "root_bound_used_for_anchor_counts": ROOT_BOUND,
        "discriminant_sturm_count_on_interval": no_discriminant_roots,
        "discriminant_positive_at_left_endpoint": bool(discriminant_core.subs(e4, INTERVAL[0]) > 0),
        "discriminant_positive_at_right_endpoint": bool(discriminant_core.subs(e4, INTERVAL[1]) > 0),
        "positive_root_count_anchor_e4_300": positive_roots_anchor,
        "negative_root_count_anchor_e4_300": negative_roots_anchor,
        "positive_real_root_count_is_constant_on_interval": no_discriminant_roots["root_count_open_interval"] == 0 and positive_roots_anchor["root_count_open_interval"] == 4 and negative_roots_anchor["root_count_open_interval"] == 0,
        "endpoint_support_snapshots": endpoint_snapshots,
        "vieta_power_sums_fixed_on_interval": {
            "p1_internal_second_shell": INTERNAL_SECOND_SHELL,
            "p2_internal_fourth_shell": INTERNAL_FOURTH_SHELL,
            "p3_internal_sixth_shell": INTERNAL_SIXTH_SHELL,
        },
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2590_payload = load_json(SOURCE_FILES["P2590_FINITE_EVEN_MOMENT_SHELL_INTERVAL"])
    p2590 = theorem(p2590_payload, "apd_finite_even_moment_shell_interval_nonuniqueness_audit")
    certificate = sturm_interval_certificate()
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2591_T1_apd_product_parameter_sturm_interval_certificate",
        "audited_chain": ["P2590/S1540"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "promote finite even-moment shell interval existence into an APD support selector/source",
        "p2590_finite_even_moment_shell_interval_inherited": p2590.get("finite_even_moment_shell_prefix_does_not_select_apd_support") is True,
        "apd_product_parameter_sturm_interval_certificate": certificate,
        "continuous_interval_of_valid_supports_certified": certificate["positive_real_root_count_is_constant_on_interval"],
        "finite_even_moment_shell_interval_still_does_not_select_apd_support": True,
        "fixed_support_full_moments_remain_conditional": True,
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Do not promote the Sturm-certified product-parameter interval into an APD support source. P2591 proves a continuous interval of valid mirror finite supports for the same finite even-moment shell prefix, so the next honest step is a strict nadsoliton-derived density/support law or an internal APD dynamics theorem, not another finite moment-shell selector."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2590_finite_even_moment_shell_interval_inherited": theorem_export["p2590_finite_even_moment_shell_interval_inherited"],
        "discriminant_has_no_roots_on_interval": certificate["discriminant_sturm_count_on_interval"]["root_count_open_interval"] == 0,
        "discriminant_positive_at_endpoints": certificate["discriminant_positive_at_left_endpoint"] and certificate["discriminant_positive_at_right_endpoint"],
        "anchor_has_four_positive_roots": certificate["positive_root_count_anchor_e4_300"]["root_count_open_interval"] == 4,
        "anchor_has_zero_negative_roots": certificate["negative_root_count_anchor_e4_300"]["root_count_open_interval"] == 0,
        "continuous_interval_of_valid_supports_certified": theorem_export["continuous_interval_of_valid_supports_certified"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2591",
        "stage_id": "S1541",
        "status": "P2591_APD_PRODUCT_PARAMETER_STURM_INTERVAL_CERTIFICATE_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_product_parameter_sturm_interval_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2590_FINITE_EVEN_MOMENT_SHELL_INTERVAL": sha256_json(p2590_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_product_parameter_sturm_interval_certificate"]["theorem_export"]
    certificate = t["apd_product_parameter_sturm_interval_certificate"]
    lines = [
        "# P2591/S1541 APD product-parameter Sturm interval certificate", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Certified product-parameter interval: `{certificate['certified_product_parameter_interval']}`.",
        f"- Discriminant roots in interval: `{certificate['discriminant_sturm_count_on_interval']['root_count_open_interval']}`.",
        f"- Anchor positive roots: `{certificate['positive_root_count_anchor_e4_300']['root_count_open_interval']}`.",
        f"- Continuous interval certified: `{t['continuous_interval_of_valid_supports_certified']}`.",
        f"- Strict APD dynamic source exported: `{t['strict_dynamical_source_for_A_P_D_exported']}`.", "",
        "## Interpretation", "",
        "P2591 turns P2590's grid into a Sturm/discriminant interval certificate.  The offset quartic has no discriminant zero on `[300, 576]`, and an anchor Sturm count gives four positive roots and no negative roots; therefore the same finite even-moment shell prefix admits a continuous product-parameter family of mirror supports.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No product-parameter interval source, Sturm selector source, finite even-moment shell source, support-selection source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_product_parameter_sturm_interval_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2591/S1541 APD product-parameter Sturm interval certificate guard

`P2591/S1541` upgrades P2590's finite grid to a Sturm/discriminant interval certificate: the offset quartic has no discriminant root on `[300, 576]`, and the anchor count has four positive roots, so the same finite even-moment shell prefix supports a continuous mirror-support family.  This strengthens nonuniqueness; it still does not source the APD support law.
""".strip()
    lag_section = """
## P2591/S1541 APD product-parameter Sturm interval certificate Ltotal guard

`P2591/S1541` blocks a role-bearing APD Gram term in `L_total` from being justified by a Sturm-certified finite moment-shell interval.  The interval proves selector freedom, not a strict support/density source; strict nadsoliton dynamics must supply the actual support law.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2591/S1541 APD product-parameter Sturm interval certificate guard", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2591/S1541 APD product-parameter Sturm interval certificate Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
