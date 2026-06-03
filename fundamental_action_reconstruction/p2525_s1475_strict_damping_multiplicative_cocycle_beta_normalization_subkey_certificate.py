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
OUT = GEN / "p2525_s1475_strict_damping_multiplicative_cocycle_beta_normalization_subkey_certificate.json"
MD = GEN / "p2525_s1475_strict_damping_multiplicative_cocycle_beta_normalization_subkey_certificate.md"

SOURCE_FILES = {
    "P2524_AFFINE_CONSISTENCY_CONTINUUM": GEN / "p2524_s1474_strict_damping_affine_consistency_continuum_nonidentifiability_certificate.json",
}

NODE_DOMAIN = list(range(1, 12))
INTERCEPT_CANDIDATES = [Fraction(-1, 1), Fraction(-1, 2), Fraction(0, 1), Fraction(1, 2), Fraction(1, 1)]
SLOPE_CANDIDATES = [Fraction(-1, 1), Fraction(0, 1), Fraction(1, 2), Fraction(4, 5), Fraction(1, 1), Fraction(9, 5), Fraction(2, 1)]
STRICT_INTERCEPT = Fraction(0, 1)
STRICT_DELTA = Fraction(4, 5)
STRICT_ETA = Fraction(9, 5)
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
    return {"count": len(lines), "samples": lines[:40]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2525|S1475|multiplicative cocycle beta normalization|semigroup character beta normalization|multiplicative character subkey|beta-normalization cocycle",
        "intended_research_nonduplication": "multiplicative cocycle|multiplicative character|semigroup character|additive log character|log-character|y\\(de\\)|beta normalization.*multiplicative",
        "precursor_packets": "P2524|S1474|affine consistency continuum|P2523|pairwise secant consistency|P2520|endpoint anchor subkey lattice",
        "subkey_language": "beta normalization anchor|left normalization|slope value source|numeric-source subkey|proper subset obstruction",
        "source_blockers": "beta_eta_numeric_source|strict damping source theorem|m2_operator_signature_source|node data source|anchor basis source",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def frac_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def ell(d: int) -> float:
    return math.log(d)


def y_affine(d: int, intercept: Fraction, slope: Fraction) -> float:
    return float(intercept) + float(slope) * ell(d)


def multiplicative_pairs() -> list[tuple[int, int, int]]:
    return [(d, e, d * e) for d, e in product(NODE_DOMAIN, repeat=2) if d * e in NODE_DOMAIN]


def multiplicative_cocycle_rows() -> list[dict[str, Any]]:
    pairs = multiplicative_pairs()
    rows = []
    for intercept in INTERCEPT_CANDIDATES:
        for slope in SLOPE_CANDIDATES:
            defects = [
                y_affine(de, intercept, slope) - y_affine(d, intercept, slope) - y_affine(e, intercept, slope)
                for d, e, de in pairs
            ]
            max_abs_defect = max(abs(value) for value in defects)
            rows.append({
                "intercept_log_beta": frac_text(intercept),
                "slope_delta": frac_text(slope),
                "eta_if_slope_delta": frac_text(Fraction(1, 1) + slope),
                "is_strict_numeric_target": intercept == STRICT_INTERCEPT and slope == STRICT_DELTA,
                "multiplicative_pair_count": len(pairs),
                "defect_min": min(defects),
                "defect_max": max(defects),
                "defect_mean": mean(defects),
                "defect_population_std": pstdev(defects),
                "max_abs_multiplicative_cocycle_defect": max_abs_defect,
                "defect_equals_minus_intercept_on_all_pairs": all(abs(value + float(intercept)) < TOL for value in defects),
                "multiplicative_character_accepts": max_abs_defect < TOL,
                "accepted_intercept_is_beta_normalization": max_abs_defect < TOL and intercept == STRICT_INTERCEPT,
            })
    return rows


def unital_rows() -> list[dict[str, Any]]:
    return [{
        "intercept_log_beta": frac_text(intercept),
        "y_1": float(intercept),
        "unital_character_y_1_zero_accepts": abs(float(intercept)) < TOL,
        "accepted_intercept_is_beta_normalization": intercept == STRICT_INTERCEPT,
    } for intercept in INTERCEPT_CANDIDATES]


def multiplicative_filter(rows: list[dict[str, Any]]) -> dict[str, Any]:
    accepted = [row for row in rows if row["multiplicative_character_accepts"]]
    accepted_strict_slope = [row for row in accepted if row["slope_delta"] == "4/5"]
    return {
        "candidate_grid_row_count": len(rows),
        "multiplicative_character_accepting_row_count": len(accepted),
        "accepted_intercepts": sorted({row["intercept_log_beta"] for row in accepted}),
        "accepted_slopes": [row["slope_delta"] for row in accepted],
        "multiplicative_law_pins_log_beta_zero_on_grid": sorted({row["intercept_log_beta"] for row in accepted}) == ["0"],
        "multiplicative_law_leaves_slope_continuum_on_grid": len({row["slope_delta"] for row in accepted}) == len(SLOPE_CANDIDATES),
        "strict_slope_filter_after_multiplicative_law_count": len(accepted_strict_slope),
        "strict_slope_filter_after_multiplicative_law_unique_target": len(accepted_strict_slope) == 1 and accepted_strict_slope[0]["is_strict_numeric_target"],
        "source_obligation_normal_form": "beta_eta_numeric_source = multiplicative_character_law_source AND slope_value_source; strict_damping_beta_eta_source additionally requires m2_operator_signature_source",
    }


def append_doc_sections() -> None:
    eq_section = """
`P2525/S1475` adds a conditional beta-normalization subkey audit.  For finite node data of affine form `y_d=b+a log d`, the multiplicative character law `y_{de}=y_d+y_e` on audited products `de<=11` has defect `-b`; hence it pins `b=log beta=0` and recovers the left-normalization subkey if that law is supplied.  It deliberately leaves the slope continuum `a` untouched, so `eta=9/5` still requires an independent slope-value source and the strict damping source still requires the `m=2` operator-signature source.
""".strip()
    lag_section = """
`P2525/S1475` narrows the numeric-source obligation by showing that a sourced multiplicative log-character law would supply the beta-normalization subkey `log beta=0`.  This is only a conditional subkey: the law itself is not sourced here, the slope `delta=4/5` is not selected, and no nonlinear compression-flow source or role-bearing `L_total` term is licensed.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "P2525/S1475", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "P2525/S1475", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2524 = theorem(sources["P2524_AFFINE_CONSISTENCY_CONTINUUM"], "strict_damping_affine_consistency_continuum_nonidentifiability_certificate")
    rows = multiplicative_cocycle_rows()
    unital = unital_rows()
    filt = multiplicative_filter(rows)
    pairs = multiplicative_pairs()
    theorem_export = {
        "frontier_atom_under_attack": "beta_normalization_left_anchor_source_subkey_after_affine_consistency_continuum_nonidentifiability",
        "p2524_continuum_nonidentifiability_inherited": bool(p2524.get("basis_independent_affine_consistency_nonidentifiability_exported", False)),
        "multiplicative_pair_count_on_domain_1_to_11": len(pairs),
        "multiplicative_pairs": [{"d": d, "e": e, "de": de} for d, e, de in pairs],
        "candidate_grid_row_count": filt["candidate_grid_row_count"],
        "multiplicative_character_accepting_row_count": filt["multiplicative_character_accepting_row_count"],
        "multiplicative_law_pins_log_beta_zero_on_grid": filt["multiplicative_law_pins_log_beta_zero_on_grid"],
        "multiplicative_law_leaves_slope_continuum_on_grid": filt["multiplicative_law_leaves_slope_continuum_on_grid"],
        "strict_slope_filter_after_multiplicative_law_unique_target": filt["strict_slope_filter_after_multiplicative_law_unique_target"],
        "conditional_beta_normalization_subkey_exported": True,
        "multiplicative_character_law_source_exported": False,
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
        "strict_damping_multiplicative_cocycle_beta_normalization_subkey_certificate": {
            "multiplicative_cocycle_rows": rows,
            "unital_character_rows": unital,
            "multiplicative_filter_audit": filt,
            "proof_sketch": "For y_d=b+a log(d), y_{de}-y_d-y_e = b+a log(de)-b-a log(d)-b-a log(e) = -b. Therefore the multiplicative character law pins b=0 but cancels all dependence on a.",
        },
    }
    gatekeepers = {
        "p2524_inherited": theorem_export["p2524_continuum_nonidentifiability_inherited"],
        "multiplicative_cocycle_pins_only_intercept": theorem_export["multiplicative_law_pins_log_beta_zero_on_grid"] and theorem_export["multiplicative_law_leaves_slope_continuum_on_grid"],
        "strict_target_requires_extra_slope_filter": theorem_export["strict_slope_filter_after_multiplicative_law_unique_target"] and theorem_export["multiplicative_character_accepting_row_count"] > 1,
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "multiplicative_character_law_source_exported",
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
        "packet_id": "P2525",
        "stage_id": "S1475",
        "status": "STRICT_DAMPING_MULTIPLICATIVE_COCYCLE_BETA_NORMALIZATION_SUBKEY_CERTIFICATE_CONDITIONAL_SUBKEY_ONLY_NO_SLOPE_SOURCE_NO_OPERATOR_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_multiplicative_cocycle_beta_normalization_subkey_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_multiplicative_cocycle_beta_normalization_subkey_certificate"]["theorem_export"]
    lines = [
        "# P2525/S1475 strict damping multiplicative-cocycle beta-normalization subkey certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2524 continuum nonidentifiability inherited: `{t['p2524_continuum_nonidentifiability_inherited']}`.",
        f"- Multiplicative product pairs on `d=1..11`: `{t['multiplicative_pair_count_on_domain_1_to_11']}`.",
        f"- Candidate grid row count: `{t['candidate_grid_row_count']}`.",
        f"- Multiplicative-character accepting rows: `{t['multiplicative_character_accepting_row_count']}`.",
        f"- Multiplicative law pins `log beta=0`: `{t['multiplicative_law_pins_log_beta_zero_on_grid']}`.",
        f"- Multiplicative law leaves slope continuum: `{t['multiplicative_law_leaves_slope_continuum_on_grid']}`.",
        f"- Conditional beta-normalization subkey exported: `{t['conditional_beta_normalization_subkey_exported']}`.",
        f"- Beta/eta numeric source exported: `{t['beta_eta_numeric_source_exported']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports only a conditional beta-normalization subkey from an assumed multiplicative log-character law. It does not source that law, does not source the slope value, does not source the m=2 operator signature, and exports no bridge completion, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_multiplicative_cocycle_beta_normalization_subkey_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_multiplicative_cocycle_beta_normalization_subkey_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
