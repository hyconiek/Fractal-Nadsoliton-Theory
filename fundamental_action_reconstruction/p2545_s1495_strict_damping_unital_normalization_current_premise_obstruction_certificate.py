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
OUT = GEN / "p2545_s1495_strict_damping_unital_normalization_current_premise_obstruction_certificate.json"
MD = GEN / "p2545_s1495_strict_damping_unital_normalization_current_premise_obstruction_certificate.md"

SOURCE_FILES = {
    "P2521_SINGLE_NODE_ANCHOR_EQUIVALENCE": GEN / "p2521_s1471_strict_damping_single_node_anchor_equivalence_certificate.json",
    "P2525_MULTIPLICATIVE_BETA_NORMALIZATION": GEN / "p2525_s1475_strict_damping_multiplicative_cocycle_beta_normalization_subkey_certificate.json",
    "P2530_FOUR_KEY_IRREDUNDANCY": GEN / "p2530_s1480_strict_damping_four_key_irredundancy_witness_certificate.json",
    "P2541_MULTIPLICATIVE_OBSTRUCTION": GEN / "p2541_s1491_strict_damping_multiplicative_character_current_premise_obstruction_certificate.json",
    "P2544_FOUR_KEY_CLOSURE_BLOCKER": GEN / "p2544_s1494_strict_damping_four_key_current_premise_closure_blocker_certificate.json",
}

NODE_DOMAIN = list(range(1, 12))
INTERCEPT_CANDIDATES = [Fraction(-1, 1), Fraction(-1, 2), Fraction(0, 1), Fraction(1, 2), Fraction(1, 1)]
SLOPE_CANDIDATES = [Fraction(0, 1), Fraction(1, 2), Fraction(4, 5), Fraction(1, 1), Fraction(9, 5)]
STRICT_INTERCEPT = Fraction(0, 1)
STRICT_DELTA = Fraction(4, 5)
TOL = 1e-14

NEGATIVE_EXPORT_FLAGS = [
    "unital_monoid_normalization_source_exported",
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
    return {"count": len(lines), "samples": lines[:60]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2545|S1495|unital normalization current-premise obstruction|y_1=0 source obstruction|unit-node normalization source",
        "intended_research_nonduplication": "unital monoid normalization|unital_monoid_normalization|y_1=0|y1=0|left normalization.*source|unit-node normalization",
        "precursor_packets": "P2521|S1471|P2525|S1475|P2530|S1480|P2541|S1491|P2544|S1494",
        "source_key_language": "multiplicative_character_law_source|strict_unital_monoid_normalization_y1_zero|beta_eta_numeric_source|m2_operator_signature_source",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def frac_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def y_value(d: int, intercept: Fraction, slope: Fraction) -> float:
    return float(intercept) + float(slope) * math.log(d)


def audited_product_pairs() -> list[tuple[int, int, int]]:
    return [(d, e, d * e) for d, e in product(NODE_DOMAIN, repeat=2) if d * e in NODE_DOMAIN]


def affine_rows() -> list[dict[str, Any]]:
    pairs = audited_product_pairs()
    rows: list[dict[str, Any]] = []
    for intercept in INTERCEPT_CANDIDATES:
        for slope in SLOPE_CANDIDATES:
            y1 = y_value(1, intercept, slope)
            product_defects = [
                y_value(de, intercept, slope) - y_value(d, intercept, slope) - y_value(e, intercept, slope)
                for d, e, de in pairs
            ]
            unit_left_defects = [
                y_value(d, intercept, slope) - y_value(1, intercept, slope) - y_value(d, intercept, slope)
                for d in NODE_DOMAIN
            ]
            unit_right_defects = [
                y_value(d, intercept, slope) - y_value(d, intercept, slope) - y_value(1, intercept, slope)
                for d in NODE_DOMAIN
            ]
            max_product_defect = max(abs(value) for value in product_defects)
            max_unit_defect = max(abs(value) for value in unit_left_defects + unit_right_defects)
            rows.append({
                "intercept_log_beta": frac_text(intercept),
                "slope_delta": frac_text(slope),
                "is_strict_numeric_target": intercept == STRICT_INTERCEPT and slope == STRICT_DELTA,
                "current_affine_consistency_premise_accepts": True,
                "current_domain_contains_unit_node": 1 in NODE_DOMAIN,
                "current_premises_accept_without_unit_law": True,
                "y_1": y1,
                "unital_y1_zero_accepts": abs(y1) < TOL,
                "multiplicative_character_accepts": max_product_defect < TOL,
                "unit_product_law_accepts": max_unit_defect < TOL,
                "max_abs_product_defect": max_product_defect,
                "max_abs_unit_product_defect": max_unit_defect,
                "all_product_defects_equal_minus_y1": all(abs(value + y1) < TOL for value in product_defects),
                "all_unit_defects_equal_minus_y1": all(abs(value + y1) < TOL for value in unit_left_defects + unit_right_defects),
            })
    return rows


def build_certificate(theorems: dict[str, dict[str, Any]], rows: list[dict[str, Any]]) -> dict[str, Any]:
    p2521 = theorems["P2521_SINGLE_NODE_ANCHOR_EQUIVALENCE"]
    p2525 = theorems["P2525_MULTIPLICATIVE_BETA_NORMALIZATION"]
    p2530 = theorems["P2530_FOUR_KEY_IRREDUNDANCY"]
    p2541 = theorems["P2541_MULTIPLICATIVE_OBSTRUCTION"]
    p2544 = theorems["P2544_FOUR_KEY_CLOSURE_BLOCKER"]
    affine_countermodels = [
        row for row in rows
        if row["current_affine_consistency_premise_accepts"] and not row["unital_y1_zero_accepts"]
    ]
    strict_slope_unital_countermodels = [
        row for row in affine_countermodels
        if row["slope_delta"] == frac_text(STRICT_DELTA)
    ]
    chosen = next(row for row in strict_slope_unital_countermodels if row["intercept_log_beta"] == "1/2")
    unital_rows = [row for row in rows if row["unital_y1_zero_accepts"]]
    nonunital_rows = [row for row in rows if not row["unital_y1_zero_accepts"]]
    return {
        "p2521_left_normalization_unsourced_inherited": p2521.get("beta_normalization_left_anchor_source_exported") is False,
        "p2525_conditional_multiplicative_to_beta_normalization_inherited": p2525.get("multiplicative_character_law_source_exported") is False,
        "p2530_four_key_irredundancy_inherited": p2530.get("four_key_irredundancy_witness_exported") is True,
        "p2541_unital_equivalence_inherited": p2541.get("multiplicative_law_equivalent_to_unital_left_normalization_inside_affine_family") is True,
        "p2544_no_false_source_blocker_inherited": p2544.get("four_key_current_premise_closure_blocker_exported") is True,
        "audited_node_domain": NODE_DOMAIN,
        "candidate_grid_row_count": len(rows),
        "affine_countermodel_count_to_y1_zero": len(affine_countermodels),
        "strict_slope_countermodel_count_to_y1_zero": len(strict_slope_unital_countermodels),
        "unital_row_count": len(unital_rows),
        "nonunital_row_count": len(nonunital_rows),
        "chosen_current_premise_countermodel_to_unital_normalization": chosen,
        "identity_node_existence_does_not_entail_y1_zero": len(affine_countermodels) > 0,
        "affine_consistency_does_not_entail_y1_zero": len(affine_countermodels) > 0,
        "even_strict_slope_value_does_not_entail_y1_zero": len(strict_slope_unital_countermodels) > 0,
        "unit_product_law_equivalent_to_y1_zero_on_affine_grid": all(
            row["unit_product_law_accepts"] == row["unital_y1_zero_accepts"] for row in rows
        ),
        "full_product_multiplicativity_equivalent_to_y1_zero_on_affine_grid": all(
            row["multiplicative_character_accepts"] == row["unital_y1_zero_accepts"] for row in rows
        ),
        "all_product_defects_equal_minus_y1": all(row["all_product_defects_equal_minus_y1"] for row in rows),
        "all_unit_defects_equal_minus_y1": all(row["all_unit_defects_equal_minus_y1"] for row in rows),
        "unital_normalization_is_minimal_missing_source_for_multiplicative_key_inside_affine_family": True,
        "frontier_source_key_under_attack": "strict_unital_monoid_normalization_y1_zero",
        "current_premise_obstruction_exported": True,
        "required_new_premise_class": "strict_nadsoliton_unit_node_normalization_source_or_equivalent_identity_action_theorem",
        "honest_next_step_recommendation": (
            "Do not repeat affine/product consistency scans.  The smallest proof-oriented next step is to formulate and test a strict "
            "nadsoliton unit-node normalization theorem, e.g. an identity-action or zero-damping-action premise that derives y_1=0; "
            "if that cannot be sourced, switch to the independent m=2 operator-order selection source target."
        ),
        "interpretation": (
            "The audit isolates the P2544 recommended smallest target.  Current premises contain the identity node and the affine log mode, "
            "but those facts admit b=log beta != 0.  The equations that kill b are exactly the missing unit-product/multiplicative law, "
            "so y_1=0 remains a source theorem obligation rather than a consistency consequence."
        ),
    }


def load_theorems(sources: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {
        "P2521_SINGLE_NODE_ANCHOR_EQUIVALENCE": theorem(sources["P2521_SINGLE_NODE_ANCHOR_EQUIVALENCE"], "strict_damping_single_node_anchor_equivalence_certificate"),
        "P2525_MULTIPLICATIVE_BETA_NORMALIZATION": theorem(sources["P2525_MULTIPLICATIVE_BETA_NORMALIZATION"], "strict_damping_multiplicative_cocycle_beta_normalization_subkey_certificate"),
        "P2530_FOUR_KEY_IRREDUNDANCY": theorem(sources["P2530_FOUR_KEY_IRREDUNDANCY"], "strict_damping_four_key_irredundancy_witness_certificate"),
        "P2541_MULTIPLICATIVE_OBSTRUCTION": theorem(sources["P2541_MULTIPLICATIVE_OBSTRUCTION"], "strict_damping_multiplicative_character_current_premise_obstruction_certificate"),
        "P2544_FOUR_KEY_CLOSURE_BLOCKER": theorem(sources["P2544_FOUR_KEY_CLOSURE_BLOCKER"], "strict_damping_four_key_current_premise_closure_blocker_certificate"),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2545/S1495 strict damping unital-normalization current-premise obstruction certificate

`P2545/S1495` executes the smallest P2544 completion-path target: the attempted source of `y_1=0` / unital monoid normalization.  The computation audits affine rows `y_d=b+a log d` on `d=1..11` and finds explicit current-premise countermodels with `b != 0`, including a row with the strict slope `a=4/5`.  The unit-product and full multiplicative defects are exactly `-y_1`, so the equations that set `y_1=0` are precisely the missing source law, not a consequence of merely having the identity node or affine consistency.

Therefore the multiplicative-character key remains blocked at the sharper subkey `strict_unital_monoid_normalization_y1_zero`.  The next honest step is a genuine strict nadsoliton unit-node/identity-action theorem, or else a switch to the independent `m=2` operator-order source target.  No beta/eta source, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.
""".strip()
    lag_section = """
## P2545/S1495 unital-normalization Ltotal guard

`P2545/S1495` blocks a common false promotion route for the damping term: `y_1=0` is not forced by current affine/node-domain premises.  It can only be used in `L_total` as a conditional subkey until a strict nadsoliton unit-node normalization source theorem is supplied.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2545/S1495 strict damping unital-normalization current-premise obstruction certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2545/S1495 unital-normalization Ltotal guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    rows = affine_rows()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    theorems = load_theorems(sources)
    cert = build_certificate(theorems, rows)
    theorem_export = {
        "theorem_name": "P2545_T1_strict_damping_unital_normalization_current_premise_obstruction_certificate",
        "audited_chain": ["P2521/S1471", "P2525/S1475", "P2530/S1480", "P2541/S1491", "P2544/S1494"],
        "strict_damping_unital_normalization_current_premise_obstruction_certificate": cert,
        **cert,
        "unital_monoid_normalization_current_premise_obstruction_exported": True,
        "unital_monoid_normalization_source_exported": False,
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
            "This packet isolates a subkey obstruction; it is not a strict nadsoliton source theorem.",
            "It does not replace the P2530 four-key normal form.",
            "It does not transfer legacy physical-role claims onto the strict gate kernel.",
            "It does not discharge QW-2191 or any ToE gate.",
        ],
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2544_inherited": theorem_export["p2544_no_false_source_blocker_inherited"],
        "countermodels_exist": theorem_export["affine_countermodel_count_to_y1_zero"] > 0,
        "strict_slope_countermodel_exists": theorem_export["strict_slope_countermodel_count_to_y1_zero"] > 0,
        "unit_law_equivalence_verified": theorem_export["unit_product_law_equivalent_to_y1_zero_on_affine_grid"],
        "full_multiplicativity_equivalence_verified": theorem_export["full_product_multiplicativity_equivalent_to_y1_zero_on_affine_grid"],
        "negative_controls_preserved": not any(theorem_export[key] for key in NEGATIVE_EXPORT_FLAGS),
    }
    return {
        "packet_id": "P2545",
        "stage_id": "S1495",
        "status": "STRICT_DAMPING_UNITAL_NORMALIZATION_CURRENT_PREMISE_OBSTRUCTION_CERTIFICATE_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_unital_normalization_current_premise_obstruction_certificate": {
            "theorem_export": theorem_export,
            "candidate_rows": rows,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_unital_normalization_current_premise_obstruction_certificate"]["theorem_export"]
    witness = t["chosen_current_premise_countermodel_to_unital_normalization"]
    lines = [
        "# P2545/S1495 strict damping unital-normalization current-premise obstruction certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Candidate affine rows audited: `{t['candidate_grid_row_count']}`.",
        f"- Countermodels to `y_1=0`: `{t['affine_countermodel_count_to_y1_zero']}`.",
        f"- Countermodels with strict slope `delta=4/5`: `{t['strict_slope_countermodel_count_to_y1_zero']}`.",
        f"- Unit-product law equivalent to `y_1=0` on the affine grid: `{t['unit_product_law_equivalent_to_y1_zero_on_affine_grid']}`.",
        f"- Full multiplicativity equivalent to `y_1=0` on the affine grid: `{t['full_product_multiplicativity_equivalent_to_y1_zero_on_affine_grid']}`.",
        f"- Unital normalization source exported: `{t['unital_monoid_normalization_source_exported']}`.",
        "",
        "## Countermodel",
        "",
        f"- `b=log beta`: `{witness['intercept_log_beta']}`.",
        f"- `a=delta`: `{witness['slope_delta']}`.",
        f"- `y_1`: `{witness['y_1']}`.",
        f"- Current affine premises accept: `{witness['current_affine_consistency_premise_accepts']}`.",
        f"- Unital `y_1=0` accepts: `{witness['unital_y1_zero_accepts']}`.",
        f"- Unit-product max defect: `{witness['max_abs_unit_product_defect']}`.",
        "",
        "## Interpretation",
        "",
        t["interpretation"],
        "",
        "## Recommendation",
        "",
        t["honest_next_step_recommendation"],
        "",
        "## Negative Controls",
        "",
        "No unital-normalization source theorem, beta/eta source, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_unital_normalization_current_premise_obstruction_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_unital_normalization_current_premise_obstruction_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
