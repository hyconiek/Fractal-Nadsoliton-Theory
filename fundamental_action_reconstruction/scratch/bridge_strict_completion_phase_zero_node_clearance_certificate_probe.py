#!/usr/bin/env python3
"""Scratch probe: rational node-clearance certificate for phase zeros.

The rational phase-zero certificate located the strict and legacy phase zeros by
interval/parity.  The margin certificate then measured robustness in omega/phi.
This probe records a different finite proof object: every integer node d=0..11
has a positive rational clearance from every relevant phase zero.

Thus the phase sign pattern is not only obtained by odd zero parity on edges;
it is also protected against accidental zero-on-node degeneracy by explicit
node-clearance lower bounds.  This remains a placement/separation certificate,
not a derivation of strict phase transport from nadsoliton dynamics.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_zero_node_clearance_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_zero_node_clearance_certificate_report.md"
RATIONAL_ZERO_REPORT = HERE / "bridge_strict_completion_phase_zero_rational_interval_certificate_report.json"
MARGIN_REPORT = HERE / "bridge_strict_completion_phase_zero_margin_certificate_report.json"
DAMPING_REPORT = HERE / "bridge_strict_completion_damping_continuous_monotonicity_certificate_report.json"

PI_LOWER = Fraction(333, 106)
PI_UPPER = Fraction(355, 113)
STRICT_OMEGA = Fraction(743, 4000)
STRICT_PHI = Fraction(13, 80)
DOMAIN = list(range(12))
LEGACY_ZEROS = [Fraction(4, 3), Fraction(16, 3), Fraction(28, 3)]
STRICT_ZERO_KS = [-1, 0, 1]
EXPECTED_SIGN_PATTERN = [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1]
EXPECTED_FLIP_EDGES = ["1->2", "5->6", "7->8", "9->10"]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def fraction_payload(value: Fraction) -> dict[str, Any]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
        "text": f"{value.numerator}/{value.denominator}",
    }


def strict_zero_interval(k: int) -> tuple[Fraction, Fraction]:
    """Certified interval for x_k=((k+1/2)*pi-phi)/omega."""
    coeff = Fraction(2 * k + 1, 2)
    if coeff >= 0:
        lower_num = coeff * PI_LOWER - STRICT_PHI
        upper_num = coeff * PI_UPPER - STRICT_PHI
    else:
        lower_num = coeff * PI_UPPER - STRICT_PHI
        upper_num = coeff * PI_LOWER - STRICT_PHI
    return lower_num / STRICT_OMEGA, upper_num / STRICT_OMEGA


def interval_distance_lower_bound(point: Fraction, lower: Fraction, upper: Fraction) -> Fraction:
    if point < lower:
        return lower - point
    if point > upper:
        return point - upper
    return Fraction(0, 1)


def edge_label_for_point_interval(lower: Fraction, upper: Fraction) -> str:
    for d in range(11):
        if Fraction(d, 1) < lower and upper < Fraction(d + 1, 1):
            return f"{d}->{d + 1}"
    if upper < 0:
        return "left-of-domain"
    if lower > 11:
        return "right-of-domain"
    return "not-contained-in-single-open-integer-edge"


def legacy_zero_rows() -> list[dict[str, Any]]:
    rows = []
    for idx, zero in enumerate(LEGACY_ZEROS):
        node_distances = [abs(Fraction(d, 1) - zero) for d in DOMAIN]
        min_distance = min(node_distances)
        rows.append({
            "k": idx,
            "zero": fraction_payload(zero),
            "edge_containing_zero": edge_label_for_point_interval(zero, zero),
            "min_integer_node_clearance": fraction_payload(min_distance),
            "nearest_integer_nodes": [d for d, dist in zip(DOMAIN, node_distances) if dist == min_distance],
            "zero_is_not_integer_node": min_distance > 0,
        })
    return rows


def strict_zero_rows() -> list[dict[str, Any]]:
    rows = []
    for k in STRICT_ZERO_KS:
        lower, upper = strict_zero_interval(k)
        node_clearances = [interval_distance_lower_bound(Fraction(d, 1), lower, upper) for d in DOMAIN]
        min_clearance = min(node_clearances)
        rows.append({
            "k": k,
            "x_interval": {"lower": fraction_payload(lower), "upper": fraction_payload(upper)},
            "edge_or_domain_location": edge_label_for_point_interval(lower, upper),
            "min_integer_node_clearance_lower_bound": fraction_payload(min_clearance),
            "nearest_integer_nodes_by_bound": [d for d, dist in zip(DOMAIN, node_clearances) if dist == min_clearance],
            "all_integer_nodes_certified_away_from_zero_interval": min_clearance > 0,
        })
    return rows


def integer_node_rows(strict_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    strict_intervals = [
        (row["k"], Fraction(row["x_interval"]["lower"]["numerator"], row["x_interval"]["lower"]["denominator"]), Fraction(row["x_interval"]["upper"]["numerator"], row["x_interval"]["upper"]["denominator"]))
        for row in strict_rows
    ]
    rows = []
    for d in DOMAIN:
        point = Fraction(d, 1)
        legacy_distances = [abs(point - zero) for zero in LEGACY_ZEROS]
        strict_distances = [interval_distance_lower_bound(point, lower, upper) for _, lower, upper in strict_intervals]
        combined = min(legacy_distances + strict_distances)
        rows.append({
            "d": d,
            "legacy_nearest_clearance": fraction_payload(min(legacy_distances)),
            "strict_nearest_clearance_lower_bound": fraction_payload(min(strict_distances)),
            "combined_phase_zero_clearance_lower_bound": fraction_payload(combined),
            "node_certified_not_phase_zero": combined > 0,
        })
    return rows


def build_payload() -> dict[str, Any]:
    rational_zero = load_json(RATIONAL_ZERO_REPORT)
    margin = load_json(MARGIN_REPORT)
    damping = load_json(DAMPING_REPORT)
    legacy_rows = legacy_zero_rows()
    strict_rows = strict_zero_rows()
    node_rows = integer_node_rows(strict_rows)

    min_legacy_clearance = min(Fraction(row["min_integer_node_clearance"]["numerator"], row["min_integer_node_clearance"]["denominator"]) for row in legacy_rows)
    min_strict_clearance = min(Fraction(row["min_integer_node_clearance_lower_bound"]["numerator"], row["min_integer_node_clearance_lower_bound"]["denominator"]) for row in strict_rows)
    min_combined_node_clearance = min(Fraction(row["combined_phase_zero_clearance_lower_bound"]["numerator"], row["combined_phase_zero_clearance_lower_bound"]["denominator"]) for row in node_rows)

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_ZERO_NODE_CLEARANCE_CERTIFICATE__RATIONAL_NO_NODE_DEGENERACY",
        "status": "all-integer-phase-nodes-have-positive-rational-clearance-from-legacy-and-strict-zeros",
        "source_reports": {
            "rational_zero_certificate": str(RATIONAL_ZERO_REPORT.relative_to(ROOT)),
            "phase_zero_margin_certificate": str(MARGIN_REPORT.relative_to(ROOT)),
            "damping_monotonicity_certificate": str(DAMPING_REPORT.relative_to(ROOT)),
        },
        "rational_inputs": {
            "pi_lower_bound": fraction_payload(PI_LOWER),
            "pi_upper_bound": fraction_payload(PI_UPPER),
            "strict_omega": fraction_payload(STRICT_OMEGA),
            "strict_phi": fraction_payload(STRICT_PHI),
            "legacy_zero_positions_exact": [fraction_payload(z) for z in LEGACY_ZEROS],
            "audited_integer_domain": DOMAIN,
        },
        "legacy_zero_clearance_rows": legacy_rows,
        "strict_zero_clearance_rows": strict_rows,
        "integer_node_clearance_rows": node_rows,
        "clearance_summary": {
            "all_legacy_zeros_off_integer_nodes": all(row["zero_is_not_integer_node"] for row in legacy_rows),
            "all_strict_zero_intervals_off_integer_nodes": all(row["all_integer_nodes_certified_away_from_zero_interval"] for row in strict_rows),
            "all_integer_nodes_certified_not_phase_zeros": all(row["node_certified_not_phase_zero"] for row in node_rows),
            "min_legacy_integer_node_clearance": fraction_payload(min_legacy_clearance),
            "min_strict_integer_node_clearance_lower_bound": fraction_payload(min_strict_clearance),
            "min_combined_phase_zero_node_clearance_lower_bound": fraction_payload(min_combined_node_clearance),
            "certified_phase_sign_flip_edges_preserved": rational_zero["interval_summary"]["phase_sign_flip_edges_from_rational_intervals"],
            "certified_phase_sign_pattern_preserved": rational_zero["interval_summary"]["phase_transport_sign_pattern_from_rational_intervals"],
            "matches_expected_phase_flips": rational_zero["interval_summary"]["phase_sign_flip_edges_from_rational_intervals"] == EXPECTED_FLIP_EDGES,
            "matches_expected_sign_pattern": rational_zero["interval_summary"]["phase_transport_sign_pattern_from_rational_intervals"] == EXPECTED_SIGN_PATTERN,
            "damping_positive_so_node_clearance_controls_sign_degeneracy": damping["monotonicity_summary"]["continuous_positive_certificate"],
            "parameter_margin_report_still_passes": margin["robustness_summary"]["all_worst_case_inequalities_hold_at_epsilon"],
        },
        "blocker_context": {
            "what_this_refines": "adds rational no-node-degeneracy clearances to the phase-zero parity and parameter-margin certificates",
            "rational_zero_status": rational_zero["status"],
            "margin_status": margin["status"],
            "damping_status": damping["status"],
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_transport_derivation_from_nadsoliton_dynamics",
                "orientation_chi11_source",
                "chi11_uniqueness",
                "role_transfer_theorem",
            ],
        },
        "proof_certificate": {
            "legacy_clearance_step": "Legacy zeros are exact rationals 4/3, 16/3, 28/3, hence their minimum integer-node clearance is exactly positive.",
            "strict_clearance_step": "Each relevant strict zero is enclosed using 333/106 < pi < 355/113 and omega=743/4000, phi=13/80; integer-node distance to the whole interval gives a rational lower bound.",
            "node_step": "For each d=0..11, the minimum of all legacy and strict zero clearances is positive, so no audited integer node is a phase zero.",
            "sign_step": "With damping already positive, sign changes in the completion factor remain controlled by phase-zero parity rather than damping or node degeneracy.",
            "nonduplication": "This is a node-clearance/no-degeneracy certificate, not another parameter-robustness, damping-monotonicity, transport-cocycle, or chain-integrity audit.",
            "theoretical_limit": "The certificate proves finite rational separation of selected phase zeros from integer nodes; it does not derive omega/phi or transport from strict nadsoliton dynamics.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in solitonic state; this audit only checks finite phase-zero separation.",
            "forbidden_reading": "No separate informational layer below the nadsoliton is introduced.",
        },
        "hard_limits": [
            "K_strict_gate remains the current live/full operational kernel.",
            "No unqualified identity K_legacy_ont == K_strict_gate is claimed.",
            "No proof derives strict omega/phi or phase transport from strict nadsoliton dynamics.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No legacy physical-role transfer to K_strict_gate is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    summary = payload["clearance_summary"]
    lines = [
        "# Strict completion phase-zero node-clearance certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This audit gives rational lower bounds proving that no audited integer",
        "node `d=0..11` is accidentally sitting on a legacy or strict phase zero.",
        "It complements edge-parity and parameter-margin certificates without",
        "claiming a strict dynamical derivation of the phase parameters.",
        "",
        "## Clearance summary",
        "",
        f"- All legacy zeros off integer nodes: `{summary['all_legacy_zeros_off_integer_nodes']}`",
        f"- All strict zero intervals off integer nodes: `{summary['all_strict_zero_intervals_off_integer_nodes']}`",
        f"- All audited integer nodes certified nonzero: `{summary['all_integer_nodes_certified_not_phase_zeros']}`",
        f"- Minimum legacy clearance: `{summary['min_legacy_integer_node_clearance']['text']}` = `{summary['min_legacy_integer_node_clearance']['decimal']:.12e}`",
        f"- Minimum strict clearance lower bound: `{summary['min_strict_integer_node_clearance_lower_bound']['text']}` = `{summary['min_strict_integer_node_clearance_lower_bound']['decimal']:.12e}`",
        f"- Minimum combined node clearance lower bound: `{summary['min_combined_phase_zero_node_clearance_lower_bound']['text']}` = `{summary['min_combined_phase_zero_node_clearance_lower_bound']['decimal']:.12e}`",
        f"- Preserved flip edges: `{summary['certified_phase_sign_flip_edges_preserved']}`",
        f"- Preserved sign pattern: `{summary['certified_phase_sign_pattern_preserved']}`",
        "",
        "## Strict zero clearance rows",
        "",
        "| k | location | lower | upper | min node clearance lower bound | nearest nodes |",
        "|---:|---|---:|---:|---:|---|",
    ]
    for row in payload["strict_zero_clearance_rows"]:
        lines.append(
            "| {k} | {loc} | {lo} | {hi} | {clearance} | {nodes} |".format(
                k=row["k"],
                loc=row["edge_or_domain_location"],
                lo=row["x_interval"]["lower"]["text"],
                hi=row["x_interval"]["upper"]["text"],
                clearance=row["min_integer_node_clearance_lower_bound"]["text"],
                nodes=row["nearest_integer_nodes_by_bound"],
            )
        )
    lines.extend([
        "",
        "## Integer node rows",
        "",
        "| d | legacy clearance | strict clearance lower bound | combined lower bound | certified nonzero |",
        "|---:|---:|---:|---:|---:|",
    ])
    for row in payload["integer_node_clearance_rows"]:
        lines.append(
            f"| {row['d']} | {row['legacy_nearest_clearance']['text']} | {row['strict_nearest_clearance_lower_bound']['text']} | {row['combined_phase_zero_clearance_lower_bound']['text']} | {row['node_certified_not_phase_zero']} |"
        )
    lines.extend([
        "",
        "## Proof certificate",
        "",
    ])
    for key, value in payload["proof_certificate"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend([
        "",
        "## Hard limits",
        "",
    ])
    for limit in payload["hard_limits"]:
        lines.append(f"- {limit}")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
