#!/usr/bin/env python3
"""Scratch probe: rational cell-partition certificate for phase zeros.

The node-clearance certificate proved that no audited integer node d=0..11 is
itself a phase zero.  This probe records the next finite separation object: the
phase-zero carriers inside the audited domain are ordered, pairwise disjoint,
separated from integer edge boundaries, and induce the same edge parity/sign
pattern as the earlier rational-zero certificates.

It is still a placement/cell-partition certificate for the selected phase
parameters, not a strict dynamical derivation of those parameters.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_zero_cell_partition_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_zero_cell_partition_certificate_report.md"
RATIONAL_ZERO_REPORT = HERE / "bridge_strict_completion_phase_zero_rational_interval_certificate_report.json"
NODE_CLEARANCE_REPORT = HERE / "bridge_strict_completion_phase_zero_node_clearance_certificate_report.json"
MARGIN_REPORT = HERE / "bridge_strict_completion_phase_zero_margin_certificate_report.json"

PI_LOWER = Fraction(333, 106)
PI_UPPER = Fraction(355, 113)
STRICT_OMEGA = Fraction(743, 4000)
STRICT_PHI = Fraction(13, 80)
DOMAIN_LEFT = Fraction(0, 1)
DOMAIN_RIGHT = Fraction(11, 1)
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
    coeff = Fraction(2 * k + 1, 2)
    if coeff >= 0:
        lower_num = coeff * PI_LOWER - STRICT_PHI
        upper_num = coeff * PI_UPPER - STRICT_PHI
    else:
        lower_num = coeff * PI_UPPER - STRICT_PHI
        upper_num = coeff * PI_LOWER - STRICT_PHI
    return lower_num / STRICT_OMEGA, upper_num / STRICT_OMEGA


def edge_label_for_interval(lower: Fraction, upper: Fraction) -> str:
    for d in range(int(DOMAIN_RIGHT)):
        if Fraction(d, 1) < lower and upper < Fraction(d + 1, 1):
            return f"{d}->{d + 1}"
    if upper < DOMAIN_LEFT:
        return "left-of-domain"
    if lower > DOMAIN_RIGHT:
        return "right-of-domain"
    return "not-contained-in-single-open-integer-edge"


def carrier_payload(label: str, source: str, k: int, lower: Fraction, upper: Fraction) -> dict[str, Any]:
    edge = edge_label_for_interval(lower, upper)
    in_domain = DOMAIN_LEFT < lower and upper < DOMAIN_RIGHT
    return {
        "label": label,
        "source": source,
        "k": k,
        "lower": fraction_payload(lower),
        "upper": fraction_payload(upper),
        "width": fraction_payload(upper - lower),
        "edge_or_domain_location": edge,
        "strictly_inside_audit_domain": in_domain,
        "strictly_inside_single_open_integer_edge": "->" in edge,
    }


def all_carriers() -> list[dict[str, Any]]:
    carriers: list[dict[str, Any]] = []
    for k, zero in enumerate(LEGACY_ZEROS):
        carriers.append(carrier_payload(f"legacy_z{k}", "legacy_exact", k, zero, zero))
    for k in STRICT_ZERO_KS:
        lower, upper = strict_zero_interval(k)
        carriers.append(carrier_payload(f"strict_z{k}", "strict_interval", k, lower, upper))
    return carriers


def as_interval(carrier: dict[str, Any]) -> tuple[Fraction, Fraction]:
    lower = carrier["lower"]
    upper = carrier["upper"]
    return Fraction(lower["numerator"], lower["denominator"]), Fraction(upper["numerator"], upper["denominator"])


def domain_carriers(carriers: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(
        [carrier for carrier in carriers if carrier["strictly_inside_audit_domain"]],
        key=lambda carrier: as_interval(carrier)[0],
    )


def boundary_clearance_rows(carriers: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for carrier in carriers:
        lower, upper = as_interval(carrier)
        edge_left = Fraction(int(lower), 1)
        if edge_left == lower:
            edge_left -= 1
        edge_right = edge_left + 1
        left_clearance = lower - edge_left
        right_clearance = edge_right - upper
        rows.append({
            "label": carrier["label"],
            "edge": carrier["edge_or_domain_location"],
            "left_boundary_clearance": fraction_payload(left_clearance),
            "right_boundary_clearance": fraction_payload(right_clearance),
            "min_boundary_clearance": fraction_payload(min(left_clearance, right_clearance)),
            "strictly_inside_open_edge": left_clearance > 0 and right_clearance > 0,
        })
    return rows


def adjacent_separation_rows(carriers: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for left, right in zip(carriers, carriers[1:]):
        _, left_upper = as_interval(left)
        right_lower, _ = as_interval(right)
        separation = right_lower - left_upper
        rows.append({
            "left_label": left["label"],
            "right_label": right["label"],
            "separation": fraction_payload(separation),
            "strictly_disjoint_and_ordered": separation > 0,
        })
    return rows


def edge_inventory_rows(carriers: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for d in range(int(DOMAIN_RIGHT)):
        edge = f"{d}->{d + 1}"
        labels = [carrier["label"] for carrier in carriers if carrier["edge_or_domain_location"] == edge]
        rows.append({
            "edge": edge,
            "zero_carrier_labels": labels,
            "zero_carrier_count": len(labels),
            "odd_parity_flip": len(labels) % 2 == 1,
        })
    return rows


def sign_pattern_from_edge_inventory(rows: list[dict[str, Any]]) -> list[int]:
    signs = [1]
    current = 1
    for row in rows:
        if row["odd_parity_flip"]:
            current *= -1
        signs.append(current)
    return signs


def cell_partition_rows(carriers: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    current_left = DOMAIN_LEFT
    for carrier in carriers:
        lower, upper = as_interval(carrier)
        rows.append({
            "cell_left": fraction_payload(current_left),
            "cell_right": fraction_payload(lower),
            "right_zero_carrier": carrier["label"],
            "positive_length": lower - current_left > 0,
        })
        current_left = upper
    rows.append({
        "cell_left": fraction_payload(current_left),
        "cell_right": fraction_payload(DOMAIN_RIGHT),
        "right_zero_carrier": "domain_right_boundary",
        "positive_length": DOMAIN_RIGHT - current_left > 0,
    })
    return rows


def build_payload() -> dict[str, Any]:
    rational_zero = load_json(RATIONAL_ZERO_REPORT)
    node_clearance = load_json(NODE_CLEARANCE_REPORT)
    margin = load_json(MARGIN_REPORT)
    carriers = all_carriers()
    in_domain = domain_carriers(carriers)
    boundary_rows = boundary_clearance_rows(in_domain)
    separation_rows = adjacent_separation_rows(in_domain)
    inventory_rows = edge_inventory_rows(in_domain)
    partition_rows = cell_partition_rows(in_domain)
    derived_flip_edges = [row["edge"] for row in inventory_rows if row["odd_parity_flip"]]
    derived_sign_pattern = sign_pattern_from_edge_inventory(inventory_rows)
    min_boundary = min(Fraction(row["min_boundary_clearance"]["numerator"], row["min_boundary_clearance"]["denominator"]) for row in boundary_rows)
    min_separation = min(Fraction(row["separation"]["numerator"], row["separation"]["denominator"]) for row in separation_rows)
    min_cell_length = min(Fraction(row["cell_right"]["numerator"], row["cell_right"]["denominator"]) - Fraction(row["cell_left"]["numerator"], row["cell_left"]["denominator"]) for row in partition_rows)

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_ZERO_CELL_PARTITION_CERTIFICATE__RATIONAL_ORDERED_ZERO_CARRIERS",
        "status": "phase-zero-carriers-form-positive-rational-cell-partition-on-audited-domain",
        "source_reports": {
            "rational_zero_certificate": str(RATIONAL_ZERO_REPORT.relative_to(ROOT)),
            "node_clearance_certificate": str(NODE_CLEARANCE_REPORT.relative_to(ROOT)),
            "phase_zero_margin_certificate": str(MARGIN_REPORT.relative_to(ROOT)),
        },
        "rational_inputs": {
            "pi_lower_bound": fraction_payload(PI_LOWER),
            "pi_upper_bound": fraction_payload(PI_UPPER),
            "strict_omega": fraction_payload(STRICT_OMEGA),
            "strict_phi": fraction_payload(STRICT_PHI),
            "legacy_zero_positions_exact": [fraction_payload(zero) for zero in LEGACY_ZEROS],
            "audited_domain": [int(DOMAIN_LEFT), int(DOMAIN_RIGHT)],
        },
        "all_zero_carriers": carriers,
        "domain_zero_carriers_ordered": in_domain,
        "domain_zero_boundary_clearance_rows": boundary_rows,
        "adjacent_domain_zero_separation_rows": separation_rows,
        "edge_zero_inventory_rows": inventory_rows,
        "cell_partition_rows": partition_rows,
        "cell_partition_summary": {
            "domain_zero_carrier_order": [carrier["label"] for carrier in in_domain],
            "all_domain_zero_carriers_inside_open_edges": all(row["strictly_inside_open_edge"] for row in boundary_rows),
            "all_adjacent_domain_zero_carriers_strictly_ordered_and_disjoint": all(row["strictly_disjoint_and_ordered"] for row in separation_rows),
            "all_cells_have_positive_rational_length": all(row["positive_length"] for row in partition_rows),
            "min_boundary_clearance": fraction_payload(min_boundary),
            "min_adjacent_zero_separation": fraction_payload(min_separation),
            "min_cell_length": fraction_payload(min_cell_length),
            "derived_phase_sign_flip_edges": derived_flip_edges,
            "derived_phase_transport_sign_pattern": derived_sign_pattern,
            "matches_rational_zero_flip_edges": derived_flip_edges == rational_zero["interval_summary"]["phase_sign_flip_edges_from_rational_intervals"],
            "matches_rational_zero_sign_pattern": derived_sign_pattern == rational_zero["interval_summary"]["phase_transport_sign_pattern_from_rational_intervals"],
            "matches_node_clearance_pattern": derived_sign_pattern == node_clearance["clearance_summary"]["certified_phase_sign_pattern_preserved"],
            "matches_expected_phase_flips": derived_flip_edges == EXPECTED_FLIP_EDGES,
            "matches_expected_sign_pattern": derived_sign_pattern == EXPECTED_SIGN_PATTERN,
            "margin_report_still_passes": margin["robustness_summary"]["all_worst_case_inequalities_hold_at_epsilon"],
            "node_clearance_report_still_passes": node_clearance["clearance_summary"]["all_integer_nodes_certified_not_phase_zeros"],
        },
        "blocker_context": {
            "what_this_refines": "adds an ordered rational cell partition for phase-zero carriers after node-clearance and margin certificates",
            "rational_zero_status": rational_zero["status"],
            "node_clearance_status": node_clearance["status"],
            "margin_status": margin["status"],
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_transport_derivation_from_nadsoliton_dynamics",
                "orientation_chi11_source",
                "chi11_uniqueness",
                "role_transfer_theorem",
            ],
        },
        "proof_certificate": {
            "carrier_step": "Legacy exact zeros and strict interval zeros are represented as rational zero-carriers on the audited line.",
            "edge_step": "Every in-domain zero-carrier is strictly inside a single open integer edge with positive rational boundary clearance.",
            "separation_step": "Adjacent in-domain zero-carriers have positive rational separation, so no two phase-zero events collide or overlap.",
            "cell_step": "The ordered carriers cut [0,11] into positive-length rational cells; signs can change only across odd-parity carrier edges.",
            "parity_step": "The induced odd-parity edge inventory recovers the established flip edges and sign pattern.",
            "nonduplication": "This is a cell-partition/order/separation certificate, not another node-clearance, margin, damping, cocycle, or chain-integrity audit.",
            "theoretical_limit": "The certificate proves finite rational phase-zero carrier separation; it does not derive omega/phi or transport from strict nadsoliton dynamics.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in solitonic state; this audit only checks finite phase-zero carrier separation.",
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
    summary = payload["cell_partition_summary"]
    lines = [
        "# Strict completion phase-zero cell-partition certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This audit proves that the in-domain phase-zero carriers are ordered,",
        "pairwise separated, and strictly inside open integer edges.  The resulting",
        "cell partition recovers the established phase flip edges and sign pattern,",
        "without claiming a strict dynamical derivation of the phase parameters.",
        "",
        "## Cell-partition summary",
        "",
        f"- Domain zero carrier order: `{summary['domain_zero_carrier_order']}`",
        f"- All carriers inside open edges: `{summary['all_domain_zero_carriers_inside_open_edges']}`",
        f"- Adjacent carriers ordered/disjoint: `{summary['all_adjacent_domain_zero_carriers_strictly_ordered_and_disjoint']}`",
        f"- All cells positive length: `{summary['all_cells_have_positive_rational_length']}`",
        f"- Minimum boundary clearance: `{summary['min_boundary_clearance']['text']}` = `{summary['min_boundary_clearance']['decimal']:.12e}`",
        f"- Minimum adjacent zero separation: `{summary['min_adjacent_zero_separation']['text']}` = `{summary['min_adjacent_zero_separation']['decimal']:.12e}`",
        f"- Minimum cell length: `{summary['min_cell_length']['text']}` = `{summary['min_cell_length']['decimal']:.12e}`",
        f"- Derived flip edges: `{summary['derived_phase_sign_flip_edges']}`",
        f"- Derived sign pattern: `{summary['derived_phase_transport_sign_pattern']}`",
        "",
        "## Edge zero inventory",
        "",
        "| edge | zero carriers | count | odd parity flip |",
        "|---|---|---:|---:|",
    ]
    for row in payload["edge_zero_inventory_rows"]:
        lines.append(f"| {row['edge']} | {row['zero_carrier_labels']} | {row['zero_carrier_count']} | {row['odd_parity_flip']} |")
    lines.extend([
        "",
        "## Adjacent zero separation",
        "",
        "| left | right | separation | ordered/disjoint |",
        "|---|---|---:|---:|",
    ])
    for row in payload["adjacent_domain_zero_separation_rows"]:
        lines.append(f"| {row['left_label']} | {row['right_label']} | {row['separation']['text']} | {row['strictly_disjoint_and_ordered']} |")
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
