#!/usr/bin/env python3
"""Scratch probe: D12-invariant selector no-go for the anchored d5 orbit.

The previous d5 selector symmetry audit found that the fully anchored
configuration is not a canonical object: the 24 source/orientation choices form
one free dihedral D_12 orbit.  This probe turns that observation into a sharper,
more theorem-like obstruction.

If a score S is D_12-invariant, then S(g.x)=S(x) for every dihedral action g.
On a transitive orbit that makes S constant on all 24 anchored representatives.
Therefore a D_12-invariant selector cannot choose a unique source/orientation
representative from the anchored d5 orbit.  Any unique selector must introduce
non-invariant data (a source, orientation, external phase, boundary, or another
strictly derived symmetry-breaking term).  The computation below verifies the
orbit replay and audits several invariant feature packets that all collapse to a
single value over the 24 representatives, while fixed-frame non-invariant moments
vary and thus only select after an anchor is supplied.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SYMMETRY_AUDIT = HERE / "bridge_strict_alpha_d5_selector_symmetry_orbit_audit_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_d12_invariant_selector_no_go_report.json"
OUT_MD = HERE / "bridge_strict_alpha_d12_invariant_selector_no_go_report.md"

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
DISTANCE_SELECTED = 5
TARGET_BINARY_EXPONENT = 8
DENOMINATOR = 3
STRICT_TARGET_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
FORWARD_ASSIGNMENT = (2, 2, 2, 1, 1)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def eta_from_product(product: Fraction, branch_count: int) -> float:
    correction = float(product) ** (1.0 / branch_count)
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def ordered_d5_support(source: int, orientation: int) -> tuple[int, ...]:
    if orientation not in {-1, 1}:
        raise ValueError("orientation must be +/-1")
    step = (orientation * DISTANCE_SELECTED) % Z12_NODE_COUNT
    return tuple((source + index * step) % Z12_NODE_COUNT for index in range(SUPPORT_SIZE))


def value_configuration(source: int, orientation: int, assignment: tuple[int, ...] = FORWARD_ASSIGNMENT) -> tuple[int, ...]:
    values = [0] * Z12_NODE_COUNT
    for node, value in zip(ordered_d5_support(source, orientation), assignment):
        values[node] = value
    return tuple(values)


def dihedral_action_on_config(config: tuple[int, ...], shift: int, reflect: bool) -> tuple[int, ...]:
    out = [0] * Z12_NODE_COUNT
    for node, value in enumerate(config):
        target = ((-node if reflect else node) + shift) % Z12_NODE_COUNT
        out[target] = value
    return tuple(out)


def cyclic_distance(left: int, right: int) -> int:
    raw = abs(left - right) % Z12_NODE_COUNT
    return min(raw, Z12_NODE_COUNT - raw)


def support_distance_histogram(config: tuple[int, ...]) -> list[int]:
    support = [index for index, value in enumerate(config) if value != 0]
    histogram = [0] * 7
    for left_index, left in enumerate(support):
        for right in support[left_index + 1 :]:
            histogram[cyclic_distance(left, right)] += 1
    return histogram[1:]


def cyclic_autocorrelation(config: tuple[int, ...]) -> list[int]:
    return [sum(config[node] * config[(node + distance) % Z12_NODE_COUNT] for node in range(Z12_NODE_COUNT)) for distance in range(7)]


def unordered_pair_product_by_distance(config: tuple[int, ...]) -> list[int]:
    by_distance = [0] * 7
    for left in range(Z12_NODE_COUNT):
        for right in range(left + 1, Z12_NODE_COUNT):
            distance = cyclic_distance(left, right)
            by_distance[distance] += config[left] * config[right]
    return by_distance[1:]


def value_count_vector(config: tuple[int, ...]) -> list[int]:
    max_value = max(config)
    return [sum(1 for value in config if value == candidate) for candidate in range(max_value + 1)]


def invariant_feature_packet(config: tuple[int, ...]) -> dict[str, Any]:
    return {
        "value_count_vector": value_count_vector(config),
        "sorted_nonzero_values_desc": sorted((value for value in config if value), reverse=True),
        "support_distance_histogram_d1_to_d6": support_distance_histogram(config),
        "cyclic_autocorrelation_d0_to_d6": cyclic_autocorrelation(config),
        "unordered_pair_product_by_distance_d1_to_d6": unordered_pair_product_by_distance(config),
        "total_mass": sum(config),
        "total_square_mass": sum(value * value for value in config),
    }


def fixed_frame_moments(config: tuple[int, ...]) -> dict[str, int]:
    return {
        "linear_node_moment": sum(node * value for node, value in enumerate(config)),
        "quadratic_node_moment": sum(node * node * value for node, value in enumerate(config)),
        "left_half_mass_nodes_0_to_5": sum(config[:6]),
        "right_half_mass_nodes_6_to_11": sum(config[6:]),
    }


def all_anchored_rows() -> list[dict[str, Any]]:
    rows = []
    for source in range(Z12_NODE_COUNT):
        for orientation in (-1, 1):
            config = value_configuration(source, orientation)
            rows.append(
                {
                    "source": source,
                    "orientation": orientation,
                    "ordered_support": list(ordered_d5_support(source, orientation)),
                    "canonical_support": sorted(ordered_d5_support(source, orientation)),
                    "value_configuration": list(config),
                    "invariant_features": invariant_feature_packet(config),
                    "fixed_frame_noninvariant_moments": fixed_frame_moments(config),
                }
            )
    return rows


def unique_json_packets(rows: list[dict[str, Any]], key: str) -> list[str]:
    return sorted({json.dumps(row[key], sort_keys=True) for row in rows})


def invariant_feature_class_counts(rows: list[dict[str, Any]]) -> dict[str, int]:
    feature_names = rows[0]["invariant_features"].keys()
    return {
        name: len({json.dumps(row["invariant_features"][name], sort_keys=True) for row in rows})
        for name in feature_names
    }


def fixed_frame_score_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    moments = rows[0]["fixed_frame_noninvariant_moments"].keys()
    summary: dict[str, Any] = {}
    for name in moments:
        values = [row["fixed_frame_noninvariant_moments"][name] for row in rows]
        minimum = min(values)
        maximum = max(values)
        summary[name] = {
            "unique_value_count": len(set(values)),
            "min": minimum,
            "max": maximum,
            "minimizers": [
                {"source": row["source"], "orientation": row["orientation"], "value": row["fixed_frame_noninvariant_moments"][name]}
                for row in rows
                if row["fixed_frame_noninvariant_moments"][name] == minimum
            ],
            "maximizers": [
                {"source": row["source"], "orientation": row["orientation"], "value": row["fixed_frame_noninvariant_moments"][name]}
                for row in rows
                if row["fixed_frame_noninvariant_moments"][name] == maximum
            ],
        }
    return summary


def orbit_replay() -> dict[str, Any]:
    base = value_configuration(0, -1)
    orbit = {
        dihedral_action_on_config(base, shift, reflect)
        for shift in range(Z12_NODE_COUNT)
        for reflect in (False, True)
    }
    stabilizer = [
        {"shift": shift, "reflect": reflect}
        for shift in range(Z12_NODE_COUNT)
        for reflect in (False, True)
        if dihedral_action_on_config(base, shift, reflect) == base
    ]
    return {
        "base_configuration": list(base),
        "dihedral_group_order": 2 * Z12_NODE_COUNT,
        "computed_orbit_size": len(orbit),
        "computed_stabilizer_size": len(stabilizer),
        "computed_stabilizer": stabilizer,
    }


def main() -> None:
    symmetry_audit = load_json(SYMMETRY_AUDIT)
    product = Fraction(2 ** TARGET_BINARY_EXPONENT, DENOMINATOR ** SUPPORT_SIZE)
    rows = all_anchored_rows()
    orbit = orbit_replay()
    invariant_packets = unique_json_packets(rows, "invariant_features")
    feature_counts = invariant_feature_class_counts(rows)
    fixed_frame_summary = fixed_frame_score_summary(rows)

    report = {
        "status": "OPEN_STRICT_ALPHA_D12_INVARIANT_SELECTOR_NO_GO_NO_STRICT_SELECTOR_DISCHARGE",
        "result_kind": "SCRATCH_STRICT_ALPHA_D12_INVARIANT_SELECTOR_NO_GO_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "d5_selector_symmetry_orbit_audit": str(SYMMETRY_AUDIT.relative_to(ROOT)),
        },
        "previous_symmetry_audit_replay": {
            "result_kind": symmetry_audit["result_kind"],
            "anchored_rows_equal_full_dihedral_orbit": symmetry_audit["symmetry_verdict"]["anchored_rows_equal_full_dihedral_orbit"],
            "value_configuration_has_trivial_stabilizer": symmetry_audit["symmetry_verdict"]["value_configuration_has_trivial_stabilizer"],
            "reported_value_configuration_orbit_size": symmetry_audit["dihedral_orbit_audit"]["value_configuration_orbit_size"],
            "reported_unoriented_support_orbit_size": symmetry_audit["dihedral_orbit_audit"]["support_orbit_size_without_orientation_assignment"],
        },
        "target_identity_replay": {
            "q_power_product": f"{product.numerator}/{product.denominator}",
            "eta_residual_vs_9_5": eta_from_product(product, SUPPORT_SIZE) - STRICT_TARGET_ETA,
            "forward_assignment": list(FORWARD_ASSIGNMENT),
        },
        "d12_orbit_replay": orbit,
        "finite_no_go_statement": {
            "premise": "A selector score S is D12-invariant on the anchored d5 configuration orbit: S(g.x)=S(x).",
            "transitivity_certificate": "The 24 anchored source/orientation configurations are exactly one D12 orbit with trivial stabilizer.",
            "conclusion": "Every D12-invariant score is constant on this orbit and cannot uniquely select a source/orientation representative.",
            "what_can_select": "A unique representative requires non-invariant source/orientation/external-phase data or a newly exported strict symmetry-breaking term.",
        },
        "invariant_feature_audit": {
            "row_count": len(rows),
            "unique_invariant_feature_packet_count": len(invariant_packets),
            "all_checked_invariant_features_constant": len(invariant_packets) == 1 and all(count == 1 for count in feature_counts.values()),
            "feature_class_unique_counts": feature_counts,
            "representative_invariant_feature_packet": json.loads(invariant_packets[0]),
        },
        "noninvariant_contrast": {
            "fixed_frame_score_summary": fixed_frame_summary,
            "interpretation": "Fixed-frame moments vary across the orbit, so they can break the degeneracy only by importing a frame/source bias; they are not D12-invariant selectors.",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": True,
            "content": "The remaining source/orientation ambiguity is a finite group-theoretic no-go for D12-invariant selectors, not merely an observed scan degeneracy.",
            "why_this_is_more_proof_like": "Transitivity plus invariance is enough to prove constancy for every D12-invariant score; the computation verifies the finite orbit and invariant feature collapse.",
            "why_this_is_not_enough": "The no-go identifies the missing premise; it does not derive the required strict symmetry-breaking source term.",
            "status": "candidate-supported-but-D12-invariant-selectors-provably-insufficient",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "No D12-invariant selector can choose a unique source/orientation representative on the free anchored d5 orbit.",
            "Non-invariant fixed-frame moments are diagnostics only; they are not strict-core selector theorems.",
            "No theorem derives the source/orientation/external-phase symmetry breaker from strict nadsoliton geometry.",
            "No theorem derives branch count m=5, total exponent n=8, denominator 3, or binary-rescale quotient.",
            "No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged by this no-go audit.",
            "No theorem derives eta=9/5 without the branch model, denominator/quotient convention, support premise, ledger selector, source/orientation premise, and assignment premise.",
            "No QW-2191 discharge and no ToE closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Search for, or prove absence of, a strict-source symmetry-breaking term that is not D12-invariant on the anchored d5 orbit; without it the source/orientation remains an explicit non-strict premise.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha D12-invariant selector no-go probe\n\n"
        "Status: finite group-theoretic no-go for D12-invariant source/orientation selection; no strict selector discharge.\n\n"
        f"- D12 orbit replay: group order `{orbit['dihedral_group_order']}`, orbit size `{orbit['computed_orbit_size']}`, stabilizer size `{orbit['computed_stabilizer_size']}`.\n"
        f"- Invariant feature packets over the 24 anchored representatives: `{len(invariant_packets)}` unique packet.\n"
        f"- Checked invariant feature class counts: `{feature_counts}`.\n"
        f"- Target replay: `q^5={product.numerator}/{product.denominator}`, eta residual `{eta_from_product(product, SUPPORT_SIZE)-STRICT_TARGET_ETA:.3e}`.\n"
        "- Honest read: every D12-invariant score is constant on the transitive anchored orbit, so unique source/orientation selection needs non-invariant strict-source data.\n"
        "- No false pass: no strict symmetry-breaking theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
