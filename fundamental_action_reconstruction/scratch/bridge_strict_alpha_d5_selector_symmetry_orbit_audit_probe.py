#!/usr/bin/env python3
"""Scratch probe: symmetry orbit audit for the d5 support/assignment selector stack.

The anchored-endpoint packet showed that an explicit source and orientation can
select one endpoint-block assignment after the distance-5 support and balanced
ledger are supplied.  This packet asks the remaining symmetry question: how much
selector data is that actually adding?

Result: the final anchored configurations form one free dihedral D_12 orbit of
size 24.  The unanchored distance-5 support orbit has 12 supports, and adding an
orientation/source endpoint doubles the ordered representatives to 24.  The
chosen value configuration has trivial D_12 stabilizer, so no internal symmetry
quotient selects a unique representative.  A source/orientation anchor is exactly
extra selector data, not a strict-core discharge of QW-2191.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ANCHORED = HERE / "bridge_strict_alpha_d5_anchored_endpoint_assignment_selector_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_d5_selector_symmetry_orbit_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_d5_selector_symmetry_orbit_audit_report.md"

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


def dihedral_action_on_support(support: tuple[int, ...], shift: int, reflect: bool) -> tuple[int, ...]:
    return tuple(sorted(((-node if reflect else node) + shift) % Z12_NODE_COUNT for node in support))


def all_anchored_config_rows() -> list[dict[str, Any]]:
    rows = []
    for source in range(Z12_NODE_COUNT):
        for orientation in (-1, 1):
            support = ordered_d5_support(source, orientation)
            rows.append(
                {
                    "source": source,
                    "orientation": orientation,
                    "ordered_support": list(support),
                    "canonical_support": list(sorted(support)),
                    "assignment": list(FORWARD_ASSIGNMENT),
                    "value_configuration": list(value_configuration(source, orientation)),
                }
            )
    return rows


def orbit_and_stabilizer() -> dict[str, Any]:
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
    support_base = tuple(sorted(ordered_d5_support(0, -1)))
    support_orbit = {
        dihedral_action_on_support(support_base, shift, reflect)
        for shift in range(Z12_NODE_COUNT)
        for reflect in (False, True)
    }
    return {
        "dihedral_group_order": 2 * Z12_NODE_COUNT,
        "base_configuration": list(base),
        "value_configuration_orbit_size": len(orbit),
        "value_configuration_stabilizer": stabilizer,
        "value_configuration_stabilizer_size": len(stabilizer),
        "support_orbit_size_without_orientation_assignment": len(support_orbit),
        "support_orbit": [list(support) for support in sorted(support_orbit)],
    }


def main() -> None:
    anchored_report = load_json(ANCHORED)
    product = Fraction(2 ** TARGET_BINARY_EXPONENT, DENOMINATOR ** SUPPORT_SIZE)
    rows = all_anchored_config_rows()
    orbit = orbit_and_stabilizer()
    unique_configs = {tuple(row["value_configuration"]) for row in rows}
    unique_supports = {tuple(row["canonical_support"]) for row in rows}

    report = {
        "status": "OPEN_STRICT_ALPHA_D5_SELECTOR_SYMMETRY_ORBIT_AUDIT_NO_STRICT_SELECTOR_DISCHARGE",
        "result_kind": "SCRATCH_STRICT_ALPHA_D5_SELECTOR_SYMMETRY_ORBIT_AUDIT_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "d5_anchored_endpoint_assignment_selector": str(ANCHORED.relative_to(ROOT)),
        },
        "previous_anchored_selector_replay": {
            "result_kind": anchored_report["result_kind"],
            "all_forward_assignments_unique": anchored_report["orbit_scan"]["all_forward_assignments_unique"],
            "forward_selected_assignment_set": anchored_report["orbit_scan"]["forward_selected_assignment_set"],
        },
        "target_identity_replay": {
            "q_power_product": f"{product.numerator}/{product.denominator}",
            "eta_residual_vs_9_5": eta_from_product(product, SUPPORT_SIZE) - STRICT_TARGET_ETA,
            "forward_assignment": list(FORWARD_ASSIGNMENT),
        },
        "anchored_configuration_enumeration": {
            "rows": rows,
            "row_count": len(rows),
            "unique_value_configuration_count": len(unique_configs),
            "unique_unoriented_support_count": len(unique_supports),
            "source_count": Z12_NODE_COUNT,
            "orientation_count_per_source": 2,
        },
        "dihedral_orbit_audit": orbit,
        "symmetry_verdict": {
            "anchored_rows_equal_full_dihedral_orbit": len(unique_configs) == orbit["value_configuration_orbit_size"] == orbit["dihedral_group_order"],
            "value_configuration_has_trivial_stabilizer": orbit["value_configuration_stabilizer_size"] == 1,
            "unoriented_support_orbit_has_12_members": orbit["support_orbit_size_without_orientation_assignment"] == Z12_NODE_COUNT,
            "interpretation": "The source/orientation anchor chooses one representative from a free D12 orbit; it is not obtained by quotienting an internal symmetry invariant.",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": True,
            "content": "The complete conditional d5 selector stack reduces the remaining ambiguity to a free D12 orbit of 24 anchored representatives.",
            "why_this_is_more_proof_like": "The probe explicitly enumerates all source/orientation value configurations and verifies the D12 orbit size and stabilizer.",
            "why_this_is_not_enough": "A free orbit has no canonical representative without extra source/orientation data, so QW-2191 is not discharged.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "The source/orientation anchor is confirmed as extra selector data, not a strict-core consequence of D12 symmetry.",
            "No theorem derives the endpoint-bias/source-moment action from strict nadsoliton geometry.",
            "No theorem derives branch count m=5, total exponent n=8, denominator 3, or binary-rescale quotient.",
            "No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged; this packet quantifies the remaining free orbit.",
            "No theorem derives eta=9/5 without the branch model, denominator/quotient convention, support premise, ledger selector, source/orientation premise, and assignment premise.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "A strict selector would need an internal source/orientation term that chooses one representative from the free D12 orbit; otherwise the anchor remains an explicit non-strict premise.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha d5 selector symmetry orbit audit probe\n\n"
        "Status: symmetry orbit audit for the conditional d5 selector stack; no strict selector discharge.\n\n"
        f"- Anchored rows: `{len(rows)}`; unique value configurations: `{len(unique_configs)}`; unoriented supports: `{len(unique_supports)}`.\n"
        f"- D12 orbit size of base value configuration: `{orbit['value_configuration_orbit_size']}` with stabilizer size `{orbit['value_configuration_stabilizer_size']}`.\n"
        f"- Support orbit size without orientation/assignment: `{orbit['support_orbit_size_without_orientation_assignment']}`.\n"
        f"- Target replay: `q^5={product.numerator}/{product.denominator}`, eta residual `{eta_from_product(product, SUPPORT_SIZE)-STRICT_TARGET_ETA:.3e}`.\n"
        "- Honest read: source/orientation picks a representative from a free D12 orbit; no symmetry-internal canonical representative is exported.\n"
        "- No false pass: no strict source/orientation theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
