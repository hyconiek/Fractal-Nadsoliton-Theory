#!/usr/bin/env python3
"""Scratch probe: strict-gate Dirichlet support obstruction certificate.

The previous strict-gate Dirichlet audit found a useful but dangerous fact: the
current strict gate tuple can support a convex balanced-ledger selector, but the
same candidate action prefers clustered supports rather than the fifth-step
resonance window.  This packet turns that warning into a more explicit finite
certificate.

We keep the candidate strict-side action

    E(x)=sum_{i<j} K_strict_gate(dist_Z12(i,j)) * (x_i-x_j)^2

with QW-2049 strict tuple omega=0.18575, phi=0.16250, beta=1, eta=9/5.  Then we
separate two support questions:

1. support-only indicator energy, where x_i=1 on a 5-node support and 0 outside;
2. balanced-ledger support energy, where x carries the ledger (2,2,2,1,1) and
   we minimize over assignments on each support.

Both finite scans give the same obstruction pattern: cyclic-contiguous 5-blocks
are global minima, while the fifth-step resonance support is a global maximum
orbit for the balanced-ledger support scan and the support-only scan.  Therefore
this Dirichlet candidate cannot by itself justify the phase-resonant support; it
requires an extra support/phase-placement theorem or a different action term.
"""
from __future__ import annotations

import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DIRICHLET_AUDIT = HERE / "bridge_strict_alpha_strict_gate_dirichlet_action_audit_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_dirichlet_support_obstruction_certificate_report.json"
OUT_MD = HERE / "bridge_strict_alpha_dirichlet_support_obstruction_certificate_report.md"

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
STRICT_OMEGA = 0.18575
STRICT_PHI = 0.16250
STRICT_BETA = 1.0
STRICT_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
DENOMINATOR = 3
TARGET_BINARY_EXPONENT = 8
BALANCED_LEDGER = (2, 2, 2, 1, 1)
FIFTH_STEP_SUPPORT_ORDERED = (0, 7, 2, 9, 4)
FIFTH_STEP_SUPPORT = tuple(sorted(FIFTH_STEP_SUPPORT_ORDERED))
CONTIGUOUS_EXAMPLE = (0, 1, 2, 3, 4)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def strict_gate_weight(distance: int) -> float:
    return math.cos(STRICT_OMEGA * distance + STRICT_PHI) / (1.0 + STRICT_BETA * (distance**STRICT_ETA))


def cyclic_distance(left: int, right: int) -> int:
    raw = abs(left - right) % Z12_NODE_COUNT
    return min(raw, Z12_NODE_COUNT - raw)


def eta_from_product(product: Fraction, branch_count: int) -> float:
    correction = float(product) ** (1.0 / branch_count)
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def unique_permutations(values: tuple[int, ...]) -> list[tuple[int, ...]]:
    return sorted(set(itertools.permutations(values)))


def support_gap_signature(support: tuple[int, ...]) -> tuple[int, ...]:
    nodes = sorted(support)
    gaps = []
    for index, node in enumerate(nodes):
        nxt = nodes[(index + 1) % len(nodes)]
        gaps.append((nxt - node) % Z12_NODE_COUNT)
    return tuple(sorted(gaps))


def is_cyclic_contiguous_support(support: tuple[int, ...]) -> bool:
    node_set = set(support)
    return any({(start + offset) % Z12_NODE_COUNT for offset in range(len(support))} == node_set for start in range(Z12_NODE_COUNT))


def values_from_support_assignment(support: tuple[int, ...], assignment: tuple[int, ...]) -> list[float]:
    values = [0.0] * Z12_NODE_COUNT
    for node, value in zip(support, assignment):
        values[node] = float(value)
    return values


def dirichlet_energy_values(values: list[float]) -> float:
    return sum(
        strict_gate_weight(cyclic_distance(left, right)) * (values[left] - values[right]) ** 2
        for left in range(Z12_NODE_COUNT)
        for right in range(left + 1, Z12_NODE_COUNT)
    )


def indicator_energy(support: tuple[int, ...]) -> float:
    values = [1.0 if node in set(support) else 0.0 for node in range(Z12_NODE_COUNT)]
    return dirichlet_energy_values(values)


def balanced_support_energy_summary(support: tuple[int, ...]) -> dict[str, Any]:
    rows = []
    for assignment in unique_permutations(BALANCED_LEDGER):
        values = values_from_support_assignment(support, assignment)
        rows.append({"assignment": list(assignment), "energy": dirichlet_energy_values(values)})
    best = min(rows, key=lambda row: row["energy"])
    worst = max(rows, key=lambda row: row["energy"])
    return {
        "support": list(support),
        "gap_signature": list(support_gap_signature(support)),
        "cyclic_contiguous": is_cyclic_contiguous_support(support),
        "permutation_count": len(rows),
        "min_energy": best["energy"],
        "max_energy": worst["energy"],
        "best_assignment": best,
        "worst_assignment": worst,
    }


def extremal_support_scan() -> dict[str, Any]:
    supports = list(itertools.combinations(range(Z12_NODE_COUNT), SUPPORT_SIZE))
    indicator_rows = []
    balanced_rows = []
    for support in supports:
        indicator_rows.append(
            {
                "support": list(support),
                "gap_signature": list(support_gap_signature(support)),
                "cyclic_contiguous": is_cyclic_contiguous_support(support),
                "energy": indicator_energy(support),
            }
        )
        balanced_rows.append(balanced_support_energy_summary(support))

    indicator_min = min(row["energy"] for row in indicator_rows)
    indicator_max = max(row["energy"] for row in indicator_rows)
    balanced_min = min(row["min_energy"] for row in balanced_rows)
    balanced_max = max(row["min_energy"] for row in balanced_rows)
    return {
        "support_count": len(supports),
        "balanced_assignment_rows_scanned": len(supports) * len(unique_permutations(BALANCED_LEDGER)),
        "indicator_energy": {
            "min_energy": indicator_min,
            "max_energy": indicator_max,
            "minimizer_count": sum(abs(row["energy"] - indicator_min) < 1e-12 for row in indicator_rows),
            "maximizer_count": sum(abs(row["energy"] - indicator_max) < 1e-12 for row in indicator_rows),
            "minimizers": [row for row in indicator_rows if abs(row["energy"] - indicator_min) < 1e-12],
            "maximizers": [row for row in indicator_rows if abs(row["energy"] - indicator_max) < 1e-12],
            "fifth_step_row": next(row for row in indicator_rows if tuple(row["support"]) == FIFTH_STEP_SUPPORT),
            "contiguous_example_row": next(row for row in indicator_rows if tuple(row["support"]) == CONTIGUOUS_EXAMPLE),
        },
        "balanced_ledger_min_assignment_energy": {
            "ledger": list(BALANCED_LEDGER),
            "min_energy": balanced_min,
            "max_of_support_minima": balanced_max,
            "minimizer_count": sum(abs(row["min_energy"] - balanced_min) < 1e-12 for row in balanced_rows),
            "maximizer_count": sum(abs(row["min_energy"] - balanced_max) < 1e-12 for row in balanced_rows),
            "minimizers": [row for row in balanced_rows if abs(row["min_energy"] - balanced_min) < 1e-12],
            "maximizers": [row for row in balanced_rows if abs(row["min_energy"] - balanced_max) < 1e-12],
            "fifth_step_row": next(row for row in balanced_rows if tuple(row["support"]) == FIFTH_STEP_SUPPORT),
            "contiguous_example_row": next(row for row in balanced_rows if tuple(row["support"]) == CONTIGUOUS_EXAMPLE),
        },
    }


def main() -> None:
    dirichlet_report = load_json(DIRICHLET_AUDIT)
    scan = extremal_support_scan()
    product = Fraction(2 ** TARGET_BINARY_EXPONENT, DENOMINATOR ** SUPPORT_SIZE)
    indicator = scan["indicator_energy"]
    balanced = scan["balanced_ledger_min_assignment_energy"]

    report = {
        "status": "OPEN_STRICT_ALPHA_DIRICHLET_SUPPORT_OBSTRUCTION_CERTIFICATE_NO_SUPPORT_SELECTOR_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_DIRICHLET_SUPPORT_OBSTRUCTION_CERTIFICATE_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "strict_gate_dirichlet_action_audit": str(DIRICHLET_AUDIT.relative_to(ROOT)),
        },
        "previous_dirichlet_audit_replay": {
            "result_kind": dirichlet_report["result_kind"],
            "global_best_ledger": dirichlet_report["canonical_ledger_assignment_scan"]["global_best_among_canonical_ledgers"]["ledger"],
            "global_best_support": dirichlet_report["canonical_ledger_assignment_scan"]["global_best_among_canonical_ledgers"]["support"],
            "global_best_support_is_cyclic_contiguous": dirichlet_report["canonical_ledger_assignment_scan"]["global_best_support_is_cyclic_contiguous"],
        },
        "strict_gate_parameters": {
            "kernel": "K_strict_gate(d)=cos(omega*d+phi)/(1+beta*d^eta)",
            "omega": STRICT_OMEGA,
            "phi": STRICT_PHI,
            "beta": STRICT_BETA,
            "eta": STRICT_ETA,
            "note": "strict-side working tuple only; no K_legacy_ont identification is used",
        },
        "target_identity_replay": {
            "q_power_product": f"{product.numerator}/{product.denominator}",
            "eta_residual_vs_9_5": eta_from_product(product, SUPPORT_SIZE) - STRICT_ETA,
            "balanced_ledger": list(BALANCED_LEDGER),
            "fifth_step_support_ordered": list(FIFTH_STEP_SUPPORT_ORDERED),
            "fifth_step_support_canonical": list(FIFTH_STEP_SUPPORT),
        },
        "support_obstruction_scan": scan,
        "obstruction_summary": {
            "indicator_contiguous_minimum": indicator["contiguous_example_row"]["energy"] == indicator["min_energy"],
            "indicator_fifth_step_is_maximum": indicator["fifth_step_row"]["energy"] == indicator["max_energy"],
            "balanced_contiguous_minimum": balanced["contiguous_example_row"]["min_energy"] == balanced["min_energy"],
            "balanced_fifth_step_is_maximum_of_support_minima": balanced["fifth_step_row"]["min_energy"] == balanced["max_of_support_minima"],
            "energy_gap_indicator_fifth_minus_contiguous": indicator["fifth_step_row"]["energy"] - indicator["contiguous_example_row"]["energy"],
            "energy_gap_balanced_fifth_minus_contiguous": balanced["fifth_step_row"]["min_energy"] - balanced["contiguous_example_row"]["min_energy"],
            "verdict": "The strict-gate Dirichlet candidate has a structural support-pressure obstruction: it clusters support instead of selecting the fifth-step resonance support.",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": True,
            "content": "The previous Dirichlet action candidate is sharpened into a negative support certificate: it can support the balanced ledger only after a separate support/phase-placement rule is supplied.",
            "why_this_is_more_proof_like": "The probe separates support-only and balanced-ledger scans, checks all 792 five-node supports, and shows the fifth-step support is extremal in the wrong direction for this action.",
            "why_this_is_not_enough": "This is not a derivation of the desired resonance support; it is evidence that the current Dirichlet candidate alone cannot be the missing strict-core support selector.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "The strict-gate Dirichlet quadratic action remains a candidate construction, not a derived strict nadsoliton action theorem.",
            "This packet proves a support obstruction for that candidate; it does not derive a replacement support/phase-placement theorem.",
            "No theorem derives branch count m=5, total exponent n=8, denominator 3, or binary-rescale quotient.",
            "No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged.",
            "No theorem derives eta=9/5 without the branch model, denominator/quotient convention, convex ledger selector, and support/phase-placement premise.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Either add a genuinely new strict support/phase-placement source that reverses the clustering pressure, or reject the plain strict-gate Dirichlet action as the standalone selector source.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Dirichlet support obstruction certificate probe\n\n"
        "Status: support obstruction certificate for the strict-gate Dirichlet candidate; no support selector theorem.\n\n"
        f"- Supports scanned: `{scan['support_count']}`; balanced assignment rows scanned: `{scan['balanced_assignment_rows_scanned']}`.\n"
        f"- Indicator support scan: contiguous example energy `{indicator['contiguous_example_row']['energy']:.12f}`, fifth-step energy `{indicator['fifth_step_row']['energy']:.12f}`.\n"
        f"- Balanced ledger scan: contiguous example min energy `{balanced['contiguous_example_row']['min_energy']:.12f}`, fifth-step min energy `{balanced['fifth_step_row']['min_energy']:.12f}`.\n"
        f"- Target replay: `q^5={product.numerator}/{product.denominator}`, eta residual `{eta_from_product(product, SUPPORT_SIZE)-STRICT_ETA:.3e}`.\n"
        "- Honest read: plain strict-gate Dirichlet energy clusters support; fifth-step resonance support needs an additional source or a different action.\n"
        "- No false pass: no support/phase-placement theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
