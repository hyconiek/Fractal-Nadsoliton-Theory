#!/usr/bin/env python3
"""Scratch probe: strict-gate Dirichlet action audit for the alpha ledger selector.

The previous majorization packet reduced the selector obligation to a source for
a symmetric strictly convex / Schur-convex branch action.  This packet asks a
more physical and computational question: does the current strict gate kernel
itself provide a plausible positive Dirichlet quadratic action on the Z_12 ring?

Candidate action audited here, with strict-side QW-2049 parameters only:

    K_strict_gate(d) = cos(omega*d+phi)/(1+d^eta)
    omega=0.18575, phi=0.16250, beta=1, eta=9/5

and

    E(x)=sum_{i<j} K_strict_gate(dist_Z12(i,j)) * (x_i-x_j)^2.

What this proves and does not prove:

* The six cyclic-distance weights are positive for the current strict tuple, so
  the associated weighted graph Laplacian is positive semidefinite with one
  constant zero mode and positive nonconstant modes.
* A bounded exact placement/assignment scan for the three canonical eta=9/5
  ledgers shows the global minimum among those ledger candidates is balanced
  (2,2,2,1,1).
* However, the same scan also selects clustered contiguous supports, not the
  fifth-step resonance window.  Thus the strict-gate Dirichlet candidate can
  support a convex ledger selector only with an additional branch-support/
  phase-placement theorem.  It does not discharge QW-2191 or close ToE.
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
MAJORIZATION = HERE / "bridge_strict_alpha_symmetric_convex_selector_majorization_certificate_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_strict_gate_dirichlet_action_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_strict_gate_dirichlet_action_audit_report.md"

Z12_NODE_COUNT = 12
TARGET_BRANCH_COUNT = 5
TARGET_BINARY_EXPONENT = 8
DENOMINATOR = 3
STRICT_OMEGA = 0.18575
STRICT_PHI = 0.16250
STRICT_BETA = 1.0
STRICT_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
CANONICAL_LEDGERS = [
    (4, 1, 1, 1, 1),
    (3, 2, 1, 1, 1),
    (2, 2, 2, 1, 1),
]
FIFTH_WINDOW_SUPPORT = (0, 7, 2, 9, 4)
CONTIGUOUS_SUPPORT = (0, 1, 2, 3, 4)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def strict_gate_weight(distance: int) -> float:
    return math.cos(STRICT_OMEGA * distance + STRICT_PHI) / (1.0 + STRICT_BETA * (distance**STRICT_ETA))


def cyclic_distance(left: int, right: int) -> int:
    raw = abs(left - right) % Z12_NODE_COUNT
    return min(raw, Z12_NODE_COUNT - raw)



def is_cyclic_contiguous_support(support: list[int] | tuple[int, ...]) -> bool:
    nodes = sorted(support)
    if not nodes:
        return False
    node_set = set(nodes)
    size = len(nodes)
    for start in range(Z12_NODE_COUNT):
        if {(start + offset) % Z12_NODE_COUNT for offset in range(size)} == node_set:
            return True
    return False

def unique_permutations(values: tuple[int, ...]) -> list[tuple[int, ...]]:
    return sorted(set(itertools.permutations(values)))


def eta_from_product(product: Fraction, branch_count: int) -> float:
    correction = float(product) ** (1.0 / branch_count)
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def weight_rows() -> list[dict[str, Any]]:
    return [
        {
            "distance": distance,
            "weight": strict_gate_weight(distance),
            "positive": strict_gate_weight(distance) > 0.0,
        }
        for distance in range(1, Z12_NODE_COUNT // 2 + 1)
    ]


def laplacian_eigenvalue(mode: int) -> float:
    # Circulant weighted complete graph on Z12 using cyclic-distance weights.
    # For N even, distances 1..N/2-1 occur twice and distance N/2 occurs once.
    total = 0.0
    for distance in range(1, Z12_NODE_COUNT // 2):
        total += 2.0 * strict_gate_weight(distance) * (1.0 - math.cos(2.0 * math.pi * mode * distance / Z12_NODE_COUNT))
    half = Z12_NODE_COUNT // 2
    total += strict_gate_weight(half) * (1.0 - math.cos(math.pi * mode))
    return total


def laplacian_spectrum_rows() -> list[dict[str, Any]]:
    return [
        {"mode": mode, "eigenvalue": laplacian_eigenvalue(mode)}
        for mode in range(Z12_NODE_COUNT)
    ]


def dirichlet_energy(nodes: tuple[int, ...], weights_on_nodes: tuple[int, ...]) -> float:
    values = [0.0] * Z12_NODE_COUNT
    for node, weight in zip(nodes, weights_on_nodes):
        values[node] = float(weight)
    total = 0.0
    for left in range(Z12_NODE_COUNT):
        for right in range(left + 1, Z12_NODE_COUNT):
            total += strict_gate_weight(cyclic_distance(left, right)) * (values[left] - values[right]) ** 2
    return total


def support_energy_summary(support: tuple[int, ...], ledger: tuple[int, ...]) -> dict[str, Any]:
    rows = []
    for perm in unique_permutations(ledger):
        rows.append({"weights": list(perm), "energy": dirichlet_energy(support, perm)})
    best = min(rows, key=lambda row: row["energy"])
    worst = max(rows, key=lambda row: row["energy"])
    return {
        "support": list(support),
        "ledger": list(ledger),
        "permutation_count": len(rows),
        "min_energy": best["energy"],
        "max_energy": worst["energy"],
        "mean_energy": sum(row["energy"] for row in rows) / len(rows),
        "best_assignment": best,
        "worst_assignment": worst,
    }


def global_assignment_scan() -> dict[str, Any]:
    ledger_summaries = []
    global_rows = []
    total_rows = 0
    for ledger in CANONICAL_LEDGERS:
        best = None
        worst = None
        energies = []
        for support in itertools.combinations(range(Z12_NODE_COUNT), TARGET_BRANCH_COUNT):
            for perm in unique_permutations(ledger):
                energy = dirichlet_energy(support, perm)
                row = {"ledger": list(ledger), "support": list(support), "weights": list(perm), "energy": energy}
                total_rows += 1
                energies.append(energy)
                if best is None or energy < best["energy"]:
                    best = row
                if worst is None or energy > worst["energy"]:
                    worst = row
        ledger_summaries.append(
            {
                "ledger": list(ledger),
                "rows_scanned": len(energies),
                "min_energy": best["energy"],
                "max_energy": worst["energy"],
                "mean_energy": sum(energies) / len(energies),
                "best_row": best,
                "worst_row": worst,
            }
        )
        global_rows.append(best)
    global_best = min(global_rows, key=lambda row: row["energy"])
    return {
        "total_rows_scanned": total_rows,
        "ledger_summaries": ledger_summaries,
        "global_best_among_canonical_ledgers": global_best,
        "global_best_is_balanced_ledger": global_best["ledger"] == [2, 2, 2, 1, 1],
        "global_best_support_is_cyclic_contiguous": is_cyclic_contiguous_support(global_best["support"]),
        "global_best_support_is_contiguous_example": global_best["support"] == list(CONTIGUOUS_SUPPORT),
        "global_best_support_is_fifth_window": global_best["support"] == list(FIFTH_WINDOW_SUPPORT),
    }


def main() -> None:
    majorization_report = load_json(MAJORIZATION)
    weights = weight_rows()
    spectrum = laplacian_spectrum_rows()
    nonzero_modes = [row for row in spectrum if row["mode"] != 0]
    assignment_scan = global_assignment_scan()
    product = Fraction(2 ** TARGET_BINARY_EXPONENT, DENOMINATOR ** TARGET_BRANCH_COUNT)
    fifth_window_rows = [support_energy_summary(FIFTH_WINDOW_SUPPORT, ledger) for ledger in CANONICAL_LEDGERS]
    contiguous_rows = [support_energy_summary(CONTIGUOUS_SUPPORT, ledger) for ledger in CANONICAL_LEDGERS]

    report = {
        "status": "OPEN_STRICT_ALPHA_STRICT_GATE_DIRICHLET_ACTION_AUDIT_NO_SUPPORT_SELECTOR_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_STRICT_GATE_DIRICHLET_ACTION_AUDIT_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "symmetric_convex_selector_majorization_certificate": str(MAJORIZATION.relative_to(ROOT)),
        },
        "previous_majorization_replay": {
            "result_kind": majorization_report["result_kind"],
            "balanced_ledger": majorization_report["target_eta_9_5_certificate"]["balanced_ledger"],
            "q_power_product": majorization_report["target_eta_9_5_certificate"]["q_power_product"],
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
            "eta_residual_vs_9_5": eta_from_product(product, TARGET_BRANCH_COUNT) - STRICT_ETA,
            "canonical_ledgers": [list(ledger) for ledger in CANONICAL_LEDGERS],
        },
        "positive_weight_audit": {
            "distance_weight_rows": weights,
            "all_cyclic_distance_weights_positive": all(row["positive"] for row in weights),
        },
        "laplacian_psd_audit": {
            "spectrum_rows": spectrum,
            "zero_mode_eigenvalue": spectrum[0]["eigenvalue"],
            "minimum_nonconstant_eigenvalue": min(row["eigenvalue"] for row in nonzero_modes),
            "all_nonconstant_modes_positive": all(row["eigenvalue"] > 0.0 for row in nonzero_modes),
            "interpretation": "Positive strict-gate distance weights give a connected weighted Z12 Dirichlet quadratic form with only the constant zero mode.",
        },
        "canonical_ledger_assignment_scan": assignment_scan,
        "support_comparison": {
            "fifth_window_support": list(FIFTH_WINDOW_SUPPORT),
            "contiguous_support": list(CONTIGUOUS_SUPPORT),
            "fifth_window_rows": fifth_window_rows,
            "contiguous_rows": contiguous_rows,
            "verdict": "The strict-gate Dirichlet candidate supports the balanced ledger in the bounded scan but prefers clustered supports over the fifth-step resonance window, so a separate support/phase-placement theorem is still required.",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": bool(
                all(row["positive"] for row in weights)
                and all(row["eigenvalue"] > 0.0 for row in nonzero_modes)
                and assignment_scan["global_best_is_balanced_ledger"]
            ),
            "content": "The current strict-gate tuple can be used as a candidate positive Dirichlet quadratic source for a convex ledger selector, and the bounded scan selects (2,2,2,1,1) among the canonical eta=9/5 ledgers.",
            "why_this_is_more_proof_like": "It tests the actual strict-gate parameters, verifies positivity/PSD by the Z12 circulant spectrum, and exhaustively scans all 5-node supports and assignments for the three canonical ledgers.",
            "why_this_is_not_enough": "The scan also exposes a support-selection problem: the strict-gate Dirichlet energy prefers contiguous/clustered supports, not the fifth-step resonance support, and no theorem derives this action or support rule from strict nadsoliton geometry.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "The strict-gate Dirichlet quadratic action is audited as a candidate construction, not derived as a strict nadsoliton action theorem.",
            "The bounded assignment scan selects the balanced ledger among canonical ledgers but does not derive branch count m=5, total exponent n=8, denominator 3, or binary-rescale quotient.",
            "The support scan prefers clustered supports over the fifth-step resonance window; no support/phase-placement theorem is exported.",
            "No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged.",
            "No theorem derives eta=9/5 without the branch model, denominator/quotient convention, Dirichlet action premise, and support/phase-placement premise.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "If pursuing a physical selector, prove or reject a strict nadsoliton action/source whose reduced branch-sector Dirichlet form and support rule are both exported; otherwise keep this as a candidate action audit with an explicit support-selection obstruction.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha strict-gate Dirichlet action audit probe\n\n"
        "Status: strict-gate Dirichlet candidate action audit; no support selector theorem.\n\n"
        f"- Strict tuple: `omega={STRICT_OMEGA}`, `phi={STRICT_PHI}`, `beta={STRICT_BETA}`, `eta={STRICT_ETA}`.\n"
        f"- Positivity: all six cyclic-distance weights are positive: `{all(row['positive'] for row in weights)}`; nonconstant Laplacian modes are positive: `{all(row['eigenvalue'] > 0.0 for row in nonzero_modes)}`.\n"
        f"- Assignment scan: `{assignment_scan['total_rows_scanned']}` support/assignment rows; global best ledger is `{assignment_scan['global_best_among_canonical_ledgers']['ledger']}` with support `{assignment_scan['global_best_among_canonical_ledgers']['support']}`.\n"
        f"- Target replay: `q^5={product.numerator}/{product.denominator}`, eta residual `{eta_from_product(product, TARGET_BRANCH_COUNT)-STRICT_ETA:.3e}`.\n"
        "- Honest read: strict-gate Dirichlet energy supports the balanced ledger as a candidate, but its support preference is clustered rather than fifth-step resonant.\n"
        "- No false pass: no derived strict action theorem, no support/phase-placement theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
