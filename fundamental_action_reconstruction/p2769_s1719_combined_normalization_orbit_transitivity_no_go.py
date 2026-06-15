#!/usr/bin/env python3
"""P2769/S1719: combined-normalization orbit-transitivity no-go.

P2768 ruled out monomial ratio invariants of the P1562 coefficients under the
combined open length/field/curvature normalization action.  This follow-up
attacks the next honest escape route: a genuinely non-monomial invariant.  The
same exponent action is a log-linear action of three positive gauge parameters
on the three positive coefficient coordinates.  Its coefficient/gauge matrix is
invertible, so the action is transitive on the positive coefficient octant: any
positive target triple is reached from the P1562 triple by a unique positive
gauge.  Therefore every invariant function on this open positive orbit is
constant.  This blocks non-monomial invariant rescue, but still does not export
a canonical normalization source or physical Lagrangian coupling theorem.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
P1562 = GEN / "p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json"
P2762 = GEN / "p2762_s1712_reference_cell_action_density_normalization_obstruction.json"
P2764 = GEN / "p2764_s1714_field_curvature_normalization_compatibility_obstruction.json"
P2768 = GEN / "p2768_s1718_combined_normalization_monomial_invariant_no_go.json"
OUT = GEN / "p2769_s1719_combined_normalization_orbit_transitivity_no_go.json"
MD = GEN / "p2769_s1719_combined_normalization_orbit_transitivity_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

COEFFICIENT_ORDER = ["lambda_sm_eff", "kappa_gr_eff", "epsilon_mix_eff"]
GAUGE_ORDER = ["ell_reference_length", "Z_phi_scalar_field", "Z_R_curvature"]
# Rows are coefficients, columns are gauges; log(c') = log(c) + ACTION_MATRIX @ log(g).
ACTION_MATRIX = [
    [0.0, -4.0, 0.0],
    [-2.0, 0.0, -1.0],
    [-1.0, -2.0, -1.0],
]

CONTENT_PATTERNS = {
    "p2768_next_step": r"P2768|non-monomial invariant|combined open normalization action|monomial-invariant",
    "p2762_length_open": r"P2762|reference-cell/action-density|positive reference length|canonical_reference_cell_theorem_exported",
    "p2764_field_curvature_open": r"P2764|field/curvature normalization|Z_phi|Z_R|field_curvature_normalization_theorem_exported",
    "invariant_boundary": r"invariant|normalization invariant|physical-coupling provenance|canonical normalization theorem",
    "closure_boundary": r"role-bearing L_total|ToE closure|bridge closure|role transfer|selector closure",
}

NEGATIVE_EXPORT_FLAGS = [
    "nonconstant_invariant_rescue_exported",
    "canonical_normalization_theorem_exported",
    "physical_coupling_provenance_theorem_exported",
    "role_bearing_ltotal_promoted",
    "selector_closure_exported",
    "bridge_closure_exported",
    "role_transfer_started",
    "toe_closure_exported",
]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def run_rg(pattern: str) -> list[str]:
    cmd = ["rg", "-n", "--glob", "!generated/*.json", pattern, "AGENTS.md", "fundamental_action_reconstruction"]
    proc = subprocess.run(cmd, cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if proc.returncode not in (0, 1):
        raise RuntimeError(f"rg failed for {pattern!r}: {proc.stderr}")
    return [line for line in proc.stdout.splitlines() if line.strip()]


def evidence_scan() -> dict[str, Any]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        hits = run_rg(pattern)
        rows.append({"lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return {"row_count": len(rows), "rows": rows, "hit_counts": {r["lane"]: r["hit_count"] for r in rows}, "all_patterns_have_hits": all(r["hit_count"] > 0 for r in rows)}


def det3(matrix: list[list[float]]) -> float:
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def replace_column(matrix: list[list[float]], col: int, values: list[float]) -> list[list[float]]:
    return [[values[r] if c == col else matrix[r][c] for c in range(3)] for r in range(3)]


def solve3(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    d = det3(matrix)
    if abs(d) < 1e-12:
        raise ValueError("singular 3x3 matrix")
    return [det3(replace_column(matrix, col, rhs)) / d for col in range(3)]


def transform(base: dict[str, float], gauges: dict[str, float]) -> dict[str, float]:
    ell = gauges["ell_reference_length"]
    z_phi = gauges["Z_phi_scalar_field"]
    z_r = gauges["Z_R_curvature"]
    return {
        "lambda_sm_eff": base["lambda_sm_eff"] * z_phi ** -4,
        "kappa_gr_eff": base["kappa_gr_eff"] * ell ** -2 * z_r ** -1,
        "epsilon_mix_eff": base["epsilon_mix_eff"] * ell ** -1 * z_phi ** -2 * z_r ** -1,
    }


def orbit_transitivity_witness(p1562: dict[str, Any]) -> dict[str, Any]:
    base = {key: float(p1562["derived_lagrangian_coefficients"][key]) for key in COEFFICIENT_ORDER}
    determinant = det3(ACTION_MATRIX)
    target_multipliers = [
        {"lambda_sm_eff": 2.0, "kappa_gr_eff": 0.5, "epsilon_mix_eff": 3.0},
        {"lambda_sm_eff": 0.25, "kappa_gr_eff": 4.0, "epsilon_mix_eff": 0.75},
        {"lambda_sm_eff": 1.5, "kappa_gr_eff": 1.25, "epsilon_mix_eff": 0.4},
        {"lambda_sm_eff": 8.0, "kappa_gr_eff": 0.125, "epsilon_mix_eff": 2.0},
    ]
    rows = []
    max_relative_error = 0.0
    for multipliers in target_multipliers:
        target = {key: base[key] * multipliers[key] for key in COEFFICIENT_ORDER}
        rhs = [math.log(target[key] / base[key]) for key in COEFFICIENT_ORDER]
        log_gauge = solve3(ACTION_MATRIX, rhs)
        gauges = {key: math.exp(value) for key, value in zip(GAUGE_ORDER, log_gauge)}
        reconstructed = transform(base, gauges)
        rel_errors = {key: abs(reconstructed[key] - target[key]) / target[key] for key in COEFFICIENT_ORDER}
        max_relative_error = max(max_relative_error, max(rel_errors.values()))
        rows.append({
            "target_multipliers": multipliers,
            "target_coefficients": target,
            "solved_log_gauge": {key: value for key, value in zip(GAUGE_ORDER, log_gauge)},
            "solved_positive_gauge": gauges,
            "reconstructed_coefficients": reconstructed,
            "relative_errors": rel_errors,
            "max_relative_error": max(rel_errors.values()),
        })
    return {
        "coefficient_order": COEFFICIENT_ORDER,
        "gauge_order": GAUGE_ORDER,
        "action_matrix_coefficients_by_gauges": ACTION_MATRIX,
        "determinant": determinant,
        "full_rank_over_R": abs(determinant) > 1e-12,
        "base_coefficients": base,
        "target_rows": rows,
        "target_row_count": len(rows),
        "all_sampled_targets_reached": max_relative_error < 1e-10,
        "max_relative_error": max_relative_error,
        "theorem_statement": "Because the log-linear coefficient/gauge action matrix is invertible, the positive gauge group acts transitively on the positive coefficient octant; every invariant function on this orbit is constant.",
    }


def acceptance_matrix(scan: dict[str, Any], witness: dict[str, Any], p2762: dict[str, Any], p2764: dict[str, Any], p2768: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_evidence_present": scan["all_patterns_have_hits"],
        "p2762_reference_cell_still_open": p2762.get("acceptance_matrix", {}).get("accepted_as_reference_cell_action_density_theorem") is False,
        "p2764_field_curvature_still_open": p2764.get("acceptance_matrix", {}).get("accepted_as_field_curvature_normalization_theorem") is False,
        "p2768_monomial_rescue_blocked": p2768.get("acceptance_matrix", {}).get("accepted_as_monomial_ratio_rescue") is False,
        "action_matrix_full_rank": witness["full_rank_over_R"],
        "sampled_positive_targets_reached": witness["all_sampled_targets_reached"],
        "nonconstant_invariant_rescue_exported": False,
        "canonical_normalization_theorem_exported": False,
        "physical_coupling_provenance_theorem_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_nonconstant_invariant_rescue": False,
        "accepted_as_canonical_normalization_theorem": False,
        "accepted_as_physical_coupling_provenance_theorem": False,
        "accepted_as_ltotal_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "The combined positive normalization action is transitive on the positive coefficient octant, so every invariant function of the three P1562 coefficients alone is constant on the open orbit.  This blocks non-monomial invariant rescue but does not choose a canonical gauge or export physical-coupling provenance.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    witness = payload["combined_normalization_orbit_transitivity_witness"]
    lines = [
        "# P2769/S1719 combined-normalization orbit-transitivity no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Log-linear action",
        f"- coefficient_order={witness['coefficient_order']}",
        f"- gauge_order={witness['gauge_order']}",
        f"- action_matrix={witness['action_matrix_coefficients_by_gauges']}",
        f"- determinant={witness['determinant']}",
        f"- full_rank_over_R={witness['full_rank_over_R']}",
        "",
        "## Target reachability witness",
        f"- target_row_count={witness['target_row_count']}",
        f"- all_sampled_targets_reached={witness['all_sampled_targets_reached']}",
        f"- max_relative_error={witness['max_relative_error']}",
        "",
        "## Decision",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p1562 = read_json(P1562)
    p2762 = read_json(P2762)
    p2764 = read_json(P2764)
    p2768 = read_json(P2768)
    scan = evidence_scan()
    witness = orbit_transitivity_witness(p1562)
    acceptance = acceptance_matrix(scan, witness, p2762, p2764, p2768)
    payload = {
        "status": "P2769_COMBINED_NORMALIZATION_ORBIT_TRANSITIVITY_NO_GO_NO_CLOSURE",
        "input_hashes": {"P1562": sha(P1562), "P2762": sha(P2762), "P2764": sha(P2764), "P2768": sha(P2768)},
        "input_statuses": {"P1562": p1562.get("status"), "P2762": p2762.get("status"), "P2764": p2764.get("status"), "P2768": p2768.get("status")},
        "audited_object": "non-monomial invariant rescue under the combined positive normalization action",
        "content_evidence_scan": scan,
        "combined_normalization_orbit_transitivity_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not try to rescue physical-coupling provenance with any invariant function of only lambda_sm_eff, kappa_gr_eff, and epsilon_mix_eff under the still-open positive normalization gauges: P2769 shows that all such invariants are constant on the open positive orbit.  The next honest move must supply an external canonical normalization source/theorem, e.g. a selector-free reference-cell/action-density source or field/curvature normalization law, and then rerun bounded acceptance; otherwise preserve the P2697-P2769 no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2769/S1719 combined-normalization orbit-transitivity no-go", "## P2769/S1719 combined-normalization orbit-transitivity no-go\n\n`P2769/S1719` attacks the non-monomial invariant escape route left after P2768.  The log-linear action of positive reference length, scalar-field normalization, and curvature normalization on the three P1562 coefficients has an invertible coefficient/gauge matrix, so it is transitive on the positive coefficient octant.  Sampled positive target triples are reached by unique positive gauges with negligible reconstruction error.  Therefore every invariant function of `lambda_sm_eff`, `kappa_gr_eff`, and `epsilon_mix_eff` alone is constant on the open orbit.  This blocks non-monomial invariant rescue, but it is not a canonical normalization theorem and exports no physical-coupling provenance, role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2769/S1719 orbit-transitivity Ltotal guard", "## P2769/S1719 orbit-transitivity Ltotal guard\n\n`P2769/S1719` adds no variational source term.  It shows the open positive normalization gauges act transitively on the positive triple of P1562 coefficients, so invariant-taking cannot promote `lambda_sm_eff`, `kappa_gr_eff`, or `epsilon_mix_eff` to physical Lagrangian couplings or role-bearing `L_total`.\n")
    append_once(AGENTS, "Current combined-normalization orbit-transitivity no-go guardrail (P2769/S1719, 2026-06-15)", "## Current combined-normalization orbit-transitivity no-go guardrail (P2769/S1719, 2026-06-15)\n\n- P2769 attacks the non-monomial invariant escape route after P2768: the positive length/field/curvature normalization action on the three P1562 coefficients is log-linearly invertible and transitive on the positive coefficient octant.\n- Every invariant function of `lambda_sm_eff`, `kappa_gr_eff`, and `epsilon_mix_eff` alone is therefore constant on the open orbit; invariant-taking cannot rescue physical-coupling provenance while P2762/P2764 normalizations remain open.\n- Do not promote this no-go to canonical normalization, physical-coupling provenance, role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.  A next admissible move must supply an external canonical normalization source/theorem and then rerun bounded acceptance; otherwise preserve the P2697-P2769 no-closure certificate.\n")
    return payload


if __name__ == "__main__":
    main()
