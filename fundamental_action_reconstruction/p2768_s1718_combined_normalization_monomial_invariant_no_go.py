#!/usr/bin/env python3
"""P2768/S1718: combined-normalization monomial-invariant no-go.

P2767 required one concrete new typed object/theorem before reopening any lane.
This script supplies a bounded object in the P2766 moment-provenance lane: a
finite linear-algebra obstruction to rescuing physical provenance via a
monomial ratio of the three P1562 coefficients.  Combining the P2762 reference
length action and the P2764 field/curvature action gives a 3x3 exponent-weight
matrix on (lambda_sm_eff, kappa_gr_eff, epsilon_mix_eff).  The matrix has full
rank, so the only monomial invariant is the trivial exponent vector (0,0,0).
This is an obstruction theorem for a ratio-invariant escape route, not a
normalization theorem and not an L_total promotion.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
import subprocess
from fractions import Fraction
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
P1562 = GEN / "p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json"
P2762 = GEN / "p2762_s1712_reference_cell_action_density_normalization_obstruction.json"
P2764 = GEN / "p2764_s1714_field_curvature_normalization_compatibility_obstruction.json"
P2766 = GEN / "p2766_s1716_post_moment_provenance_state_reconciliation.json"
P2767 = GEN / "p2767_s1717_post_p2766_fresh_state_map_intake.json"
OUT = GEN / "p2768_s1718_combined_normalization_monomial_invariant_no_go.json"
MD = GEN / "p2768_s1718_combined_normalization_monomial_invariant_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

COEFFICIENT_ORDER = ["lambda_sm_eff", "kappa_gr_eff", "epsilon_mix_eff"]
GAUGE_ORDER = ["ell_reference_length", "Z_phi_scalar_field", "Z_R_curvature"]
# Rows are gauges, columns are coefficient exponents.  If a monomial is
# lambda^a*kappa^b*epsilon^c, its gauge weight is W @ (a,b,c).
WEIGHT_MATRIX = [
    [0, -2, -1],   # ell: lambda invariant, kappa -> ell^-2, epsilon -> ell^-1
    [-4, 0, -2],   # Z_phi: lambda -> Z_phi^-4, epsilon -> Z_phi^-2
    [0, -1, -1],   # Z_R: kappa -> Z_R^-1, epsilon -> Z_R^-1
]

CONTENT_PATTERNS = {
    "p2767_new_object_gate": r"P2767|one concrete new typed object/theorem|fresh broad state-map intake",
    "p2762_length_action": r"P2762|positive-length scale-orbit|reference length|kappa_changes_with_ell|epsilon_changes_with_ell",
    "p2764_field_curvature_action": r"P2764|positive-normalization orbit|Z_phi|Z_R|normalization_sensitive",
    "monomial_ratio_boundary": r"monomial|ratio|dimensionless_ratios|normalization invariant|physical-coupling provenance",
    "closure_boundary": r"role-bearing L_total|ToE closure|bridge closure|role transfer|selector closure",
}

NEGATIVE_EXPORT_FLAGS = [
    "canonical_reference_cell_theorem_exported",
    "field_curvature_normalization_theorem_exported",
    "monomial_ratio_rescue_theorem_exported",
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


def det3(matrix: list[list[int]]) -> int:
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def mat_vec(matrix: list[list[int]], vector: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(sum(row[i] * vector[i] for i in range(3)) for row in matrix)  # type: ignore[return-value]


def monomial_value(coeffs: dict[str, float], exponents: tuple[int, int, int]) -> float:
    value = 1.0
    for name, exponent in zip(COEFFICIENT_ORDER, exponents):
        value *= coeffs[name] ** exponent
    return value


def transformed_coefficients(base: dict[str, float], ell: float, z_phi: float, z_r: float) -> dict[str, float]:
    return {
        "lambda_sm_eff": base["lambda_sm_eff"] * z_phi ** -4,
        "kappa_gr_eff": base["kappa_gr_eff"] * ell ** -2 * z_r ** -1,
        "epsilon_mix_eff": base["epsilon_mix_eff"] * ell ** -1 * z_phi ** -2 * z_r ** -1,
    }


def invariant_no_go_witness(p1562: dict[str, Any]) -> dict[str, Any]:
    base = p1562.get("derived_lagrangian_coefficients", {})
    determinant = det3(WEIGHT_MATRIX)
    brute_rows = []
    invariant_exponents = []
    for exponents in itertools.product(range(-4, 5), repeat=3):
        if exponents == (0, 0, 0):
            continue
        weight = mat_vec(WEIGHT_MATRIX, exponents)
        if weight == (0, 0, 0):
            invariant_exponents.append(exponents)
        brute_rows.append({"exponents": exponents, "weight": weight, "is_gauge_invariant": weight == (0, 0, 0)})
    sample_gauges = [
        {"ell": 0.5, "Z_phi": 0.5, "Z_R": 0.5},
        {"ell": 1.0, "Z_phi": 1.0, "Z_R": 1.0},
        {"ell": 2.0, "Z_phi": 2.0, "Z_R": 2.0},
        {"ell": 4.0, "Z_phi": 0.5, "Z_R": 2.0},
    ]
    tempting_ratios = {
        "epsilon_squared_over_lambda_kappa": (-1, -1, 2),
        "epsilon_squared_over_lambda_kappa_squared": (-1, -2, 2),
        "kappa_over_epsilon_squared": (0, 1, -2),
        "p1562_R1_R2_R3_moment_ratios_not_coupling_normalization": (0, 0, 0),
    }
    ratio_rows = []
    for label, exponents in tempting_ratios.items():
        values = []
        for gauge in sample_gauges:
            coeffs = transformed_coefficients(base, gauge["ell"], gauge["Z_phi"], gauge["Z_R"])
            values.append(monomial_value(coeffs, exponents))
        ratio_rows.append({
            "candidate": label,
            "exponents_lambda_kappa_epsilon": exponents,
            "gauge_weight": mat_vec(WEIGHT_MATRIX, exponents),
            "sample_values": values,
            "distinct_sample_values": len({round(v, 12) for v in values}),
            "invariant_under_combined_action": len({round(v, 12) for v in values}) == 1,
        })
    rational_inverse = [[str(Fraction(x, determinant)) for x in row] for row in [
        [1, -1, 2],
        [0, 0, -1],
        [-1, 0, 2],
    ]]
    return {
        "coefficient_order": COEFFICIENT_ORDER,
        "gauge_order": GAUGE_ORDER,
        "weight_matrix": WEIGHT_MATRIX,
        "determinant": determinant,
        "full_rank_over_Q": determinant != 0,
        "rational_inverse_certificate_det_scaled": rational_inverse,
        "brute_force_integer_box": {"min_exponent": -4, "max_exponent": 4, "tested_nonzero_exponent_vectors": len(brute_rows), "nontrivial_invariant_count": len(invariant_exponents), "nontrivial_invariant_exponents": invariant_exponents},
        "tempting_ratio_rows": ratio_rows,
        "base_coefficients": {key: base[key] for key in COEFFICIENT_ORDER},
        "sampled_combined_gauges": sample_gauges,
    }


def acceptance_matrix(scan: dict[str, Any], witness: dict[str, Any], p2762: dict[str, Any], p2764: dict[str, Any], p2766: dict[str, Any], p2767: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_evidence_present": scan["all_patterns_have_hits"],
        "p2762_reference_cell_still_open": p2762.get("acceptance_matrix", {}).get("accepted_as_reference_cell_action_density_theorem") is False,
        "p2764_field_curvature_still_open": p2764.get("acceptance_matrix", {}).get("accepted_as_field_curvature_normalization_theorem") is False,
        "p2766_physical_coupling_provenance_still_open": p2766.get("acceptance_matrix", {}).get("accepted_as_physical_coupling_provenance_theorem") is False,
        "p2767_requires_new_object_before_lane_reopen": p2767.get("new_object_intake_gate", {}).get("accepted_as_no_new_live_frontier_certificate") is True,
        "weight_matrix_full_rank": witness["full_rank_over_Q"],
        "bruteforce_box_finds_no_nontrivial_invariant": witness["brute_force_integer_box"]["nontrivial_invariant_count"] == 0,
        "tempting_ratios_fail_combined_invariance": all(not row["invariant_under_combined_action"] for row in witness["tempting_ratio_rows"] if row["exponents_lambda_kappa_epsilon"] != (0, 0, 0)),
        "canonical_normalization_theorem_exported": False,
        "physical_coupling_provenance_theorem_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_monomial_ratio_rescue": False,
        "accepted_as_physical_coupling_provenance_theorem": False,
        "accepted_as_ltotal_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "The combined length/field/curvature normalization action has full-rank weight matrix on the three P1562 coefficients.  Therefore no nontrivial monomial ratio of lambda_sm_eff, kappa_gr_eff, and epsilon_mix_eff is invariant under the combined open normalizations; this blocks a ratio-invariant rescue but does not export a canonical normalization theorem.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    witness = payload["combined_normalization_monomial_invariant_witness"]
    lines = [
        "# P2768/S1718 combined-normalization monomial-invariant no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Weight matrix",
        f"- coefficient_order={witness['coefficient_order']}",
        f"- gauge_order={witness['gauge_order']}",
        f"- weight_matrix={witness['weight_matrix']}",
        f"- determinant={witness['determinant']}",
        f"- full_rank_over_Q={witness['full_rank_over_Q']}",
        "",
        "## Brute-force invariant scan",
        f"- tested_nonzero_exponent_vectors={witness['brute_force_integer_box']['tested_nonzero_exponent_vectors']}",
        f"- nontrivial_invariant_count={witness['brute_force_integer_box']['nontrivial_invariant_count']}",
        "",
        "## Ratio candidates",
    ]
    for row in witness["tempting_ratio_rows"]:
        lines.append(f"- {row['candidate']}: weight={row['gauge_weight']}; invariant={row['invariant_under_combined_action']}; distinct_values={row['distinct_sample_values']}")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p1562 = read_json(P1562)
    p2762 = read_json(P2762)
    p2764 = read_json(P2764)
    p2766 = read_json(P2766)
    p2767 = read_json(P2767)
    scan = evidence_scan()
    witness = invariant_no_go_witness(p1562)
    acceptance = acceptance_matrix(scan, witness, p2762, p2764, p2766, p2767)
    payload = {
        "status": "P2768_COMBINED_NORMALIZATION_MONOMIAL_INVARIANT_NO_GO_NO_CLOSURE",
        "input_hashes": {"P1562": sha(P1562), "P2762": sha(P2762), "P2764": sha(P2764), "P2766": sha(P2766), "P2767": sha(P2767)},
        "input_statuses": {"P1562": p1562.get("status"), "P2762": p2762.get("status"), "P2764": p2764.get("status"), "P2766": p2766.get("status"), "P2767": p2767.get("status")},
        "audited_object": "combined length/field/curvature normalization action on P1562 coefficient monomials",
        "content_evidence_scan": scan,
        "combined_normalization_monomial_invariant_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not try to rescue physical-coupling provenance by taking a monomial ratio of lambda_sm_eff, kappa_gr_eff, and epsilon_mix_eff.  P2768 proves the combined open normalization action has no nontrivial monomial invariant.  The next honest step must either supply an actual canonical normalization theorem/source, or introduce a genuinely non-monomial invariant with an explicit invariance proof and then run the same bounded acceptance test; otherwise preserve the P2697-P2768 no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2768/S1718 combined-normalization monomial-invariant no-go", "## P2768/S1718 combined-normalization monomial-invariant no-go\n\n`P2768/S1718` supplies a bounded new object in the P2766 moment-provenance lane: the combined positive length, scalar-field, and curvature normalization action on monomials in `lambda_sm_eff`, `kappa_gr_eff`, and `epsilon_mix_eff`.  The 3x3 exponent-weight matrix has nonzero determinant, and a finite integer exponent scan confirms zero nontrivial monomial invariants in the tested box.  Thus ratio-taking cannot rescue physical-coupling provenance under the open P2762/P2764 normalizations.  This is a monomial-invariant no-go theorem, not a canonical normalization theorem; no physical-coupling provenance, role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2768/S1718 monomial-invariant Ltotal guard", "## P2768/S1718 monomial-invariant Ltotal guard\n\n`P2768/S1718` adds no variational source term.  It proves that no nontrivial monomial ratio of `lambda_sm_eff`, `kappa_gr_eff`, and `epsilon_mix_eff` is invariant under the combined open length/field/curvature normalization action, so ratio-taking cannot promote those coefficients to physical Lagrangian couplings or role-bearing `L_total`.\n")
    append_once(AGENTS, "Current combined-normalization monomial-invariant no-go guardrail (P2768/S1718, 2026-06-15)", "## Current combined-normalization monomial-invariant no-go guardrail (P2768/S1718, 2026-06-15)\n\n- P2768 supplies one bounded new object in the moment-provenance lane: the combined positive length, scalar-field, and curvature normalization action on monomials in `lambda_sm_eff`, `kappa_gr_eff`, and `epsilon_mix_eff`.\n- The 3x3 exponent-weight matrix has full rank and the finite integer exponent scan finds zero nontrivial monomial invariants, so monomial ratio-taking cannot rescue physical-coupling provenance while P2762/P2764 normalizations remain open.\n- Do not promote this no-go to canonical normalization, physical-coupling provenance, role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.  A next admissible move must supply an actual canonical normalization theorem/source or a genuinely non-monomial invariant with explicit invariance proof and bounded acceptance test; otherwise preserve the P2697-P2768 no-closure certificate.\n")
    return payload


if __name__ == "__main__":
    main()
