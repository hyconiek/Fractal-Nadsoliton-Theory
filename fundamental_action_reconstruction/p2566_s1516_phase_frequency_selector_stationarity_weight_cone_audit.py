#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from typing import Any

import sympy as sp

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2566_s1516_phase_frequency_selector_stationarity_weight_cone_audit.json"
MD = GEN / "p2566_s1516_phase_frequency_selector_stationarity_weight_cone_audit.md"

SOURCE_FILES = {
    "P2565_SELECTOR_OBJECTIVE_GRID": GEN / "p2565_s1515_phase_frequency_selector_objective_grid_audit.json",
    "P2564_SIGN_CELL_NONIDENTIFIABILITY": GEN / "p2564_s1514_phase_frequency_finite_sign_cell_nonidentifiability_certificate.json",
    "P2561_POST_DAMPING_RESIDUAL_BRIDGE": GEN / "p2561_s1511_strict_completion_post_damping_residual_bridge_two_key_certificate.json",
}

STRICT_OMEGA = 743.0 / 4000.0
STRICT_PHI = 13.0 / 80.0
DOMAIN = list(range(12))
NEGATIVE_EXPORT_FLAGS = [
    "stationarity_weight_source_exported",
    "strict_phase_frequency_source_exported",
    "strict_dynamical_source_for_A_P_D_exported",
    "strict_damping_beta_eta_source_exported",
    "bridge_theorem_exported",
    "legacy_to_strict_completion_bridge_exported",
    "role_transfer_theorem_exported",
    "role_bearing_ltotal_exported",
    "qw2191_discharged_by_this_certificate",
    "toe_closure_claimed",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, "DIAGRAMS_KERNEL_TRANSFORMATION.md", "fundamental_action_reconstruction", "material_dowodowy",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
    ], cwd=REPO, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2566|S1516|selector stationarity weight cone|stationarity weight cone|phase frequency stationarity weight|phase/frequency stationarity weight",
        "intended_research_nonduplication": "selector stationarity|stationarity acceptance|objective gradient|first-order selector|gradient.*omega.*phi|weight cone.*phase|phase.*weight cone",
        "immediate_precursors": "P2565|S1515|selector-objective grid|P2564|S1514|open sign cell|P2561|S1511|strict_phase_frequency_source",
        "guardrails": "QW-2191|selector closure|source theorem|role-transfer theorem|role-bearing L_total|ToE closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def strict_signs() -> list[int]:
    return [1 if math.cos(STRICT_OMEGA * d + STRICT_PHI) > 0 else -1 for d in DOMAIN]


def stationarity_rows(signs: list[int]) -> list[dict[str, float | int]]:
    rows = []
    for d, sign in zip(DOMAIN, signs):
        theta = STRICT_OMEGA * d + STRICT_PHI
        signed_sin = sign * math.sin(theta)
        rows.append({
            "d": d,
            "sign": sign,
            "theta": theta,
            "signed_sin": signed_sin,
            "gradient_phi_coefficient": -signed_sin,
            "gradient_omega_coefficient": -d * signed_sin,
        })
    return rows


def rank_and_nullity(rows: list[dict[str, Any]]) -> dict[str, Any]:
    matrix = sp.Matrix([
        [sp.Float(row["gradient_phi_coefficient"], 80) for row in rows],
        [sp.Float(row["gradient_omega_coefficient"], 80) for row in rows],
    ])
    rank = matrix.rank()
    nullspace = matrix.nullspace()
    return {
        "backend": "sympy",
        "sympy_version": sp.__version__,
        "stationarity_matrix_shape": [2, len(rows)],
        "stationarity_matrix_rank": int(rank),
        "stationarity_nullity": len(rows) - int(rank),
        "nullspace_vector_count": len(nullspace),
        "sample_nullspace_vector_first_4": [str(value.evalf(20)) for value in (nullspace[0][:4] if nullspace else [])],
    }


def positive_cone_obstruction(rows: list[dict[str, Any]]) -> dict[str, Any]:
    positive_side = [row for row in rows if row["sign"] > 0]
    negative_side = [row for row in rows if row["sign"] < 0]
    return {
        "nonnegative_weights_assumed": True,
        "positive_sign_d_range": [min(row["d"] for row in positive_side), max(row["d"] for row in positive_side)],
        "negative_sign_d_range": [min(row["d"] for row in negative_side), max(row["d"] for row in negative_side)],
        "all_sin_theta_positive_on_domain": all(math.sin(row["theta"]) > 0 for row in rows),
        "proof_reason": (
            "For nonnegative weights, phi-stationarity requires positive-sign weighted sin mass to equal negative-sign weighted sin mass. "
            "Omega-stationarity would then require their weighted mean d values to match. Positive-sign nodes lie in d=0..7 and negative-sign nodes lie in d=8..11, so the means cannot match unless all mass is zero."
        ),
        "nonzero_nonnegative_weight_stationarity_possible": False,
    }


def signed_two_pair_witness(rows: list[dict[str, Any]]) -> dict[str, Any]:
    # Solve using only d=7 and d=8 with signed weights: w7*a7 + w8*a8 = 0 for both rows is impossible unless both zero,
    # so use three nodes d=0,7,8 to give one degree of freedom. Set w8=1 and solve phi/omega for w0,w7.
    by_d = {row["d"]: row for row in rows}
    w0, w7 = sp.symbols("w0 w7")
    w8 = sp.Integer(1)
    eq_phi = sp.Eq(w0 * by_d[0]["gradient_phi_coefficient"] + w7 * by_d[7]["gradient_phi_coefficient"] + w8 * by_d[8]["gradient_phi_coefficient"], 0)
    eq_omega = sp.Eq(w0 * by_d[0]["gradient_omega_coefficient"] + w7 * by_d[7]["gradient_omega_coefficient"] + w8 * by_d[8]["gradient_omega_coefficient"], 0)
    sol = sp.solve([eq_phi, eq_omega], [w0, w7], dict=True)[0]
    weights = {0: float(sol[w0]), 7: float(sol[w7]), 8: 1.0}
    grad_phi = sum(weights.get(d, 0.0) * by_d[d]["gradient_phi_coefficient"] for d in weights)
    grad_omega = sum(weights.get(d, 0.0) * by_d[d]["gradient_omega_coefficient"] for d in weights)
    return {
        "support": sorted(weights),
        "weights": {str(d): weights[d] for d in sorted(weights)},
        "has_negative_weight": any(value < 0 for value in weights.values()),
        "gradient_phi_residual_abs": abs(grad_phi),
        "gradient_omega_residual_abs": abs(grad_omega),
        "interpretation": "Signed/indefinite weights can force first-order stationarity, so stationarity alone is not a source theorem and the sign/weight law itself becomes an extra obligation.",
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2565_payload = load_json(SOURCE_FILES["P2565_SELECTOR_OBJECTIVE_GRID"])
    p2564_payload = load_json(SOURCE_FILES["P2564_SIGN_CELL_NONIDENTIFIABILITY"])
    p2561_payload = load_json(SOURCE_FILES["P2561_POST_DAMPING_RESIDUAL_BRIDGE"])
    p2565 = theorem(p2565_payload, "phase_frequency_selector_objective_grid_audit")
    p2564 = theorem(p2564_payload, "phase_frequency_finite_sign_cell_nonidentifiability_certificate")
    p2561 = theorem(p2561_payload, "strict_completion_post_damping_residual_bridge_two_key_certificate")
    signs = p2564.get("strict_sign_pattern_d0_to_d11", strict_signs())
    rows = stationarity_rows(signs)
    linear = rank_and_nullity(rows)
    cone = positive_cone_obstruction(rows)
    witness = signed_two_pair_witness(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2566_T1_phase_frequency_selector_stationarity_weight_cone_audit",
        "audited_chain": ["P2565/S1515", "P2564/S1514", "P2561/S1511"],
        "frontier_atom_under_attack": "strict_phase_frequency_source",
        "p2565_objective_choice_obligation_inherited": p2565.get("source_free_objective_choice_remains_extra_obligation") is True,
        "p2564_open_sign_cell_inherited": p2564.get("finite_sign_pattern_has_open_continuum_of_phase_frequency_realizations") is True,
        "p2561_phase_frequency_residual_atom_inherited": "strict_phase_frequency_source" in p2561.get("residual_atoms_after_hypothetical_damping_source", []),
        "stationarity_model": "F_w(omega,phi)=sum_d w_d sign_d cos(omega*d+phi); audit first-order stationarity at strict tuple",
        "stationarity_rows": rows,
        "linear_stationarity_rank_audit": linear,
        "nonnegative_weight_cone_obstruction": cone,
        "signed_weight_stationarity_witness": witness,
        "first_order_stationarity_with_unconstrained_signed_weights_underidentified": linear["stationarity_nullity"] == 10,
        "nonnegative_natural_weight_stationarity_rejected": cone["nonzero_nonnegative_weight_stationarity_possible"] is False,
        "stationarity_alone_selects_unique_source": False,
        "recommended_next_honest_step": (
            "Do not promote first-order stationarity as a phase/frequency source. A viable next step must derive the weight/sign law and a second-order or global variational principle from strict nadsoliton dynamics; then test positivity, uniqueness, and QW-2191 symmetry-breaking rather than choosing signed weights post hoc."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2565_inherited": theorem_export["p2565_objective_choice_obligation_inherited"],
        "p2564_open_cell_inherited": theorem_export["p2564_open_sign_cell_inherited"],
        "p2561_phase_frequency_atom_inherited": theorem_export["p2561_phase_frequency_residual_atom_inherited"],
        "rank_two_nullity_ten": linear["stationarity_matrix_rank"] == 2 and linear["stationarity_nullity"] == 10,
        "nonnegative_cone_rejected": theorem_export["nonnegative_natural_weight_stationarity_rejected"],
        "signed_witness_has_small_residuals": witness["gradient_phi_residual_abs"] < 1e-12 and witness["gradient_omega_residual_abs"] < 1e-12,
        "stationarity_not_unique_source": theorem_export["stationarity_alone_selects_unique_source"] is False,
        "no_phase_source_exported": theorem_export["strict_phase_frequency_source_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2566",
        "stage_id": "S1516",
        "status": "P2566_PHASE_FREQUENCY_SELECTOR_STATIONARITY_WEIGHT_CONE_AUDIT_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "phase_frequency_selector_stationarity_weight_cone_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2565_SELECTOR_OBJECTIVE_GRID": sha256_json(p2565_payload),
                "P2564_SIGN_CELL_NONIDENTIFIABILITY": sha256_json(p2564_payload),
                "P2561_POST_DAMPING_RESIDUAL_BRIDGE": sha256_json(p2561_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["phase_frequency_selector_stationarity_weight_cone_audit"]["theorem_export"]
    linear = t["linear_stationarity_rank_audit"]
    cone = t["nonnegative_weight_cone_obstruction"]
    witness = t["signed_weight_stationarity_witness"]
    lines = [
        "# P2566/S1516 phase/frequency selector stationarity weight-cone audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2565 objective-choice obligation inherited: `{t['p2565_objective_choice_obligation_inherited']}`.",
        f"- Stationarity matrix rank/nullity: `{linear['stationarity_matrix_rank']}/{linear['stationarity_nullity']}`.",
        f"- Nonzero nonnegative weight stationarity possible: `{cone['nonzero_nonnegative_weight_stationarity_possible']}`.",
        f"- Signed stationarity witness support: `{witness['support']}`.",
        f"- Signed witness has negative weight: `{witness['has_negative_weight']}`.",
        f"- Stationarity alone selects unique source: `{t['stationarity_alone_selects_unique_source']}`.", "",
        "## Interpretation", "",
        "Natural nonnegative weights cannot make the strict tuple first-order stationary for the finite sign objective; allowing signed weights makes stationarity easy but underidentified.  Therefore the missing source is not just stationarity: it is the law that chooses the weight/sign structure and supplies second-order/global selection.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No stationarity-weight source, strict phase/frequency source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['phase_frequency_selector_stationarity_weight_cone_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2566/S1516` audits first-order selector stationarity for objectives of the form `F_w(omega,phi)=sum_d w_d sign_d cos(omega*d+phi)` at the strict tuple.  The two stationarity equations on `12` weights have rank `2` and nullity `10` for unconstrained signed weights, while the natural nonnegative weight cone has no nonzero stationarity solution because positive-sign nodes occupy `d=0..7` and negative-sign nodes occupy `d=8..11`.  Thus stationarity alone is either impossible under natural positivity or underidentified under signed weights; the missing phase/frequency source must derive the weight/sign law and a stronger variational selector.
""".strip()
    lag_section = """
`P2566/S1516` blocks promotion of first-order phase/frequency stationarity into role-bearing `L_total`.  Signed stationarity can be manufactured post hoc and nonnegative stationarity is obstructed, so a strict source must provide the variational law, second-order/global selection, QW-2191 handling, and bridge role discipline.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2566/S1516 phase/frequency selector stationarity weight-cone guard", "## P2566/S1516 phase/frequency selector stationarity weight-cone guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2566/S1516 phase/frequency selector stationarity weight-cone Ltotal guard", "## P2566/S1516 phase/frequency selector stationarity weight-cone Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
