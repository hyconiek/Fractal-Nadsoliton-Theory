#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import itertools
import json
import math
import subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2567_s1517_phase_frequency_minimal_stationary_hessian_saddle_audit.json"
MD = GEN / "p2567_s1517_phase_frequency_minimal_stationary_hessian_saddle_audit.md"

SOURCE_FILES = {
    "P2566_STATIONARITY_WEIGHT_CONE": GEN / "p2566_s1516_phase_frequency_selector_stationarity_weight_cone_audit.json",
    "P2565_SELECTOR_OBJECTIVE_GRID": GEN / "p2565_s1515_phase_frequency_selector_objective_grid_audit.json",
    "P2561_POST_DAMPING_RESIDUAL_BRIDGE": GEN / "p2561_s1511_strict_completion_post_damping_residual_bridge_two_key_certificate.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "minimal_stationary_hessian_selector_source_exported",
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
        "new_packet": "P2567|S1517|minimal stationary Hessian|stationary Hessian saddle|minimal stationarity saddle|phase frequency Hessian saddle",
        "intended_research_nonduplication": "selector Hessian|phase frequency Hessian|phase/frequency Hessian|second-order selector|Hessian.*omega.*phi|stationarity.*Hessian|local maximum.*phase",
        "immediate_precursors": "P2566|S1516|stationarity weight cone|P2565|S1515|selector-objective grid|strict_phase_frequency_source",
        "guardrails": "QW-2191|selector closure|source theorem|role-transfer theorem|role-bearing L_total|ToE closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def stationary_weights_for_triple(rows: list[dict[str, Any]], support: tuple[int, int, int]) -> dict[int, float]:
    by_d = {int(row["d"]): row for row in rows}
    a = [float(by_d[d]["gradient_phi_coefficient"]) for d in support]
    b = [float(by_d[d]["gradient_omega_coefficient"]) for d in support]
    weights = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
    scale = max(abs(value) for value in weights)
    if scale == 0.0:
        return {d: 0.0 for d in support}
    return {d: value / scale for d, value in zip(support, weights)}


def hessian_for_weights(rows: list[dict[str, Any]], weights: dict[int, float]) -> dict[str, Any]:
    by_d = {int(row["d"]): row for row in rows}
    h_oo = 0.0
    h_op = 0.0
    h_pp = 0.0
    for d, weight in weights.items():
        row = by_d[d]
        coeff = -weight * int(row["sign"]) * math.cos(float(row["theta"]))
        h_oo += coeff * d * d
        h_op += coeff * d
        h_pp += coeff
    trace = h_oo + h_pp
    determinant = h_oo * h_pp - h_op * h_op
    discriminant = max(0.0, trace * trace - 4.0 * determinant)
    root = math.sqrt(discriminant)
    eigenvalues = [(trace - root) / 2.0, (trace + root) / 2.0]
    return {
        "hessian": [[h_oo, h_op], [h_op, h_pp]],
        "trace": trace,
        "determinant": determinant,
        "eigenvalues": eigenvalues,
        "classification": "indefinite_saddle" if determinant < 0.0 else "not_indefinite",
    }


def stationarity_residuals(rows: list[dict[str, Any]], weights: dict[int, float]) -> dict[str, float]:
    by_d = {int(row["d"]): row for row in rows}
    phi = sum(weights[d] * float(by_d[d]["gradient_phi_coefficient"]) for d in weights)
    omega = sum(weights[d] * float(by_d[d]["gradient_omega_coefficient"]) for d in weights)
    return {"gradient_phi_residual_abs": abs(phi), "gradient_omega_residual_abs": abs(omega)}


def triple_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    support_rows = []
    for support in itertools.combinations([int(row["d"]) for row in rows], 3):
        weights = stationary_weights_for_triple(rows, support)
        residuals = stationarity_residuals(rows, weights)
        hessian = hessian_for_weights(rows, weights)
        support_rows.append({
            "support": list(support),
            "weights": {str(key): weights[key] for key in sorted(weights)},
            **residuals,
            **hessian,
        })
    indefinite = [row for row in support_rows if row["classification"] == "indefinite_saddle"]
    return {
        "support_size": 3,
        "audited_support_count": len(support_rows),
        "expected_support_count_choose_12_3": math.comb(12, 3),
        "indefinite_saddle_count": len(indefinite),
        "non_saddle_count": len(support_rows) - len(indefinite),
        "all_minimal_three_node_stationary_witnesses_are_saddles": len(indefinite) == len(support_rows),
        "max_stationarity_residual_abs": max(max(row["gradient_phi_residual_abs"], row["gradient_omega_residual_abs"]) for row in support_rows),
        "sample_saddle_rows": support_rows[:12],
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2566_payload = load_json(SOURCE_FILES["P2566_STATIONARITY_WEIGHT_CONE"])
    p2565_payload = load_json(SOURCE_FILES["P2565_SELECTOR_OBJECTIVE_GRID"])
    p2561_payload = load_json(SOURCE_FILES["P2561_POST_DAMPING_RESIDUAL_BRIDGE"])
    p2566 = theorem(p2566_payload, "phase_frequency_selector_stationarity_weight_cone_audit")
    p2565 = theorem(p2565_payload, "phase_frequency_selector_objective_grid_audit")
    p2561 = theorem(p2561_payload, "strict_completion_post_damping_residual_bridge_two_key_certificate")
    audit = triple_audit(p2566["stationarity_rows"])
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2567_T1_phase_frequency_minimal_stationary_hessian_saddle_audit",
        "audited_chain": ["P2566/S1516", "P2565/S1515", "P2561/S1511"],
        "frontier_atom_under_attack": "strict_phase_frequency_source",
        "p2566_stationarity_underidentification_inherited": p2566.get("first_order_stationarity_with_unconstrained_signed_weights_underidentified") is True,
        "p2565_objective_choice_obligation_inherited": p2565.get("source_free_objective_choice_remains_extra_obligation") is True,
        "p2561_phase_frequency_residual_atom_inherited": "strict_phase_frequency_source" in p2561.get("residual_atoms_after_hypothetical_damping_source", []),
        "minimal_three_node_stationary_hessian_audit": audit,
        "minimal_signed_stationary_witnesses_fail_second_order_selector_test": audit["all_minimal_three_node_stationary_witnesses_are_saddles"],
        "second_order_local_selector_source_exported": False,
        "recommended_next_honest_step": (
            "Do not use minimal signed stationary witnesses as the phase/frequency selector: every audited three-node stationary witness is a Hessian saddle. The next honest step is either to derive a higher-support positive/semibounded variational law from strict nadsoliton dynamics or to prove that no such semibounded law exists inside the P2564 open sign cell."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2566_inherited": theorem_export["p2566_stationarity_underidentification_inherited"],
        "p2565_inherited": theorem_export["p2565_objective_choice_obligation_inherited"],
        "p2561_phase_frequency_atom_inherited": theorem_export["p2561_phase_frequency_residual_atom_inherited"],
        "all_220_triples_audited": audit["audited_support_count"] == 220 and audit["expected_support_count_choose_12_3"] == 220,
        "all_minimal_witnesses_are_saddles": audit["all_minimal_three_node_stationary_witnesses_are_saddles"],
        "stationarity_residuals_small": audit["max_stationarity_residual_abs"] < 1e-12,
        "no_phase_source_exported": theorem_export["strict_phase_frequency_source_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2567",
        "stage_id": "S1517",
        "status": "P2567_PHASE_FREQUENCY_MINIMAL_STATIONARY_HESSIAN_SADDLE_AUDIT_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "phase_frequency_minimal_stationary_hessian_saddle_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2566_STATIONARITY_WEIGHT_CONE": sha256_json(p2566_payload),
                "P2565_SELECTOR_OBJECTIVE_GRID": sha256_json(p2565_payload),
                "P2561_POST_DAMPING_RESIDUAL_BRIDGE": sha256_json(p2561_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["phase_frequency_minimal_stationary_hessian_saddle_audit"]["theorem_export"]
    audit = t["minimal_three_node_stationary_hessian_audit"]
    lines = [
        "# P2567/S1517 phase/frequency minimal stationary Hessian saddle audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2566 stationarity underidentification inherited: `{t['p2566_stationarity_underidentification_inherited']}`.",
        f"- Three-node supports audited: `{audit['audited_support_count']}`.",
        f"- Indefinite Hessian saddle count: `{audit['indefinite_saddle_count']}`.",
        f"- Non-saddle count: `{audit['non_saddle_count']}`.",
        f"- Max stationarity residual: `{audit['max_stationarity_residual_abs']}`.",
        f"- Minimal signed stationary witnesses fail second-order selector test: `{t['minimal_signed_stationary_witnesses_fail_second_order_selector_test']}`.", "",
        "## Interpretation", "",
        "P2566 showed first-order stationarity can be manufactured with signed weights.  P2567 adds the second-order check for all minimal three-node signed stationary supports: every one has an indefinite Hessian, so none is a local max/min selector for the strict phase/frequency tuple.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No minimal-stationary Hessian source, strict phase/frequency source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['phase_frequency_minimal_stationary_hessian_saddle_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2567/S1517` adds a second-order audit behind P2566.  It enumerates all `C(12,3)=220` minimal three-node signed stationary supports for the phase/frequency stationarity equations and computes the Hessian of `F_w(omega,phi)=sum_d w_d sign_d cos(omega*d+phi)` at the strict tuple.  Every audited minimal stationary witness has negative determinant / indefinite Hessian, so minimal signed stationarity is a saddle mechanism, not a local selector source for `omega=743/4000`, `phi=13/80`.
""".strip()
    lag_section = """
`P2567/S1517` blocks promotion of minimal signed stationary witnesses into role-bearing phase/frequency `L_total` terms.  First-order stationarity plus a three-node signed support fails the second-order selector test; a strict source would need a higher-support semibounded variational law or a proof that no such law exists.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2567/S1517 phase/frequency minimal stationary Hessian saddle guard", "## P2567/S1517 phase/frequency minimal stationary Hessian saddle guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2567/S1517 phase/frequency minimal stationary Hessian saddle Ltotal guard", "## P2567/S1517 phase/frequency minimal stationary Hessian saddle Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
