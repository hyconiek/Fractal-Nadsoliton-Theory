#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2565_s1515_phase_frequency_selector_objective_grid_audit.json"
MD = GEN / "p2565_s1515_phase_frequency_selector_objective_grid_audit.md"

SOURCE_FILES = {
    "P2564_SIGN_CELL_NONIDENTIFIABILITY": GEN / "p2564_s1514_phase_frequency_finite_sign_cell_nonidentifiability_certificate.json",
    "P2563_RATIONAL_WINDING_OBSTRUCTION": GEN / "p2563_s1513_phase_frequency_rational_winding_quotient_obstruction_certificate.json",
    "P2561_POST_DAMPING_RESIDUAL_BRIDGE": GEN / "p2561_s1511_strict_completion_post_damping_residual_bridge_two_key_certificate.json",
}

STRICT_OMEGA = 743.0 / 4000.0
STRICT_PHI = 13.0 / 80.0
DOMAIN = list(range(12))
GRID_STEPS_PER_AXIS = 21
NEGATIVE_EXPORT_FLAGS = [
    "selector_objective_source_exported",
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
        "new_packet": "P2565|S1515|phase frequency selector objective|selector objective grid|open sign cell objective|phase/frequency objective grid",
        "intended_research_nonduplication": "selector objective.*phase|phase.*selector objective|objective.*sign cell|sign-cell.*objective|max.*clearance|signed cosine margin|phase L2 objective",
        "immediate_precursors": "P2564|S1514|finite sign-cell nonidentifiability|P2563|S1513|P2561|S1511|strict_phase_frequency_source",
        "selector_guardrails": "QW-2191|selector closure|symmetry-breaking|source theorem|non-strict|ToE closure",
        "toe_context": "ToE|Theory of Everything|role-transfer theorem|legacy -> strict completion bridge|role-bearing L_total",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def sign(value: float) -> int:
    if value > 0.0:
        return 1
    if value < 0.0:
        return -1
    return 0


def zero_clearance(omega: float, phi: float) -> float:
    return min(abs(omega * d + phi - (math.pi / 2.0 + k * math.pi)) for d in DOMAIN for k in range(-8, 9))


def signs_for(omega: float, phi: float) -> list[int]:
    return [sign(math.cos(omega * d + phi)) for d in DOMAIN]


def objective_values(omega: float, phi: float, strict_signs: list[int]) -> dict[str, float]:
    phases = [omega * d + phi for d in DOMAIN]
    signed_cos = [strict_signs[d] * math.cos(phases[d]) for d in DOMAIN]
    return {
        "max_min_zero_clearance": zero_clearance(omega, phi),
        "max_signed_cosine_sum": sum(signed_cos),
        "max_signed_cosine_min_margin": min(signed_cos),
        "min_phase_l2": -sum(theta * theta for theta in phases),
        "min_parameter_l2": -(omega * omega + phi * phi),
    }


def scan_objectives(box: dict[str, Any], strict_signs: list[int]) -> dict[str, Any]:
    eps_omega = box["epsilon_omega"]
    eps_phi = box["epsilon_phi"]
    strict_objectives = objective_values(STRICT_OMEGA, STRICT_PHI, strict_signs)
    winners: dict[str, dict[str, Any]] = {}
    accepted_rows: list[dict[str, Any]] = []
    for i in range(GRID_STEPS_PER_AXIS):
        omega_factor = -1.0 + 2.0 * i / (GRID_STEPS_PER_AXIS - 1)
        for j in range(GRID_STEPS_PER_AXIS):
            phi_factor = -1.0 + 2.0 * j / (GRID_STEPS_PER_AXIS - 1)
            delta_omega = omega_factor * eps_omega
            delta_phi = phi_factor * eps_phi
            omega = STRICT_OMEGA + delta_omega
            phi = STRICT_PHI + delta_phi
            signs = signs_for(omega, phi)
            if signs != strict_signs:
                continue
            values = objective_values(omega, phi, strict_signs)
            row = {
                "omega_factor": omega_factor,
                "phi_factor": phi_factor,
                "delta_omega": delta_omega,
                "delta_phi": delta_phi,
                "omega": omega,
                "phi": phi,
                "objectives": values,
                "is_strict_tuple": abs(delta_omega) < 1e-15 and abs(delta_phi) < 1e-15,
            }
            accepted_rows.append(row)
            for name, value in values.items():
                if name not in winners or value > winners[name]["objective_value"]:
                    winners[name] = {**{key: row[key] for key in ["omega_factor", "phi_factor", "delta_omega", "delta_phi", "omega", "phi", "is_strict_tuple"]}, "objective_value": value, "strict_objective_value": strict_objectives[name], "improvement_over_strict": value - strict_objectives[name]}
    return {
        "grid_steps_per_axis": GRID_STEPS_PER_AXIS,
        "candidate_grid_point_count": GRID_STEPS_PER_AXIS * GRID_STEPS_PER_AXIS,
        "accepted_sign_preserving_point_count": len(accepted_rows),
        "strict_objective_values": strict_objectives,
        "objective_winners": winners,
        "objective_names": sorted(winners),
        "objective_count": len(winners),
        "objectives_selecting_strict_tuple_count": sum(1 for winner in winners.values() if winner["is_strict_tuple"]),
        "objectives_with_non_strict_winner_count": sum(1 for winner in winners.values() if not winner["is_strict_tuple"]),
        "sample_non_strict_winners": {name: winner for name, winner in winners.items() if not winner["is_strict_tuple"]},
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2564_payload = load_json(SOURCE_FILES["P2564_SIGN_CELL_NONIDENTIFIABILITY"])
    p2563_payload = load_json(SOURCE_FILES["P2563_RATIONAL_WINDING_OBSTRUCTION"])
    p2561_payload = load_json(SOURCE_FILES["P2561_POST_DAMPING_RESIDUAL_BRIDGE"])
    p2564 = theorem(p2564_payload, "phase_frequency_finite_sign_cell_nonidentifiability_certificate")
    p2563 = theorem(p2563_payload, "phase_frequency_rational_winding_quotient_obstruction_certificate")
    p2561 = theorem(p2561_payload, "strict_completion_post_damping_residual_bridge_two_key_certificate")
    strict_signs = p2564.get("strict_sign_pattern_d0_to_d11", signs_for(STRICT_OMEGA, STRICT_PHI))
    scan = scan_objectives(p2564["certified_open_sign_cell_box"], strict_signs)
    toe_potential = {
        "positive_potential": [
            "The strict kernel has a compact numerical tuple and coherent finite sign/topological certificates.",
            "P2502/P2561 isolate a small bridge frontier rather than an unstructured list of missing facts.",
            "P2564/P2565 turn phase/frequency ambiguity into a computable selector problem inside a certified open cell.",
        ],
        "current_limit": "This is ToE potential only: the bridge atoms, selector source, role-transfer theorem, and QW-2191 discharge are still not exported.",
    }
    closure_problem_diagnosis = {
        "why_closure_is_hard": [
            "Endpoint and finite sign data are underdetermined: open continua of kernels preserve the audited facts.",
            "Natural source-free objectives are not canonical; different objectives select different non-strict points in the same cell.",
            "The legacy-to-strict bridge still needs independent strict sources for phase/frequency and A/P/D dynamics after damping progress.",
            "QW-2191 blocks uniqueness/selector closure unless a real symmetry-breaking source is added.",
            "Even a bridge theorem would still require a separate role-transfer audit before legacy physical claims can become ToE claims.",
        ],
        "honest_status": "The present repository has strong finite evidence and good obstruction localization, but not a closed ToE theorem.",
    }
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2565_T1_phase_frequency_selector_objective_grid_audit",
        "audited_chain": ["P2564/S1514", "P2563/S1513", "P2561/S1511"],
        "frontier_atom_under_attack": "strict_phase_frequency_source",
        "p2564_open_sign_cell_inherited": p2564.get("finite_sign_pattern_has_open_continuum_of_phase_frequency_realizations") is True,
        "p2563_rational_winding_obstruction_inherited": p2563.get("exact_rational_winding_quotient_source_rejected") is True,
        "p2561_phase_frequency_residual_atom_inherited": "strict_phase_frequency_source" in p2561.get("residual_atoms_after_hypothetical_damping_source", []),
        "selector_objective_grid_audit": scan,
        "non_circular_objective_candidates_select_strict_tuple": scan["objectives_selecting_strict_tuple_count"] > 0,
        "all_audited_objectives_have_non_strict_winners": scan["objectives_with_non_strict_winner_count"] == scan["objective_count"],
        "source_free_objective_choice_remains_extra_obligation": True,
        "toe_potential_assessment": toe_potential,
        "closure_problem_diagnosis": closure_problem_diagnosis,
        "recommended_next_honest_step": (
            "Do not add more source-free objective scans as if objective choice were automatic. The next honest step is to derive one selector functional from strict nadsoliton dynamics, then rerun this grid/interval audit to test whether the derived functional uniquely selects omega=743/4000 and phi=13/80 while keeping QW-2191 open unless the derivation supplies real symmetry breaking."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2564_open_cell_inherited": theorem_export["p2564_open_sign_cell_inherited"],
        "p2563_inherited": theorem_export["p2563_rational_winding_obstruction_inherited"],
        "p2561_phase_frequency_atom_inherited": theorem_export["p2561_phase_frequency_residual_atom_inherited"],
        "five_objectives_audited": scan["objective_count"] == 5,
        "grid_has_441_points": scan["candidate_grid_point_count"] == 441,
        "all_objectives_have_non_strict_winners": theorem_export["all_audited_objectives_have_non_strict_winners"],
        "objective_choice_obligation_recorded": theorem_export["source_free_objective_choice_remains_extra_obligation"],
        "no_phase_source_exported": theorem_export["strict_phase_frequency_source_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2565",
        "stage_id": "S1515",
        "status": "P2565_PHASE_FREQUENCY_SELECTOR_OBJECTIVE_GRID_AUDIT_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "phase_frequency_selector_objective_grid_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2564_SIGN_CELL_NONIDENTIFIABILITY": sha256_json(p2564_payload),
                "P2563_RATIONAL_WINDING_OBSTRUCTION": sha256_json(p2563_payload),
                "P2561_POST_DAMPING_RESIDUAL_BRIDGE": sha256_json(p2561_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["phase_frequency_selector_objective_grid_audit"]["theorem_export"]
    scan = t["selector_objective_grid_audit"]
    lines = [
        "# P2565/S1515 phase/frequency selector-objective grid audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2564 open sign cell inherited: `{t['p2564_open_sign_cell_inherited']}`.",
        f"- P2563 rational-winding obstruction inherited: `{t['p2563_rational_winding_obstruction_inherited']}`.",
        f"- Grid points audited: `{scan['candidate_grid_point_count']}`.",
        f"- Sign-preserving grid points: `{scan['accepted_sign_preserving_point_count']}`.",
        f"- Objective candidates audited: `{scan['objective_count']}`.",
        f"- Objectives with non-strict winners: `{scan['objectives_with_non_strict_winner_count']}`.",
        f"- Source-free objective choice remains extra obligation: `{t['source_free_objective_choice_remains_extra_obligation']}`.", "",
        "## ToE potential", "",
        "The program has real ToE potential because the bridge frontier is now small, computable, and obstruction-localized.  But this packet does not claim ToE closure: it shows that phase/frequency selector choice is still an extra source obligation.", "",
        "## Why closure is hard", "",
        "Endpoint, sign, and finite-topological data underdetermine the exact kernel tuple; natural objective choices are not canonical; QW-2191 and role-transfer remain separate gates.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No selector-objective source, strict phase/frequency source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['phase_frequency_selector_objective_grid_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2565/S1515` audits five non-circular, source-free selector-objective candidates inside the P2564 open phase-sign cell on a `21 x 21` grid.  The audited objectives include zero-clearance margin, signed-cosine aggregate, signed-cosine minimum margin, phase `L2`, and parameter `L2`.  Every audited objective has a non-strict grid winner, so the exact tuple `omega=743/4000`, `phi=13/80` is still not selected without a sourced objective principle.  This sharpens the ToE bottleneck: the route has computable potential, but source selection remains underived.
""".strip()
    lag_section = """
`P2565/S1515` blocks promotion of source-free selector objectives into role-bearing phase/frequency dynamics in `L_total`.  The phase/frequency source must come from strict nadsoliton dynamics or a real selector premise; QW-2191, role-transfer, and ToE closure remain open.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2565/S1515 phase/frequency selector-objective grid guard", "## P2565/S1515 phase/frequency selector-objective grid guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2565/S1515 phase/frequency selector-objective grid Ltotal guard", "## P2565/S1515 phase/frequency selector-objective grid Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
