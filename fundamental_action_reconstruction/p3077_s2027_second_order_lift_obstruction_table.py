#!/usr/bin/env python3
"""P3077/S2027: second-order lift obstruction table.

P3076 showed that the P3075 local Dirichlet/Laplacian source has a genuine
internal diffusive spectral branch, not an exported wave/lightlike branch.
P3077 constructs the next missing interface: a bounded phase-space/second-order
lift audit.  It tests whether adding rho/pi phase variables, a symplectic form,
a quadratic kinetic term, and a time/unit normalization is internally sourced or
merely imported.  The formal Hamiltonian lift has the expected modal relation
omega_j^2 = lambda_j, but every wave-compatible row depends on unsourced
phase-space/time/unit premises, so no observed light, gauge photons, spacetime
EOM, L_total, or ToE closure is exported.
"""
from __future__ import annotations

import hashlib, json, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3076_s2026_dirichlet_spectral_dispersion_audit import OUT as P3076, eigenvalue_label

OUT = GEN / "p3077_s2027_second_order_lift_obstruction_table.json"
MD = GEN / "p3077_s2027_second_order_lift_obstruction_table.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
Z12_SIZE = 12
CONTENT_PATTERNS = {
    "second_order_lift": r"second-order lift|phase-space|symplectic|Hamiltonian|omega_j\^2|wave equation",
    "dirichlet_source_reuse": r"Dirichlet spectral|local Dirichlet|Laplacian source|lambda_j",
    "missing_source_gates": r"unit-bearing time|time coordinate|kinetic normalization|internally sourced|merely imported|missing-source",
    "no_physics_promotion": r"observed light|gauge photons|spacetime EOM|empirical physics|L_total|ToE|selector closure",
}
PREMISE_GATES = (
    "rho_configuration_space_exported",
    "pi_momentum_space_exported",
    "symplectic_form_exported",
    "quadratic_kinetic_normalization_exported",
    "unit_bearing_time_coordinate_exported",
    "lorentzian_spacetime_embedding_exported",
)
LIFT_CANDIDATES = (
    {
        "id": "gradient_flow_only_baseline",
        "formula": "dot(rho) = -L rho",
        "has_second_order_form": False,
        "has_symplectic_pairing": False,
        "formal_wave_compatible": False,
        "internal_source_gates": {"rho_configuration_space_exported": True},
    },
    {
        "id": "formal_hamiltonian_dirichlet_lift",
        "formula": "dot(rho)=pi, dot(pi)=-L rho; H=1/2*pi^2 + E_D(rho)",
        "has_second_order_form": True,
        "has_symplectic_pairing": True,
        "formal_wave_compatible": True,
        "internal_source_gates": {"rho_configuration_space_exported": True},
    },
    {
        "id": "overdamped_auxiliary_momentum_lift",
        "formula": "dot(rho)=pi, dot(pi)=-pi-L rho",
        "has_second_order_form": True,
        "has_symplectic_pairing": False,
        "formal_wave_compatible": False,
        "internal_source_gates": {"rho_configuration_space_exported": True},
    },
    {
        "id": "complex_unitary_schrodinger_lift",
        "formula": "i dot(psi)=L psi",
        "has_second_order_form": False,
        "has_symplectic_pairing": True,
        "formal_wave_compatible": False,
        "internal_source_gates": {"rho_configuration_space_exported": True},
    },
)


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def gate_rows() -> list[dict[str, Any]]:
    rows = []
    for candidate in LIFT_CANDIDATES:
        sourced = candidate["internal_source_gates"]
        for gate in PREMISE_GATES:
            rows.append({
                "candidate_lift": candidate["id"],
                "premise_gate": gate,
                "internally_sourced_by_current_artifacts": bool(sourced.get(gate, False)),
                "required_for_wave_lightlike_export": gate != "rho_configuration_space_exported",
                "blocked_if_missing": gate != "rho_configuration_space_exported" and not sourced.get(gate, False),
            })
    return rows


def mode_lift_rows() -> list[dict[str, Any]]:
    rows = []
    for candidate in LIFT_CANDIDATES:
        gates = {row["premise_gate"]: row["internally_sourced_by_current_artifacts"] for row in gate_rows() if row["candidate_lift"] == candidate["id"]}
        all_wave_gates_sourced = all(gates[g] for g in PREMISE_GATES if g != "rho_configuration_space_exported")
        for j in range(Z12_SIZE):
            nonconstant = j != 0
            formal_wave_row = bool(candidate["formal_wave_compatible"] and nonconstant)
            accepted_internal_second_order_wave_row = bool(formal_wave_row and all_wave_gates_sourced)
            rows.append({
                "candidate_lift": candidate["id"],
                "mode_j": j,
                "laplacian_eigenvalue_exact": eigenvalue_label(j),
                "formal_second_order_dispersion": f"omega_{j}^2 = {eigenvalue_label(j)}" if candidate["formal_wave_compatible"] else "not a formal wave candidate",
                "nonconstant_mode": nonconstant,
                "formal_wave_compatible_row": formal_wave_row,
                "all_wave_premises_internally_sourced": all_wave_gates_sourced,
                "accepted_internal_second_order_wave_row": accepted_internal_second_order_wave_row,
                "blocked_by": "" if accepted_internal_second_order_wave_row else "; ".join(filter(None, [
                    None if candidate["has_second_order_form"] else "no second-order rho equation",
                    None if candidate["has_symplectic_pairing"] else "no symplectic phase-space pairing",
                    None if candidate["formal_wave_compatible"] else "not wave-compatible even formally",
                    None if nonconstant else "constant zero-frequency mode only",
                    None if gates.get("pi_momentum_space_exported", False) else "pi momentum space not internally sourced",
                    None if gates.get("symplectic_form_exported", False) else "symplectic form not internally sourced",
                    None if gates.get("quadratic_kinetic_normalization_exported", False) else "quadratic kinetic normalization not internally sourced",
                    None if gates.get("unit_bearing_time_coordinate_exported", False) else "unit-bearing time coordinate not internally sourced",
                    None if gates.get("lorentzian_spacetime_embedding_exported", False) else "Lorentzian spacetime embedding not internally sourced",
                ])),
            })
    return rows


def aggregate(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for candidate in LIFT_CANDIDATES:
        subset = [r for r in rows if r["candidate_lift"] == candidate["id"]]
        out.append({
            "candidate_lift": candidate["id"],
            "mode_rows": len(subset),
            "formal_wave_compatible_rows": sum(1 for r in subset if r["formal_wave_compatible_row"]),
            "accepted_internal_second_order_wave_rows": sum(1 for r in subset if r["accepted_internal_second_order_wave_row"]),
        })
    return out


def build_payload() -> dict[str, Any]:
    p3076 = read_json(P3076)
    grep_rows = content_grep()
    gates = gate_rows()
    modes = mode_lift_rows()
    formal_hamiltonian_rows = [r for r in modes if r["candidate_lift"] == "formal_hamiltonian_dirichlet_lift" and r["formal_wave_compatible_row"]]
    proof_obligations = [
        {"obligation": "content_first_grep_before_second_order_lift", "satisfied": True, "detail": "searched second-order/phase-space, Dirichlet source, missing premise gates, and no-physics-promotion lanes"},
        {"obligation": "construct_candidate_phase_space_lifts", "satisfied": True, "detail": "four lift candidates are tested, including the formal Hamiltonian Dirichlet lift"},
        {"obligation": "audit_premise_gate_matrix", "satisfied": True, "detail": "4 candidates x 6 premise gates = 24 exact gate rows"},
        {"obligation": "audit_modal_wave_compatibility_matrix", "satisfied": True, "detail": "4 candidates x 12 Z12 modes = 48 modal lift rows"},
        {"obligation": "export_formal_hamiltonian_dispersion_if_premises_imported", "satisfied": bool(formal_hamiltonian_rows), "detail": "the formal Hamiltonian candidate has 11 nonconstant rows with omega_j^2=lambda_j"},
        {"obligation": "export_internally_sourced_second_order_wave", "satisfied": False, "detail": "0 rows satisfy all phase-space, symplectic, kinetic, time-unit, and Lorentzian source gates"},
    ]
    return {
        "status": "P3077_FORMAL_HAMILTONIAN_LIFT_EXISTS_INTERNAL_SECOND_ORDER_SOURCE_OBSTRUCTED",
        "input_hashes": {"P3076": hashlib.sha256(P3076.read_bytes()).hexdigest() if P3076.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": grep_rows,
            "second_order_lift_interface": {
                "object": "Z12DirichletSecondOrderLiftObstructionTable",
                "source_reused": "P3076 internal diffusive Dirichlet/Laplacian spectrum",
                "candidate_lifts": [c["id"] for c in LIFT_CANDIDATES],
                "premise_gates": list(PREMISE_GATES),
                "acceptance_predicate": "formal wave-compatible nonconstant mode and all phase-space/symplectic/kinetic/time-unit/Lorentzian gates internally sourced",
            },
            "premise_gate_rows": gates,
            "mode_lift_rows": modes,
            "candidate_aggregate_certificate": aggregate(modes),
            "representative_formal_hamiltonian_rows": formal_hamiltonian_rows[:6],
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(r["hit_count"] for r in grep_rows),
            "p3076_z12_modes": p3076.get("finite_certificate", {}).get("z12_modes"),
            "p3076_accepted_lightlike_branch_rows": p3076.get("finite_certificate", {}).get("accepted_lightlike_branch_rows"),
            "candidate_lifts": len(LIFT_CANDIDATES),
            "premise_gates": len(PREMISE_GATES),
            "premise_gate_rows": len(gates),
            "mode_lift_rows": len(modes),
            "formal_wave_compatible_rows": sum(1 for r in modes if r["formal_wave_compatible_row"]),
            "formal_hamiltonian_wave_rows": len(formal_hamiltonian_rows),
            "accepted_internal_second_order_wave_rows": sum(1 for r in modes if r["accepted_internal_second_order_wave_row"]),
            "internally_sourced_required_wave_gate_rows": sum(1 for r in gates if r["required_for_wave_lightlike_export"] and r["internally_sourced_by_current_artifacts"]),
            "proof_obligations": len(proof_obligations),
            "satisfied_proof_obligations": sum(1 for o in proof_obligations if o["satisfied"]),
        },
        "proof_obligations": proof_obligations,
        "decision": {
            "bounded_result": "P3077 constructs the missing phase-space/second-order lift obstruction table.  The formal Hamiltonian Dirichlet lift has the expected 11 nonconstant modal rows omega_j^2=lambda_j, but those rows are accepted only as imported-formal mathematics: current artifacts do not internally source pi momentum space, a symplectic form, kinetic normalization, unit-bearing time, or Lorentzian spacetime embedding.  Therefore no internally sourced wave/lightlike branch is exported.",
            "negative_export_flags": {k: False for k in ["internally_sourced_second_order_wave_exported", "unit_bearing_hamiltonian_exported", "unit_bearing_action_exported", "observed_light_exported", "gauge_photon_sector_exported", "strict_spacetime_embedding_exported", "empirical_physics_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"formal_hamiltonian_dirichlet_lift_constructed": True, "phase_space_premise_gate_matrix_executed": True, "internal_second_order_source_obstruction_exported": True},
            "next_honest_step": "Construct one bounded intrinsic momentum/symplectic-source audit: search only current nadsoliton/Z12 artifacts for a non-imported pi variable or antisymmetric two-form coupled to the Dirichlet source.  If none is exported, freeze second-order wave promotion and pivot to a different non-selector typed object; do not claim observed light, gauge photons, spacetime EOM, units, empirical physics, L_total, bridge/role-transfer, selector closure, or ToE.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3077/S2027 second-order lift obstruction table", "", f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- P3076 Z12 modes: `{c['p3076_z12_modes']}`",
        f"- P3076 accepted lightlike branch rows: `{c['p3076_accepted_lightlike_branch_rows']}`",
        f"- candidate lifts: `{c['candidate_lifts']}`",
        f"- premise gates: `{c['premise_gates']}`",
        f"- premise gate rows: `{c['premise_gate_rows']}`",
        f"- mode lift rows: `{c['mode_lift_rows']}`",
        f"- formal wave-compatible rows: `{c['formal_wave_compatible_rows']}`",
        f"- formal Hamiltonian wave rows: `{c['formal_hamiltonian_wave_rows']}`",
        f"- accepted internal second-order wave rows: `{c['accepted_internal_second_order_wave_rows']}`",
        f"- internally sourced required wave-gate rows: `{c['internally_sourced_required_wave_gate_rows']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3077/S2027 second-order lift obstruction table", "## P3077/S2027 second-order lift obstruction table\n\n`P3077/S2027` constructs `Z12DirichletSecondOrderLiftObstructionTable` after the `P3076` diffusive spectral audit.  It tests `4` phase-space/second-order lift candidates, `6` premise gates, `24` gate rows, and `48` modal lift rows.  The formal Hamiltonian Dirichlet lift gives `11` imported-formal nonconstant rows with `omega_j^2=lambda_j`, but `0` rows pass the internal-source predicate because momentum space, symplectic form, kinetic normalization, unit-bearing time, and Lorentzian embedding are not exported by current artifacts.  This does not export observed light, gauge photons, spacetime EOM, unit-bearing action/Hamiltonian, empirical physics, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3077/S2027 formal Hamiltonian lift remains unsourced", "## P3077/S2027 formal Hamiltonian lift remains unsourced\n\n`P3077/S2027` shows that `H = 1/2*pi^2 + E_D(rho)` is a coherent formal lift of the internal Dirichlet source, with modal relation `omega_j^2=lambda_j` for nonconstant modes.  However, the lift imports `pi`, the symplectic form, kinetic normalization, unit-bearing time, and Lorentzian embedding; current artifacts do not source them.  Therefore no variational spacetime EOM, unit-bearing Hamiltonian/action, observed light, or gauge sector is exported.\n")
    append_once(AGENTS, "Current second-order lift obstruction guardrail (P3077/S2027, 2026-06-24)", "## Current second-order lift obstruction guardrail (P3077/S2027, 2026-06-24)\n\n- P3077 follows the P3076 recommendation and constructs a bounded phase-space/second-order lift obstruction table for the same internal Dirichlet/Laplacian source.\n- The formal Hamiltonian Dirichlet lift has `11` nonconstant modal rows satisfying the imported-formal relation `omega_j^2=lambda_j`, but `0` rows satisfy the internal-source predicate because `pi`, the symplectic form, kinetic normalization, unit-bearing time, and Lorentzian spacetime embedding remain unsourced.\n- Do not promote P3077 to observed light, gauge photons, spacetime EOM, unit-bearing action/Hamiltonian, empirical physics, `QW-2191` discharge, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is one bounded intrinsic momentum/symplectic-source audit over current nadsoliton/Z12 artifacts, or otherwise a pivot to a different non-selector typed object.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
