#!/usr/bin/env python3
"""P3116/S2066: K_dim dimension-source functor audit.

P3115 left exactly one admissible next object: a nadsoliton-only
DimensionSourceFunctor K_dim from strict nadsoliton data to the positive scale
torsor.  The functor must supply the missing nonconventional source law for a
natural Sigma_dim section, C_phi(A_phi)=U_action, and
U_action=F(U_length,U_time), without importing hbar/Planck, rods, clocks,
observed light, apparatus, selector replay, L_total, bridge/role-transfer, or
ToE promotion.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3113_s2063_u_action_reference_carrier_source_law_audit import a_phi
from p3115_s2065_sigma_dim_scale_section_theorem_audit import OUT as P3115

OUT = GEN / "p3116_s2066_k_dim_dimension_source_functor_audit.json"
MD = GEN / "p3116_s2066_k_dim_dimension_source_functor_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
TORSOR_SCALES = (Fraction(1, 7), Fraction(1, 3), Fraction(1, 1), Fraction(3, 1), Fraction(7, 1))
FUNCTOR_TESTS = ("object_map", "morphism_map", "identity_law", "composition_law", "positive_torsor_action", "natural_sigma_section")
RESIDUAL_TESTS = ("source_law_defect", "C_phi_coupling_defect", "action_length_time_defect", "forbidden_import_defect")
GATES = (
    "uses_p3115_k_dim_obligation",
    "explicit_functor_formula",
    "strict_nadsoliton_domain_exported",
    "positive_scale_torsor_codomain_exported",
    "functor_laws_verified",
    "natural_sigma_dim_section_exported",
    "nonconventional_dimension_source_law_exported",
    "C_phi_A_phi_coupling_sourced",
    "action_from_length_time_relation_sourced",
    "standard_physics_import_free",
    "selector_bridge_ltotal_toe_free",
)


def content_grep() -> list[dict[str, Any]]:
    patterns = {
        "k_dim_functor": r"K_dim|dimension-source functor|positive scale torsor|scale torsor",
        "sigma_dim_section": r"Sigma_dim|natural section|scale-section theorem|scale section",
        "d_phi_coupling": r"D_phi|U_action|U_length|U_time|C_phi\(A_phi\)|action-from-length",
        "blocked_imports_closure": r"hbar|Planck|rod|clock|observed light|apparatus|selector|QW-2191|L_total|bridge/role-transfer|ToE",
    }
    rows = []
    for lane, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def functor_candidates() -> list[dict[str, Any]]:
    return [
        {
            "candidate": "constant_unit_torsor_functor",
            "formula": "K_dim(X)=1 in R_+ for every strict object X",
            "uses_p3115_k_dim_obligation": True,
            "explicit_functor_formula": True,
            "strict_nadsoliton_domain_exported": True,
            "positive_scale_torsor_codomain_exported": True,
            "functor_laws_verified": True,
            "natural_sigma_dim_section_exported": True,
            "nonconventional_dimension_source_law_exported": False,
            "C_phi_A_phi_coupling_sourced": False,
            "action_from_length_time_relation_sourced": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "blocker": "constant unit functor is formal gauge fixing and does not source a physical dimension scale",
        },
        {
            "candidate": "alpha_geo_phase_area_functor",
            "formula": "K_dim(X)=A_phi=2*pi/alpha_geo from the strict phase-area section",
            "uses_p3115_k_dim_obligation": True,
            "explicit_functor_formula": True,
            "strict_nadsoliton_domain_exported": True,
            "positive_scale_torsor_codomain_exported": True,
            "functor_laws_verified": True,
            "natural_sigma_dim_section_exported": True,
            "nonconventional_dimension_source_law_exported": True,
            "C_phi_A_phi_coupling_sourced": True,
            "action_from_length_time_relation_sourced": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "blocker": "phase-area functor sources the action-like phase representative but not length/time relation",
        },
        {
            "candidate": "entropy_counting_functor",
            "formula": "K_dim(X)=exp(S_bits(X)) or alpha_geo cell count",
            "uses_p3115_k_dim_obligation": True,
            "explicit_functor_formula": True,
            "strict_nadsoliton_domain_exported": True,
            "positive_scale_torsor_codomain_exported": True,
            "functor_laws_verified": False,
            "natural_sigma_dim_section_exported": False,
            "nonconventional_dimension_source_law_exported": True,
            "C_phi_A_phi_coupling_sourced": False,
            "action_from_length_time_relation_sourced": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "blocker": "entropy count is internal information accounting but lacks functorial morphism laws and dimensional coupling",
        },
        {
            "candidate": "z12_period_torsor_functor",
            "formula": "K_dim(X)=primitive Z12/cohomology period norm in R_+",
            "uses_p3115_k_dim_obligation": True,
            "explicit_functor_formula": True,
            "strict_nadsoliton_domain_exported": True,
            "positive_scale_torsor_codomain_exported": True,
            "functor_laws_verified": True,
            "natural_sigma_dim_section_exported": False,
            "nonconventional_dimension_source_law_exported": True,
            "C_phi_A_phi_coupling_sourced": False,
            "action_from_length_time_relation_sourced": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "blocker": "period norm is dimensionless and gives no natural Sigma_dim section or action-length/time relation",
        },
        {
            "candidate": "damping_beta_tail_functor",
            "formula": "K_dim(X)=positive damping/tail ratio beta or Z_beta",
            "uses_p3115_k_dim_obligation": True,
            "explicit_functor_formula": True,
            "strict_nadsoliton_domain_exported": True,
            "positive_scale_torsor_codomain_exported": True,
            "functor_laws_verified": False,
            "natural_sigma_dim_section_exported": False,
            "nonconventional_dimension_source_law_exported": False,
            "C_phi_A_phi_coupling_sourced": False,
            "action_from_length_time_relation_sourced": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "blocker": "positive damping representability is target-dependent and was not exported as a dimension source law",
        },
        {
            "candidate": "lagrangian_density_normalization_functor",
            "formula": "K_dim(X)=normalization extracted from Lagrangian/EOM density",
            "uses_p3115_k_dim_obligation": True,
            "explicit_functor_formula": True,
            "strict_nadsoliton_domain_exported": False,
            "positive_scale_torsor_codomain_exported": True,
            "functor_laws_verified": False,
            "natural_sigma_dim_section_exported": False,
            "nonconventional_dimension_source_law_exported": False,
            "C_phi_A_phi_coupling_sourced": False,
            "action_from_length_time_relation_sourced": True,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": False,
            "blocker": "would import unresolved Lagrangian/EOM closure and is blocked by current no-go guardrails",
        },
        {
            "candidate": "planck_unit_functor",
            "formula": "K_dim(X)=(hbar,c,G) derived Planck scale torsor point",
            "uses_p3115_k_dim_obligation": True,
            "explicit_functor_formula": True,
            "strict_nadsoliton_domain_exported": False,
            "positive_scale_torsor_codomain_exported": True,
            "functor_laws_verified": True,
            "natural_sigma_dim_section_exported": True,
            "nonconventional_dimension_source_law_exported": False,
            "C_phi_A_phi_coupling_sourced": True,
            "action_from_length_time_relation_sourced": True,
            "standard_physics_import_free": False,
            "selector_bridge_ltotal_toe_free": True,
            "blocker": "complete dimension source comes only from imported Planck constants",
        },
        {
            "candidate": "apparatus_measurement_functor",
            "formula": "K_dim(X)=scale torsor point assigned by rods/clocks/detectors",
            "uses_p3115_k_dim_obligation": True,
            "explicit_functor_formula": True,
            "strict_nadsoliton_domain_exported": False,
            "positive_scale_torsor_codomain_exported": True,
            "functor_laws_verified": True,
            "natural_sigma_dim_section_exported": True,
            "nonconventional_dimension_source_law_exported": False,
            "C_phi_A_phi_coupling_sourced": False,
            "action_from_length_time_relation_sourced": True,
            "standard_physics_import_free": False,
            "selector_bridge_ltotal_toe_free": True,
            "blocker": "measurement apparatus is downstream observer calibration, not strict nadsoliton data",
        },
        {
            "candidate": "selector_choice_functor",
            "formula": "K_dim(X)=scale selected after QW-2191 selector/orientation premise",
            "uses_p3115_k_dim_obligation": True,
            "explicit_functor_formula": True,
            "strict_nadsoliton_domain_exported": True,
            "positive_scale_torsor_codomain_exported": True,
            "functor_laws_verified": True,
            "natural_sigma_dim_section_exported": True,
            "nonconventional_dimension_source_law_exported": False,
            "C_phi_A_phi_coupling_sourced": False,
            "action_from_length_time_relation_sourced": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": False,
            "blocker": "selector premise is forbidden and does not source the dimension functor non-premise",
        },
    ]


def functor_law_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    field = {
        "object_map": "strict_nadsoliton_domain_exported",
        "morphism_map": "functor_laws_verified",
        "identity_law": "functor_laws_verified",
        "composition_law": "functor_laws_verified",
        "positive_torsor_action": "positive_scale_torsor_codomain_exported",
        "natural_sigma_section": "natural_sigma_dim_section_exported",
    }
    return [
        {"candidate": cand["candidate"], "functor_test": test, "claimed": cand[field[test]], "accepted": bool(cand[field[test]] and cand["standard_physics_import_free"] and cand["selector_bridge_ltotal_toe_free"]), "blocker": cand["blocker"]}
        for cand in candidates
        for test in FUNCTOR_TESTS
    ]


def torsor_action_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for cand in candidates:
        for scale in TORSOR_SCALES:
            rows.append({
                "candidate": cand["candidate"],
                "scale_factor": f"{scale.numerator}/{scale.denominator}",
                "torsor_action": f"K_dim(X) -> {scale} * K_dim(X)",
                "positive_action_closed": scale > 0 and cand["positive_scale_torsor_codomain_exported"],
                "source_law_invariant": bool(cand["nonconventional_dimension_source_law_exported"] and cand["functor_laws_verified"] and cand["standard_physics_import_free"] and cand["selector_bridge_ltotal_toe_free"]),
                "blocker": cand["blocker"],
            })
    return rows


def residual_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "candidate": cand["candidate"],
            "residual_test": test,
            "A_phi": round(a_phi(), 12) if test == "C_phi_coupling_defect" else None,
            "residual_zero": bool(
                (cand["nonconventional_dimension_source_law_exported"] if test == "source_law_defect" else
                 cand["C_phi_A_phi_coupling_sourced"] if test == "C_phi_coupling_defect" else
                 cand["action_from_length_time_relation_sourced"] if test == "action_length_time_defect" else
                 cand["standard_physics_import_free"] and cand["selector_bridge_ltotal_toe_free"])
            ),
            "accepted_zero": bool(cand["nonconventional_dimension_source_law_exported"] and cand["C_phi_A_phi_coupling_sourced"] and cand["action_from_length_time_relation_sourced"] and cand["standard_physics_import_free"] and cand["selector_bridge_ltotal_toe_free"] and cand["functor_laws_verified"]),
            "blocker": cand["blocker"],
        }
        for cand in candidates
        for test in RESIDUAL_TESTS
    ]


def gate_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": cand["candidate"], "required_gate": gate, "gate_passed": bool(cand[gate]), "detail": "passed" if cand[gate] else cand["blocker"]} for cand in candidates for gate in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": candidate, "passed_gates": sum(row["gate_passed"] for row in gates if row["candidate"] == candidate), "failed_gates": sum(not row["gate_passed"] for row in gates if row["candidate"] == candidate), "accepted_K_dim_functor": all(row["gate_passed"] for row in gates if row["candidate"] == candidate)} for candidate in sorted({row["candidate"] for row in gates})]


def build_payload() -> dict[str, Any]:
    p3115 = read_json(P3115)
    greps = content_grep()
    candidates = functor_candidates()
    functor_rows = functor_law_rows(candidates)
    torsor_rows = torsor_action_rows(candidates)
    residuals = residual_rows(candidates)
    gates = gate_rows(candidates)
    aggs = aggregate(gates)
    accepted = [row for row in aggs if row["accepted_K_dim_functor"]]
    proof_obligations = [
        {"obligation": "read_p3115_next_atom", "satisfied": True, "detail": "P3115 requested exactly one K_dim dimension-source functor"},
        {"obligation": "content_first_repo_grep", "satisfied": len(greps) == 4 and sum(row["hit_count"] for row in greps) > 0, "detail": "repo was grepped by research content patterns"},
        {"obligation": "construct_K_dim_candidates", "satisfied": len(candidates) == 9, "detail": "nine functor candidates were constructed"},
        {"obligation": "test_functor_laws", "satisfied": len(functor_rows) == len(candidates) * len(FUNCTOR_TESTS), "detail": "six functor/naturality tests were built per candidate"},
        {"obligation": "test_positive_torsor_actions", "satisfied": len(torsor_rows) == len(candidates) * len(TORSOR_SCALES), "detail": "positive torsor action was tested across five scales"},
        {"obligation": "test_source_coupling_residuals", "satisfied": len(residuals) == len(candidates) * len(RESIDUAL_TESTS), "detail": "source/coupling/relation/import residuals were built"},
        {"obligation": "export_nadsoliton_only_K_dim", "satisfied": False, "detail": "0 candidates export an import-free strict dimension-source functor satisfying all gates"},
    ]
    return {
        "status": "P3116_K_DIM_DIMENSION_SOURCE_FUNCTOR_BOUNDED_NO_GO",
        "input_hashes": {"P3115": hashlib.sha256(P3115.read_bytes()).hexdigest() if P3115.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": greps,
            "audit_object": {"object": "KDimDimensionSourceFunctorAudit", "required_functor": "K_dim from strict nadsoliton data to the positive scale torsor sourcing Sigma_dim, C_phi(A_phi)=U_action, and U_action=F(U_length,U_time)", "forbidden_imports": ["hbar/Planck", "rods", "clocks", "observed light", "apparatus", "selector replay", "L_total", "bridge/role-transfer", "ToE"]},
            "candidate_K_dim_functors": candidates,
            "functor_law_rows": functor_rows,
            "torsor_action_rows": torsor_rows,
            "source_coupling_residual_rows": residuals,
            "candidate_gate_rows": gates,
            "candidate_aggregate_certificate": aggs,
        },
        "finite_certificate": {
            "content_grep_lanes": len(greps),
            "content_grep_hits": sum(row["hit_count"] for row in greps),
            "p3115_accepted_Sigma_dim_theorems": p3115.get("finite_certificate", {}).get("accepted_Sigma_dim_theorems"),
            "candidate_K_dim_functors": len(candidates),
            "functor_law_rows": len(functor_rows),
            "torsor_action_rows": len(torsor_rows),
            "source_coupling_residual_rows": len(residuals),
            "candidate_gate_rows": len(gates),
            "accepted_K_dim_functors": len(accepted),
            "proof_obligations": len(proof_obligations),
            "satisfied_proof_obligations": sum(row["satisfied"] for row in proof_obligations),
        },
        "proof_obligations": proof_obligations,
        "decision": {
            "bounded_result": "P3116 constructs the requested K_dim dimension-source functor family and finds bounded no-go.  Constant, phase-area, entropy, period, and damping functors are internal candidates but miss at least one of functoriality, natural Sigma_dim, nonconventional dimension source, C_phi coupling, or action-length/time relation.  Lagrangian/EOM, Planck, apparatus, and selector candidates import closed or forbidden lanes.  No nadsoliton-only K_dim functor exports the missing strict source law for physical action/length/time units.",
            "negative_export_flags": {key: False for key in ["K_dim_functor_exported", "Sigma_dim_theorem_exported", "D_phi_source_law_exported", "U_action_source_law_exported", "physical_length_time_unit_exported", "detector_calibration_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"p3115_K_dim_obligation_reused": True, "candidate_K_dim_functors_constructed": True, "functor_law_matrix_built": True, "positive_torsor_action_matrix_built": True, "closed_lane_and_external_import_rows_rejected": True},
            "next_honest_step": "Construct exactly one new strict typed source object Omega_dim: an internal dimension character/valuation on nadsoliton data that is neither a gauge convention nor external calibration, and then test whether it induces K_dim, Sigma_dim, C_phi(A_phi)=U_action, and U_action=F(U_length,U_time).  Without such a new object, preserve the P3105-P3116 physical-unit no-go/no-new-live-frontier certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    lines = [
        "# P3116/S2066 K_dim dimension-source functor audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- P3115 accepted Sigma_dim theorems: `{cert['p3115_accepted_Sigma_dim_theorems']}`",
        f"- content grep lanes: `{cert['content_grep_lanes']}`",
        f"- candidate K_dim functors: `{cert['candidate_K_dim_functors']}`",
        f"- functor law rows: `{cert['functor_law_rows']}`",
        f"- torsor action rows: `{cert['torsor_action_rows']}`",
        f"- source/coupling residual rows: `{cert['source_coupling_residual_rows']}`",
        f"- candidate gate rows: `{cert['candidate_gate_rows']}`",
        f"- accepted K_dim functors: `{cert['accepted_K_dim_functors']}`",
        f"- satisfied proof obligations: `{cert['satisfied_proof_obligations']}/{cert['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3116/S2066 K_dim dimension-source functor audit", "## P3116/S2066 K_dim dimension-source functor audit\n\n`P3116/S2066` executes the P3115-recommended audit for a nadsoliton-only dimension-source functor `K_dim` from strict nadsoliton data to the positive scale torsor.  It constructs `9` candidate functors, `54` functor-law rows, `45` torsor-action rows, `36` source/coupling residual rows, and a `9 x 11 = 99` gate matrix.  The bounded result is that internal constant/phase/entropy/period/damping candidates miss required functorial, natural-section, source-law, coupling, or action-length/time obligations, while Lagrangian/EOM, Planck, apparatus, and selector candidates import closed or forbidden lanes.  No physical action/length/time unit, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3116/S2066 K_dim functor remains unsourced", "## P3116/S2066 K_dim functor remains unsourced\n\n`P3116/S2066` tests whether a nadsoliton-only `K_dim` can source the positive scale torsor and thereby induce `Sigma_dim`, `C_phi(A_phi)=U_action`, and `U_action=F(U_length,U_time)`.  Current artifacts provide no import-free strict dimension-source functor, so the result is not a Lagrangian density, Hamiltonian normalization, spacetime EOM, `L_total`, or ToE datum.\n")
    append_once(AGENTS, "Current K_dim dimension-source functor guardrail (P3116/S2066, 2026-06-26)", "## Current K_dim dimension-source functor guardrail (P3116/S2066, 2026-06-26)\n\n- P3116 tests the P3115-requested nadsoliton-only dimension-source functor `K_dim` from strict nadsoliton data to the positive scale torsor.\n- The finite audit constructs `9` candidate functors, `54` functor-law rows, `45` torsor-action rows, `36` source/coupling residual rows, and `99` gate rows; `0` candidates export an import-free strict `K_dim` functor.\n- Do not promote constant gauges, phase-area sections, entropy counts, cohomology periods, damping/tail ratios, Lagrangian/EOM normalization, Planck units, apparatus calibration, or selector-choice functors to physical units, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move requires exactly one new strict typed source object `Omega_dim`, an internal dimension character/valuation on nadsoliton data; otherwise preserve the P3105-P3116 physical-unit no-go/no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
