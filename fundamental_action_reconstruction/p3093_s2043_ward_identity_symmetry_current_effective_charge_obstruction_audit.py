#!/usr/bin/env python3
"""P3093/S2043: Ward-identity/symmetry-current effective-charge audit.

P3092 left finite RG/scale-flow-like witnesses but no sourced physical running
coupling law.  P3093 attacks exactly one adjacent standard-physics interface
atom: whether the Z12 Dirichlet/Laplacian branch internally sources a Ward
identity, conserved symmetry current, or effective charge without importing
continuum gauge theory, spacetime EOM, observed photons/light, apparatus units,
selector closure, L_total, bridge/role-transfer, or ToE.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3092_s2042_renormalization_scale_flow_effective_coupling_obstruction_audit import OUT as P3092, Z12_SIZE

OUT = GEN / "p3093_s2043_ward_identity_symmetry_current_effective_charge_obstruction_audit.json"
MD = GEN / "p3093_s2043_ward_identity_symmetry_current_effective_charge_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
SOURCE_MASKS = (1, 3, 5, 15, 31, 63, 127, 255, 511, 1023, 2047, 4095)
CHARGE_MODES = tuple(range(Z12_SIZE))

CONTENT_PATTERNS = {
    "ward_symmetry_current_atom": r"Ward|symmetry-current|conserved current|effective charge|gauge-charge",
    "predecessor_rg_atom": r"renormalization|scale-flow|running coupling|beta function|effective-coupling",
    "blocked_promotions": r"continuum gauge theory|observed photons|observed light|spacetime EOM|L_total|ToE|selector closure|bridge/role-transfer",
}

CHARGE_CANDIDATES = (
    {"id": "z12_graph_automorphism_commutant", "description": "translations/reflections commuting with the Z12 circulant Laplacian", "finite_graph_symmetry_exported": True, "laplacian_commutator_identity_exported": True, "source_variation_identity_exported": False, "conserved_current_exported": False, "ward_identity_exported": False, "gauge_charge_normalization_exported": False, "empirical_charge_readout_exported": False, "nonimported_physical_charge_source_exported": False, "blocker": "graph automorphism commutation is kinematic and does not source a current, Ward identity, gauge charge, or readout"},
    {"id": "finite_source_variation_balance", "description": "finite source-profile variation balances under translations", "finite_graph_symmetry_exported": True, "laplacian_commutator_identity_exported": True, "source_variation_identity_exported": True, "conserved_current_exported": False, "ward_identity_exported": False, "gauge_charge_normalization_exported": False, "empirical_charge_readout_exported": False, "nonimported_physical_charge_source_exported": False, "blocker": "source-balance rows are formal profile identities and lack current conservation semantics"},
    {"id": "spectral_mode_charge_labels", "description": "integer Fourier-mode labels used as formal charge labels", "finite_graph_symmetry_exported": True, "laplacian_commutator_identity_exported": True, "source_variation_identity_exported": False, "conserved_current_exported": False, "ward_identity_exported": False, "gauge_charge_normalization_exported": False, "empirical_charge_readout_exported": False, "nonimported_physical_charge_source_exported": False, "blocker": "mode labels are representational charges only and do not fix physical normalization or readout"},
    {"id": "imported_continuum_gauge_ward_template", "description": "external continuum gauge-theory Ward identity template", "finite_graph_symmetry_exported": True, "laplacian_commutator_identity_exported": True, "source_variation_identity_exported": True, "conserved_current_exported": True, "ward_identity_exported": True, "gauge_charge_normalization_exported": False, "empirical_charge_readout_exported": False, "nonimported_physical_charge_source_exported": False, "blocker": "Ward/current semantics are imported and still lack internal charge units/readout"},
    {"id": "imported_empirical_electric_charge_template", "description": "external empirically calibrated electric-charge/readout template", "finite_graph_symmetry_exported": True, "laplacian_commutator_identity_exported": True, "source_variation_identity_exported": True, "conserved_current_exported": True, "ward_identity_exported": True, "gauge_charge_normalization_exported": True, "empirical_charge_readout_exported": True, "nonimported_physical_charge_source_exported": False, "blocker": "passes only by importing electric-charge normalization and empirical apparatus semantics"},
)
REQUIRED_GATES = ("finite_graph_symmetry_exported", "laplacian_commutator_identity_exported", "source_variation_identity_exported", "conserved_current_exported", "ward_identity_exported", "gauge_charge_normalization_exported", "empirical_charge_readout_exported", "nonimported_physical_charge_source_exported")


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def laplacian_entry(i: int, j: int) -> int:
    if i == j:
        return 2
    if (i - j) % Z12_SIZE in (1, Z12_SIZE - 1):
        return -1
    return 0


def apply_affine(unit: int, shift: int, site: int) -> int:
    return (unit * site + shift) % Z12_SIZE


def symmetry_commutator_rows() -> list[dict[str, Any]]:
    rows = []
    for unit in (1, -1):
        for shift in range(Z12_SIZE):
            max_abs = 0
            trace_defect = 0
            for i in range(Z12_SIZE):
                for j in range(Z12_SIZE):
                    diff = laplacian_entry(apply_affine(unit, shift, i), apply_affine(unit, shift, j)) - laplacian_entry(i, j)
                    max_abs = max(max_abs, abs(diff))
                    if i == j:
                        trace_defect += diff
            rows.append({"affine_unit": unit, "translation_shift": shift, "commutator_max_abs_entry": max_abs, "trace_defect": trace_defect, "laplacian_commutator_identity_witness": max_abs == 0, "current_or_charge_semantics_attached": False})
    return rows


def profile(mask: int) -> list[int]:
    return [int((mask >> i) & 1) for i in range(Z12_SIZE)]


def source_variation_rows() -> list[dict[str, Any]]:
    rows = []
    for mask in SOURCE_MASKS:
        p = profile(mask)
        forward_diff = [p[(i + 1) % Z12_SIZE] - p[i] for i in range(Z12_SIZE)]
        rows.append({"source_mask": mask, "support_size": sum(p), "cyclic_variation_sum": sum(forward_diff), "total_variation_l1": sum(abs(x) for x in forward_diff), "finite_source_variation_balance_witness": sum(forward_diff) == 0, "conserved_current_attached": False, "ward_identity_attached": False})
    return rows


def spectral_charge_rows() -> list[dict[str, Any]]:
    rows = []
    for k in CHARGE_MODES:
        paired = (-k) % Z12_SIZE
        lam = 2.0 - 2.0 * math.cos(2.0 * math.pi * k / Z12_SIZE)
        rows.append({"mode": k, "opposite_mode": paired, "laplacian_eigenvalue": round(lam, 12), "integer_mode_charge_label": k, "opposite_label": paired, "degenerate_with_opposite": abs(lam - (2.0 - 2.0 * math.cos(2.0 * math.pi * paired / Z12_SIZE))) < 1e-12, "physical_charge_normalization_attached": False, "empirical_charge_readout_attached": False})
    return rows


def discrete_current_rows() -> list[dict[str, Any]]:
    p = profile(0b001111001111)
    edge_current = [p[(i + 1) % Z12_SIZE] - p[i] for i in range(Z12_SIZE)]
    rows = []
    for i in range(Z12_SIZE):
        divergence = edge_current[i] - edge_current[(i - 1) % Z12_SIZE]
        rows.append({"site": i, "formal_edge_current_out": edge_current[i], "formal_edge_current_in": edge_current[(i - 1) % Z12_SIZE], "formal_divergence": divergence, "local_conservation_witness": divergence == 0, "current_is_formal_profile_gradient": True, "physical_current_law_attached": False})
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in CHARGE_CANDIDATES for gate in REQUIRED_GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for c in CHARGE_CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        out.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_nonimported_ward_current_effective_charge_source": all(r["gate_passed"] for r in subset) and bool(c["nonimported_physical_charge_source_exported"])})
    return out


def build_payload() -> dict[str, Any]:
    p3092 = read_json(P3092)
    greps = content_grep(); syms = symmetry_commutator_rows(); vars_ = source_variation_rows(); charges = spectral_charge_rows(); currents = discrete_current_rows(); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_nonimported_ward_current_effective_charge_source"]]
    obligations = [
        {"obligation": "read_p3092_next_atom", "satisfied": True, "detail": "P3092 selected Ward-identity/symmetry-current effective-charge audit as the next interface atom"},
        {"obligation": "construct_finite_graph_symmetry_commutators", "satisfied": len(syms) == 24 and all(r["laplacian_commutator_identity_witness"] for r in syms), "detail": "translations/reflections commute with the finite Z12 Laplacian"},
        {"obligation": "construct_source_variation_balance_rows", "satisfied": len(vars_) == len(SOURCE_MASKS) and all(r["finite_source_variation_balance_witness"] for r in vars_), "detail": "cyclic finite-difference source variations have zero total variation sum"},
        {"obligation": "construct_spectral_charge_and_current_rows", "satisfied": len(charges) == Z12_SIZE and len(currents) == Z12_SIZE, "detail": "mode-label charge rows and formal edge-current rows are explicit"},
        {"obligation": "export_nonimported_physical_ward_current_charge_source", "satisfied": False, "detail": "0 candidates pass as internal non-imported Ward/current/effective-charge sources"},
    ]
    return {
        "status": "P3093_WARD_IDENTITY_SYMMETRY_CURRENT_EFFECTIVE_CHARGE_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3092": hashlib.sha256(P3092.read_bytes()).hexdigest() if P3092.exists() else None},
        "constructed_theoretical_objects": {"content_first_repo_grep": greps, "ward_current_effective_charge_audit_object": {"object": "Z12DirichletWardIdentitySymmetryCurrentEffectiveChargeObstructionAudit", "source_reused": "P3092 recommendation: bounded Ward-identity/symmetry-current effective-charge obstruction audit", "required_gates": list(REQUIRED_GATES), "candidate_charge_sources": [c["id"] for c in CHARGE_CANDIDATES], "acceptance_predicate": "finite graph symmetry plus Laplacian commutator identity, source variation identity, conserved current, Ward identity, gauge-charge normalization, empirical charge/readout, and non-imported physical charge source"}, "symmetry_commutator_rows": syms, "source_variation_balance_rows": vars_, "spectral_charge_label_rows": charges, "formal_discrete_current_rows": currents, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs},
        "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps), "p3092_accepted_nonimported_renormalization_scale_flow_sources": p3092.get("finite_certificate", {}).get("accepted_nonimported_renormalization_scale_flow_sources"), "symmetry_commutator_rows": len(syms), "symmetry_commutator_failures": sum(not r["laplacian_commutator_identity_witness"] for r in syms), "symmetry_rows_with_current_or_charge_semantics": sum(r["current_or_charge_semantics_attached"] for r in syms), "source_variation_rows": len(vars_), "source_variation_balance_failures": sum(not r["finite_source_variation_balance_witness"] for r in vars_), "source_variation_rows_with_ward_identity": sum(r["ward_identity_attached"] for r in vars_), "spectral_charge_label_rows": len(charges), "spectral_charge_rows_with_physical_normalization": sum(r["physical_charge_normalization_attached"] for r in charges), "spectral_charge_rows_with_empirical_readout": sum(r["empirical_charge_readout_attached"] for r in charges), "formal_discrete_current_rows": len(currents), "formal_current_rows_with_physical_current_law": sum(r["physical_current_law_attached"] for r in currents), "charge_candidates": len(CHARGE_CANDIDATES), "required_gates": len(REQUIRED_GATES), "candidate_gate_rows": len(gates), "accepted_nonimported_ward_current_effective_charge_sources": len(accepted), "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations)},
        "proof_obligations": obligations,
        "decision": {"bounded_result": "P3093 constructs the requested Ward-identity/symmetry-current effective-charge obstruction audit.  The Z12 Laplacian supports exact finite translation/reflection commutator identities, cyclic source-variation balance rows, formal Fourier-mode charge labels, and formal discrete edge-current rows.  These are real symmetry/current-like witnesses, but no internal artifact exports a physical conserved current, Ward identity, gauge-charge normalization, empirical charge/readout semantics, spacetime EOM, observed photons/light, or a non-imported physical effective-charge source.  Imported continuum gauge-theory and empirical electric-charge templates pass only as imported templates.  Therefore no non-imported Ward/current/effective-charge source is exported.", "negative_export_flags": {key: False for key in ["conserved_current_exported", "ward_identity_exported", "gauge_charge_normalization_exported", "empirical_charge_readout_exported", "nonimported_physical_charge_source_exported", "physical_action_functional_exported", "empirical_observable_exported", "observed_radiation_exported", "observed_photons_exported", "observed_light_exported", "spacetime_eom_exported", "physical_hamiltonian_exported", "green_response_closure_exported", "effective_action_closure_exported", "rg_scale_flow_closure_exported", "born_rule_readout_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"finite_graph_symmetry_commutators_computed": True, "source_variation_balance_rows_computed": True, "spectral_charge_label_rows_computed": True, "formal_discrete_current_rows_computed": True}, "next_honest_step": "Pivot to exactly one adjacent standard-physics interface atom outside selector replay: construct a bounded stress-energy/metric-response obstruction audit for the Z12 Dirichlet/Laplacian branch, testing whether finite Laplacian metric variations, graph-energy quadratic forms, symmetry currents, and spectral pressure-like derivatives supply a non-imported stress tensor, metric coupling, conservation law, and empirical gravitational/field-response interface without importing spacetime geometry, continuum EOM, observed radiation/light, apparatus units, L_total, bridge/role-transfer, or ToE."},
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3093/S2043 Ward-identity/symmetry-current effective-charge obstruction audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- P3092 accepted non-imported renormalization/scale-flow sources: `{c['p3092_accepted_nonimported_renormalization_scale_flow_sources']}`", f"- symmetry commutator rows: `{c['symmetry_commutator_rows']}`", f"- symmetry commutator failures: `{c['symmetry_commutator_failures']}`", f"- symmetry rows with current or charge semantics: `{c['symmetry_rows_with_current_or_charge_semantics']}`", f"- source variation rows: `{c['source_variation_rows']}`", f"- source variation balance failures: `{c['source_variation_balance_failures']}`", f"- source variation rows with Ward identity: `{c['source_variation_rows_with_ward_identity']}`", f"- spectral charge label rows: `{c['spectral_charge_label_rows']}`", f"- spectral charge rows with physical normalization: `{c['spectral_charge_rows_with_physical_normalization']}`", f"- spectral charge rows with empirical readout: `{c['spectral_charge_rows_with_empirical_readout']}`", f"- formal discrete current rows: `{c['formal_discrete_current_rows']}`", f"- formal current rows with physical current law: `{c['formal_current_rows_with_physical_current_law']}`", f"- charge candidates: `{c['charge_candidates']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted non-imported Ward/current/effective-charge sources: `{c['accepted_nonimported_ward_current_effective_charge_sources']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3093/S2043 Ward-identity/symmetry-current effective-charge obstruction audit", "## P3093/S2043 Ward-identity/symmetry-current effective-charge obstruction audit\n\n`P3093/S2043` attacks exactly one post-P3092 interface atom: a non-imported Ward-identity/symmetry-current effective-charge source for the Z12 Dirichlet/Laplacian branch.  It constructs `24` exact translation/reflection commutator rows, `12` source-variation balance rows, `12` spectral charge-label rows, `12` formal discrete-current rows, and a `5 x 8 = 40` candidate gate matrix.  The finite symmetry/current-like algebra remains formal; no physical conserved current, Ward identity, gauge-charge normalization, empirical charge/readout interface, spacetime EOM, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3093/S2043 Ward/current/effective-charge source remains unsourced", "## P3093/S2043 Ward/current/effective-charge source remains unsourced\n\n`P3093/S2043` confirms that the Z12 Laplacian supports exact graph-symmetry commutators, source-variation balances, Fourier-mode charge labels, and formal edge-current rows.  A Lagrangian/EOM, gauge-charge, or empirical current reading still needs strict sources for a physical current law, Ward identity, charge normalization, empirical readout, and spacetime/EOM interpretation; imported continuum gauge templates do not supply those strict sources.\n")
    append_once(AGENTS, "Current Ward-identity/symmetry-current effective-charge obstruction guardrail (P3093/S2043, 2026-06-25)", "## Current Ward-identity/symmetry-current effective-charge obstruction guardrail (P3093/S2043, 2026-06-25)\n\n- P3093 follows the P3092 recommendation and audits one standard-physics interface atom: a non-imported Ward-identity/symmetry-current effective-charge source for the Z12 Dirichlet/Laplacian branch.\n- The finite audit computes `24` exact translation/reflection commutator rows, `12` source-variation balance rows, `12` spectral charge-label rows, `12` formal discrete-current rows, and `40` candidate gate rows; `0` candidates export an internal non-imported Ward/current/effective-charge law.\n- Do not promote graph automorphism commutators, finite source-variation balances, Fourier-mode charge labels, formal edge-current rows, imported continuum gauge Ward templates, or imported empirical electric-charge templates to physical conserved current, gauge-charge normalization, empirical observable, observed photons/light, spacetime EOM, physical Hamiltonian, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one stress-energy/metric-response obstruction audit for the Z12 Dirichlet/Laplacian branch, unless a genuinely new Ward/current/effective-charge theorem is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
