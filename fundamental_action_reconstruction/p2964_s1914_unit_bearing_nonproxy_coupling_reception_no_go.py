#!/usr/bin/env python3
"""P2964/S1914: unit-bearing nonproxy coupling reception no-go.

P2963 left an honest continuation: try to construct an actual unit-bearing
nonproxy action-density coupling that can receive the P2938/P2961 aggregate
without requiring artifact-level K/C exchange.  This audit builds a finite menu
of coupling templates and grades each against the missing source obligations:
aggregate reception, no K/C exchange dependence, canonical unit/scale quotient,
nonproxy field variable, local action density, variational derivative, and a
strict coefficient/source law.

The finite reception side is positive for several templates, but every audited
unit-bearing/nonproxy route still misses at least one strict source obligation.
Therefore no sourced damping coefficient enters L_total from this package.
"""
from __future__ import annotations

import hashlib, json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2961_s1911_exchange_balanced_scale_quotient_source_candidate import OUT as P2961
from p2963_s1913_typed_mediator_functor_no_go import OUT as P2963

OUT = GEN / "p2964_s1914_unit_bearing_nonproxy_coupling_reception_no_go.json"
MD = GEN / "p2964_s1914_unit_bearing_nonproxy_coupling_reception_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

AGGREGATE_VECTOR = [1, 2, 2, 2, 2]
ETA = sum(AGGREGATE_VECTOR) / len(AGGREGATE_VECTOR)
DELTA = ETA - 1


def coupling_rows() -> list[dict[str, Any]]:
    candidates = [
        {
            "template": "dimensionless_aggregate_source_token",
            "formula": "lambda * sum(V) * Phi",
            "receives_aggregate": True,
            "requires_KC_artifact_exchange": False,
            "canonical_scale_quotient": False,
            "unit_bearing": False,
            "nonproxy_field_variable": False,
            "local_action_density": False,
            "variational_derivative_exported": False,
            "strict_coefficient_source_law": False,
        },
        {
            "template": "beta_scale_orbit_normalized_density",
            "formula": "beta(ell) * |Phi|^2 with ell chosen by normalization",
            "receives_aggregate": True,
            "requires_KC_artifact_exchange": False,
            "canonical_scale_quotient": False,
            "unit_bearing": True,
            "nonproxy_field_variable": True,
            "local_action_density": True,
            "variational_derivative_exported": True,
            "strict_coefficient_source_law": False,
        },
        {
            "template": "aggregate_norm_action_density",
            "formula": "||V||_1 / 5 * Phi^2 dmu",
            "receives_aggregate": True,
            "requires_KC_artifact_exchange": False,
            "canonical_scale_quotient": False,
            "unit_bearing": False,
            "nonproxy_field_variable": True,
            "local_action_density": True,
            "variational_derivative_exported": True,
            "strict_coefficient_source_law": False,
        },
        {
            "template": "formal_unit_token_pullback",
            "formula": "u_* V(x) Phi(x)^2 dVol_g",
            "receives_aggregate": True,
            "requires_KC_artifact_exchange": False,
            "canonical_scale_quotient": False,
            "unit_bearing": True,
            "nonproxy_field_variable": True,
            "local_action_density": True,
            "variational_derivative_exported": False,
            "strict_coefficient_source_law": False,
        },
        {
            "template": "KC_exchange_locked_coupling",
            "formula": "fixed(K<->C) coefficient times Phi^2",
            "receives_aggregate": True,
            "requires_KC_artifact_exchange": True,
            "canonical_scale_quotient": False,
            "unit_bearing": True,
            "nonproxy_field_variable": True,
            "local_action_density": True,
            "variational_derivative_exported": True,
            "strict_coefficient_source_law": False,
        },
        {
            "template": "strict_completed_coupling_schema",
            "formula": "c_N[V] Phi^2 dVol_N with sourced c_N and sourced unit quotient",
            "receives_aggregate": True,
            "requires_KC_artifact_exchange": False,
            "canonical_scale_quotient": True,
            "unit_bearing": True,
            "nonproxy_field_variable": True,
            "local_action_density": True,
            "variational_derivative_exported": True,
            "strict_coefficient_source_law": True,
        },
    ]
    rows = []
    for row in candidates:
        current_artifact_available = row["template"] != "strict_completed_coupling_schema"
        accepted_now = current_artifact_available and all([
            row["receives_aggregate"],
            not row["requires_KC_artifact_exchange"],
            row["canonical_scale_quotient"],
            row["unit_bearing"],
            row["nonproxy_field_variable"],
            row["local_action_density"],
            row["variational_derivative_exported"],
            row["strict_coefficient_source_law"],
        ])
        missing = [k for k in ["canonical_scale_quotient", "unit_bearing", "nonproxy_field_variable", "local_action_density", "variational_derivative_exported", "strict_coefficient_source_law"] if not row[k]]
        if row["requires_KC_artifact_exchange"]:
            missing.append("independence_from_KC_artifact_exchange")
        if not current_artifact_available:
            missing.append("current_artifact_export")
        rows.append({**row, "current_artifact_available": current_artifact_available, "missing_obligations": missing, "accepted_current_strict_coupling": accepted_now})
    return rows


def obligation_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {"obligation": "aggregate_reception_without_KC_exchange_representable", "satisfied": any(r["receives_aggregate"] and not r["requires_KC_artifact_exchange"] for r in rows if r["current_artifact_available"]), "evidence": "several templates can receive V=[1,2,2,2,2] without artifact exchange"},
        {"obligation": "canonical_scale_quotient_exported", "satisfied": any(r["canonical_scale_quotient"] and r["current_artifact_available"] for r in rows), "evidence": "only the unavailable completed schema contains this row"},
        {"obligation": "unit_bearing_nonproxy_density_exported", "satisfied": any(r["unit_bearing"] and r["nonproxy_field_variable"] and r["local_action_density"] and r["current_artifact_available"] for r in rows), "evidence": "formal templates can be dimensioned, but not canonically sourced"},
        {"obligation": "variational_derivative_exported", "satisfied": any(r["variational_derivative_exported"] and r["current_artifact_available"] for r in rows), "evidence": "some formal quadratic densities have Euler variations"},
        {"obligation": "strict_coefficient_source_law_exported", "satisfied": any(r["strict_coefficient_source_law"] and r["current_artifact_available"] for r in rows), "evidence": "no current artifact sources the coefficient/unit quotient"},
        {"obligation": "accepted_current_strict_coupling", "satisfied": any(r["accepted_current_strict_coupling"] for r in rows), "evidence": "no currently available row satisfies all obligations"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["receives_aggregate", "KC_exchange_independent", "canonical_scale_quotient", "unit_bearing", "nonproxy_local_density", "variational_derivative", "strict_coefficient_source"]
    full = (1 << len(names)) - 1
    return [{"mask": m, "present": {n: bool(m & (1 << i)) for i, n in enumerate(names)}, "accepts_strict_nonproxy_coupling": m == full} for m in range(1 << len(names))]


def build_payload(p2961: dict[str, Any], p2963: dict[str, Any]) -> dict[str, Any]:
    rows = coupling_rows()
    obligations = obligation_rows(rows)
    matrix = acceptance_matrix()
    return {
        "status": "P2964_UNIT_BEARING_NONPROXY_COUPLING_RECEPTION_NO_GO_NO_STRICT_EXPORT",
        "input_hashes": {"P2961": hashlib.sha256(P2961.read_bytes()).hexdigest() if P2961.exists() else None, "P2963": hashlib.sha256(P2963.read_bytes()).hexdigest() if P2963.exists() else None},
        "package_constants": {"aggregate_vector": AGGREGATE_VECTOR, "sum": sum(AGGREGATE_VECTOR), "eta": ETA, "delta": DELTA},
        "constructed_theoretical_objects": {"candidate_object": "AggregateReception_UnitBearingNonproxyCoupling_NoGo", "coupling_template_rows": rows, "proof_obligation_rows": obligations, "finite_acceptance_matrix": matrix},
        "coupling_certificate": {"template_count": len(rows), "currently_available_templates": [r["template"] for r in rows if r["current_artifact_available"]], "accepted_current_strict_couplings": [r["template"] for r in rows if r["accepted_current_strict_coupling"]], "strict_nonproxy_coupling_exported": any(r["accepted_current_strict_coupling"] for r in rows), "acceptance_matrix_rows": len(matrix), "accepted_rows": sum(1 for row in matrix if row["accepts_strict_nonproxy_coupling"])},
        "decision": {
            "positive_progress": "P2964 constructs the action-coupling reception interface explicitly: the P2938/P2961 aggregate can be received without assuming artifact-level K/C exchange.",
            "breakthrough": "No strict L_total coupling is exported: the current templates still lack a canonical scale quotient and strict coefficient/source law, while the only fully satisfying schema is a missing theorem object rather than a current artifact.",
            "negative_export_flags": {k: False for k in ["strict_unit_bearing_nonproxy_coupling_exported", "strict_ratio_package_source_exported", "damping_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "strict_p2938_provenance_exported", "strict_KC_exchange_source_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay formal quadratic densities, beta-scale normalization, scalar Euler insertion, K/C exchange, typed scalar mediators, target_sum cuts, or primitive equal-summand arguments.  The next proof-grade move must construct exactly one missing strict object: either a canonical scale/unit quotient with a coefficient source law for the aggregate-reception coupling, or a genuinely new typed structural object richer than the current K/C split; otherwise pivot outside the ratio-package lane while preserving the P2929-P2964 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["coupling_certificate"]
    lines = [
        "# P2964/S1914 unit-bearing nonproxy coupling reception no-go", "", f"Status: `{payload['status']}`", "", "## Coupling certificate", f"- template count: `{cert['template_count']}`", f"- currently available templates: `{cert['currently_available_templates']}`", f"- accepted current strict couplings: `{cert['accepted_current_strict_couplings']}`", f"- strict nonproxy coupling exported: `{cert['strict_nonproxy_coupling_exported']}`", f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`", "", "## Lay summary", payload["decision"]["positive_progress"], payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2961), read_json(P2963))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2964/S1914 unit-bearing nonproxy coupling reception no-go", "## P2964/S1914 unit-bearing nonproxy coupling reception no-go\n\n`P2964/S1914` attacks the post-P2963 admissible route of constructing a unit-bearing nonproxy action-density coupling receiving the P2938/P2961 aggregate without artifact-level K/C exchange.  Several formal templates receive `V=[1,2,2,2,2]`, and some have formal quadratic Euler variations, but current artifacts still do not export a canonical scale/unit quotient plus strict coefficient/source law.  The only fully satisfying row is an unavailable completed schema, so no strict ratio-package source, damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2964/S1914 aggregate-reception coupling `L_total` guard", "## P2964/S1914 aggregate-reception coupling `L_total` guard\n\n`P2964/S1914` constructs and audits aggregate-reception coupling templates.  Reception of the P2938/P2961 vector without K/C exchange is representable, but representation is not a sourced unit-bearing nonproxy `L_total` term: the canonical scale quotient and strict coefficient/source law are still missing.  Therefore no sourced damping coefficient enters `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE from P2964.\n")
    append_once(AGENTS, "Current unit-bearing nonproxy coupling reception no-go guardrail (P2964/S1914, 2026-06-20)", "## Current unit-bearing nonproxy coupling reception no-go guardrail (P2964/S1914, 2026-06-20)\n\n- P2964 constructs the post-P2963 aggregate-reception interface: the P2938/P2961 vector can be received by formal action-density templates without assuming artifact-level K/C exchange.\n- The route remains blocked because current artifacts still lack a canonical scale/unit quotient and a strict coefficient/source law; the only row satisfying all coupling obligations is an unavailable completed schema.\n- Do not promote formal quadratic densities, beta-scale normalization, scalar Euler insertion, K/C exchange, or typed scalar mediators to strict ratio-package source, damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- A next admissible move must construct exactly one missing strict object: a canonical scale/unit quotient with coefficient source law for the aggregate-reception coupling, or a genuinely new typed structural object richer than the current K/C split; otherwise pivot outside the ratio-package lane while preserving the P2929-P2964 boundary.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
