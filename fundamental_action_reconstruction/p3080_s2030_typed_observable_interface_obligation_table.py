#!/usr/bin/env python3
"""P3080/S2030: typed observable-interface obligation table.

P3079 showed that the current Z12 smoothing/dispersion branch does not source a
Lorentzian light cone, unit-normalized time, or spacetime wave EOM.  P3080 takes
the next honest non-selector step: build a finite typed interface table that
spells out exactly which additional objects would be needed before comparing the
internal nadsoliton/Z12 branch with standard theoretical physics.  The audit is
constructive but bounded: it creates typed interface candidates, obligation
gates, and a status lattice (sourced / formal-imported / absent), while blocking
observed-light, photon, empirical, L_total, bridge, selector, and ToE promotion.
"""
from __future__ import annotations

import hashlib, json, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3079_s2029_light_cone_causal_order_compatibility_audit import OUT as P3079

OUT = GEN / "p3080_s2030_typed_observable_interface_obligation_table.json"
MD = GEN / "p3080_s2030_typed_observable_interface_obligation_table.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "observable_interface": r"observable-interface|standard theoretical physics|dimension map|continuum limit|Lorentz signature|gauge representation|conserved current|empirical observable",
    "z12_smoothing_branch": r"Z12|Dirichlet|Laplacian|diffusion|smoothing|dispersion",
    "blocked_promotions": r"observed light|gauge photons|spacetime EOM|empirical physics|L_total|ToE|selector closure",
}

INTERFACE_OBJECTS = (
    {
        "id": "z12_dirichlet_energy_scalar",
        "typed_domain": "rho: Z12 -> R",
        "typed_codomain": "R_nonnegative_internal",
        "internal_source": True,
        "formula": "E_D(rho)=1/2 sum_i (rho_{i+1}-rho_i)^2",
    },
    {
        "id": "z12_laplacian_smoothing_flow",
        "typed_domain": "rho: Z12 -> R",
        "typed_codomain": "tangent vector in R^12 / constants",
        "internal_source": True,
        "formula": "dot(rho)=-L_Z12 rho",
    },
    {
        "id": "formal_spectral_dispersion_proxy",
        "typed_domain": "mode j in Z12^",
        "typed_codomain": "dimensionless lambda_j and sqrt(lambda_j)",
        "internal_source": True,
        "formula": "lambda_j=2-2 cos(2*pi*j/12)",
    },
    {
        "id": "imported_standard_physics_readout_template",
        "typed_domain": "external spacetime/gauge/measurement data",
        "typed_codomain": "standard-physics observable slots",
        "internal_source": False,
        "formula": "dimensionful spacetime, Lorentz metric, gauge bundle, detector map (imported)",
    },
)

OBLIGATIONS = (
    {
        "id": "dimension_map",
        "question": "Does the object export units mapping internal scalars to length/time/action/energy dimensions?",
        "required_for_standard_physics": True,
    },
    {
        "id": "continuum_limit",
        "question": "Does the object export a controlled Z12-to-continuum limiting functor with error bounds?",
        "required_for_standard_physics": True,
    },
    {
        "id": "lorentz_signature",
        "question": "Does the object export a non-imported Lorentzian metric/signature and finite physical cone?",
        "required_for_standard_physics": True,
    },
    {
        "id": "gauge_representation",
        "question": "Does the object export a gauge bundle/connection/representation comparable to photon or gauge fields?",
        "required_for_standard_physics": True,
    },
    {
        "id": "conserved_current",
        "question": "Does the object export a nontrivial sourced Noether/conservation current with units?",
        "required_for_standard_physics": True,
    },
    {
        "id": "empirical_observable",
        "question": "Does the object export a detector/readout map to empirical observables with calibrated units?",
        "required_for_standard_physics": True,
    },
)

STATUS_DETAIL = {
    ("z12_dirichlet_energy_scalar", "dimension_map"): ("absent", "positive internal scalar but no length/action/energy unit source"),
    ("z12_dirichlet_energy_scalar", "continuum_limit"): ("formal_imported", "quadratic Dirichlet form has a familiar continuum analogy only after importing lattice spacing and limit"),
    ("z12_dirichlet_energy_scalar", "lorentz_signature"): ("absent", "positive energy scalar has no metric signature or time coordinate"),
    ("z12_dirichlet_energy_scalar", "gauge_representation"): ("absent", "no bundle, connection, representation, or photon sector is exported"),
    ("z12_dirichlet_energy_scalar", "conserved_current"): ("absent", "no unit-bearing nontrivial Noether current follows from a static scalar"),
    ("z12_dirichlet_energy_scalar", "empirical_observable"): ("absent", "no calibrated detector/readout map is exported"),
    ("z12_laplacian_smoothing_flow", "dimension_map"): ("absent", "flow parameter is internal and not a calibrated physical time"),
    ("z12_laplacian_smoothing_flow", "continuum_limit"): ("formal_imported", "heat-equation analogy needs imported lattice spacing, time unit, and limiting procedure"),
    ("z12_laplacian_smoothing_flow", "lorentz_signature"): ("absent", "parabolic smoothing is not a Lorentzian wave equation"),
    ("z12_laplacian_smoothing_flow", "gauge_representation"): ("absent", "scalar smoothing has no gauge connection or field strength"),
    ("z12_laplacian_smoothing_flow", "conserved_current"): ("formal_imported", "total preservation is real internally but not a unit-bearing physical Noether current"),
    ("z12_laplacian_smoothing_flow", "empirical_observable"): ("absent", "no detector/readout calibration is supplied"),
    ("formal_spectral_dispersion_proxy", "dimension_map"): ("absent", "lambda_j and sqrt(lambda_j) are dimensionless modal diagnostics"),
    ("formal_spectral_dispersion_proxy", "continuum_limit"): ("formal_imported", "small-k lambda~k^2 is a formal analogy without sourced lattice spacing"),
    ("formal_spectral_dispersion_proxy", "lorentz_signature"): ("absent", "dispersion proxy lacks a sourced time and metric sign"),
    ("formal_spectral_dispersion_proxy", "gauge_representation"): ("absent", "Fourier modes are not gauge photons without a gauge representation theorem"),
    ("formal_spectral_dispersion_proxy", "conserved_current"): ("absent", "modal labels do not provide a sourced current"),
    ("formal_spectral_dispersion_proxy", "empirical_observable"): ("absent", "no measured frequency/wavelength/readout unit is exported"),
    ("imported_standard_physics_readout_template", "dimension_map"): ("formal_imported", "standard units can be written only by external convention"),
    ("imported_standard_physics_readout_template", "continuum_limit"): ("formal_imported", "continuum manifold is assumed rather than derived"),
    ("imported_standard_physics_readout_template", "lorentz_signature"): ("formal_imported", "Lorentz signature is supplied by the template, not by Z12 data"),
    ("imported_standard_physics_readout_template", "gauge_representation"): ("formal_imported", "gauge representation is supplied by the template"),
    ("imported_standard_physics_readout_template", "conserved_current"): ("formal_imported", "standard currents can be named only after importing a physical action/symmetry"),
    ("imported_standard_physics_readout_template", "empirical_observable"): ("formal_imported", "detector calibration is external"),
}


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def obligation_rows() -> list[dict[str, Any]]:
    rows = []
    for obj in INTERFACE_OBJECTS:
        for obligation in OBLIGATIONS:
            status, detail = STATUS_DETAIL[(obj["id"], obligation["id"])]
            rows.append({
                "interface_object": obj["id"],
                "obligation": obligation["id"],
                "status": status,
                "strictly_sourced": status == "sourced",
                "formal_or_imported": status == "formal_imported",
                "absent": status == "absent",
                "detail": detail,
            })
    return rows


def object_rows(obligations: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for obj in INTERFACE_OBJECTS:
        subset = [row for row in obligations if row["interface_object"] == obj["id"]]
        rows.append({
            **obj,
            "sourced_obligations": sum(1 for row in subset if row["strictly_sourced"]),
            "formal_imported_obligations": sum(1 for row in subset if row["formal_or_imported"]),
            "absent_obligations": sum(1 for row in subset if row["absent"]),
            "standard_physics_interface_accepted": all(row["strictly_sourced"] for row in subset),
        })
    return rows


def obligation_summary(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for obligation in OBLIGATIONS:
        subset = [row for row in rows if row["obligation"] == obligation["id"]]
        out.append({
            "obligation": obligation["id"],
            "sourced_rows": sum(1 for row in subset if row["strictly_sourced"]),
            "formal_imported_rows": sum(1 for row in subset if row["formal_or_imported"]),
            "absent_rows": sum(1 for row in subset if row["absent"]),
            "discharged_by_current_artifacts": any(row["strictly_sourced"] for row in subset),
        })
    return out


def build_payload() -> dict[str, Any]:
    p3079 = read_json(P3079)
    greps = content_grep()
    obligations = obligation_rows()
    objects = object_rows(obligations)
    summary = obligation_summary(obligations)
    proof_obligations = [
        {"obligation": "read_p3079_boundary", "satisfied": True, "detail": "P3079 accepted internal causal-order sources remain zero"},
        {"obligation": "construct_typed_interface_objects", "satisfied": True, "detail": "4 typed objects separate internal Z12 data from imported standard-physics template"},
        {"obligation": "construct_standard_physics_obligation_basis", "satisfied": True, "detail": "6 required comparison obligations are explicit"},
        {"obligation": "classify_each_cell_as_sourced_formal_or_absent", "satisfied": True, "detail": "24 object-obligation cells are classified"},
        {"obligation": "export_standard_physics_interface", "satisfied": False, "detail": "0 cells are strictly sourced and 0 objects satisfy every obligation"},
    ]
    return {
        "status": "P3080_TYPED_OBSERVABLE_INTERFACE_OBLIGATION_BOUNDED_NO_GO",
        "input_hashes": {"P3079": hashlib.sha256(P3079.read_bytes()).hexdigest() if P3079.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": greps,
            "interface_audit_object": {
                "object": "Z12TypedObservableInterfaceObligationTable",
                "source_reused": "P3079 smoothing/dispersion branch with causal-order promotion frozen",
                "interface_objects": [obj["id"] for obj in INTERFACE_OBJECTS],
                "obligations": [ob["id"] for ob in OBLIGATIONS],
                "acceptance_predicate": "one interface object must strictly source all six standard-physics comparison obligations without imported spacetime/gauge/measurement premises",
            },
            "interface_object_rows": objects,
            "obligation_cell_rows": obligations,
            "obligation_summary_rows": summary,
        },
        "finite_certificate": {
            "content_grep_lanes": len(greps),
            "content_grep_hits": sum(row["hit_count"] for row in greps),
            "p3079_accepted_internal_causal_order_sources": p3079.get("finite_certificate", {}).get("accepted_internal_causal_order_sources"),
            "interface_objects": len(INTERFACE_OBJECTS),
            "standard_physics_obligations": len(OBLIGATIONS),
            "obligation_cell_rows": len(obligations),
            "strictly_sourced_cells": sum(1 for row in obligations if row["strictly_sourced"]),
            "formal_imported_cells": sum(1 for row in obligations if row["formal_or_imported"]),
            "absent_cells": sum(1 for row in obligations if row["absent"]),
            "accepted_standard_physics_interface_objects": sum(1 for row in objects if row["standard_physics_interface_accepted"]),
            "fully_discharged_obligations": sum(1 for row in summary if row["discharged_by_current_artifacts"]),
            "proof_obligations": len(proof_obligations),
            "satisfied_proof_obligations": sum(1 for row in proof_obligations if row["satisfied"]),
        },
        "proof_obligations": proof_obligations,
        "decision": {
            "bounded_result": "P3080 constructs a typed observable-interface obligation table between the internal Z12 smoothing/dispersion branch and standard theoretical physics.  The internal Dirichlet energy, Laplacian smoothing flow, and spectral dispersion proxy are real typed mathematical objects, but every standard-physics comparison obligation is either formal/imported or absent.  No current object strictly sources units, a continuum limit, Lorentz signature, gauge representation, unit-bearing conserved current, or empirical readout together.",
            "negative_export_flags": {key: False for key in ["dimension_map_exported", "continuum_limit_exported", "lorentz_signature_exported", "gauge_representation_exported", "unit_bearing_conserved_current_exported", "empirical_observable_exported", "observed_light_exported", "gauge_photon_sector_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"typed_interface_table_constructed": True, "obligation_status_lattice_constructed": True, "standard_physics_comparison_boundary_made_explicit": True},
            "next_honest_step": "Attack exactly one missing interface atom rather than replaying selector or light-cone promotion: construct a bounded dimension/action-unit source audit for the internal Dirichlet energy scalar, testing whether any current nadsoliton/Z12 artifact provides a non-imported unit map from E_D(rho) to action/energy/time units.  If no such source is found, keep standard-physics comparison at formal-interface status and pivot to a different typed object.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3080/S2030 typed observable-interface obligation table", "", f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- P3079 accepted internal causal-order sources: `{c['p3079_accepted_internal_causal_order_sources']}`",
        f"- interface objects: `{c['interface_objects']}`",
        f"- standard-physics obligations: `{c['standard_physics_obligations']}`",
        f"- obligation cell rows: `{c['obligation_cell_rows']}`",
        f"- strictly sourced cells: `{c['strictly_sourced_cells']}`",
        f"- formal/imported cells: `{c['formal_imported_cells']}`",
        f"- absent cells: `{c['absent_cells']}`",
        f"- accepted standard-physics interface objects: `{c['accepted_standard_physics_interface_objects']}`",
        f"- fully discharged obligations: `{c['fully_discharged_obligations']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3080/S2030 typed observable-interface obligation table", "## P3080/S2030 typed observable-interface obligation table\n\n`P3080/S2030` constructs `Z12TypedObservableInterfaceObligationTable` after the `P3079` causal-order no-go.  It separates four typed interface objects (Dirichlet energy scalar, Laplacian smoothing flow, formal spectral dispersion proxy, and imported standard-physics readout template) across six obligations: dimension map, continuum limit, Lorentz signature, gauge representation, conserved current, and empirical observable.  The `24` obligation cells contain `0` strictly sourced cells, `10` formal/imported cells, and `14` absent cells, so no standard-physics interface object is accepted.  No observed light, gauge photons, spacetime EOM, empirical physics, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3080/S2030 observable-interface remains formal", "## P3080/S2030 observable-interface remains formal\n\n`P3080/S2030` makes the internal-to-standard-physics interface obligations explicit.  The Dirichlet/Laplacian branch remains a real internal variational/smoothing structure, but it still lacks a non-imported dimension map, controlled continuum functor with units, Lorentzian signature, gauge representation, unit-bearing conserved current, and empirical readout.  Therefore it adds no physical `L_total` term and no spacetime/gauge EOM.\n")
    append_once(AGENTS, "Current typed observable-interface guardrail (P3080/S2030, 2026-06-24)", "## Current typed observable-interface guardrail (P3080/S2030, 2026-06-24)\n\n- P3080 follows the P3079 recommendation and constructs a bounded typed observable-interface obligation table for comparing the Z12 smoothing/dispersion branch with standard theoretical physics.\n- The finite audit has `4` typed interface objects, `6` standard-physics obligations, and `24` obligation cells: `0` strictly sourced, `10` formal/imported, and `14` absent; no object discharges the full interface predicate.\n- Do not promote Dirichlet energy, Laplacian smoothing, spectral dispersion proxies, or imported standard-physics templates to observed light, gauge photons, spacetime EOM, empirical physics, `QW-2191` discharge, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one missing interface atom: a bounded dimension/action-unit source audit for the internal Dirichlet energy scalar, unless a different genuinely new typed object is introduced.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
