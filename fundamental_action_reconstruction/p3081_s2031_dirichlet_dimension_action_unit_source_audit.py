#!/usr/bin/env python3
"""P3081/S2031: Dirichlet dimension/action-unit source audit.

P3080 made the standard-physics interface obligations explicit and selected one
honest missing atom: a non-imported dimension/action-unit source for the internal
Dirichlet energy scalar E_D(rho).  P3081 constructs that atom as a finite audit.
It computes an exact binary-profile Dirichlet-energy spectrum on Z12, builds a
candidate unit-source matrix, and checks whether any current nadsoliton/Z12 datum
exports a dimensionful action/energy/time unit without target normalization,
external hbar/c/lattice spacing, selector replay, bridge transfer, or observed-
physics promotion.
"""
from __future__ import annotations

import hashlib, itertools, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3080_s2030_typed_observable_interface_obligation_table import OUT as P3080

OUT = GEN / "p3081_s2031_dirichlet_dimension_action_unit_source_audit.json"
MD = GEN / "p3081_s2031_dirichlet_dimension_action_unit_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
Z12_SIZE = 12

CONTENT_PATTERNS = {
    "unit_source": r"dimension/action-unit|action-unit|dimension map|unit-bearing|canonical unit|hbar|lattice spacing",
    "dirichlet_branch": r"Dirichlet energy|E_D\(rho\)|Laplacian|Z12|smoothing",
    "blocked_promotions": r"observed light|gauge photons|spacetime EOM|empirical physics|L_total|ToE|selector closure|bridge/role-transfer",
}

UNIT_CANDIDATES = (
    {
        "id": "binary_dirichlet_minimum_gap",
        "description": "minimum nonzero E_D over binary Z12 profiles",
        "internal_artifact": True,
        "positive_nonzero_reference": True,
        "profile_independent": True,
        "scale_invariant_under_rho_rescaling": False,
        "dimensionful_action_unit_exported": False,
        "independent_physical_calibration": False,
        "blocker": "the gap is an internal dimensionless profile-energy quantum and rescales with rho-amplitude conventions",
    },
    {
        "id": "cycle_cardinality_twelve",
        "description": "|Z12|=12 and its reciprocal normalization",
        "internal_artifact": True,
        "positive_nonzero_reference": True,
        "profile_independent": True,
        "scale_invariant_under_rho_rescaling": True,
        "dimensionful_action_unit_exported": False,
        "independent_physical_calibration": False,
        "blocker": "cardinality is a dimensionless count, not an action/energy/time unit",
    },
    {
        "id": "alpha_geo_information_amplitude",
        "description": "alpha_geo=4 ln 2 as informational amplitude",
        "internal_artifact": True,
        "positive_nonzero_reference": True,
        "profile_independent": True,
        "scale_invariant_under_rho_rescaling": True,
        "dimensionful_action_unit_exported": False,
        "independent_physical_calibration": False,
        "blocker": "alpha_geo is dimensionless and role transfer to physical constants is not licensed",
    },
    {
        "id": "p2663_one_bit_entropy_reference",
        "description": "one-bit entropy/reference-cell candidate from earlier UV-unit lane",
        "internal_artifact": True,
        "positive_nonzero_reference": True,
        "profile_independent": True,
        "scale_invariant_under_rho_rescaling": True,
        "dimensionful_action_unit_exported": False,
        "independent_physical_calibration": False,
        "blocker": "bit count remains dimensionless; P2689/P2690 left bit-to-length/action and source-localizer premises missing",
    },
    {
        "id": "spectral_gap_lambda_one",
        "description": "first nonzero Z12 Laplacian eigenvalue lambda_1",
        "internal_artifact": True,
        "positive_nonzero_reference": True,
        "profile_independent": True,
        "scale_invariant_under_rho_rescaling": True,
        "dimensionful_action_unit_exported": False,
        "independent_physical_calibration": False,
        "blocker": "lambda_1 is dimensionless without sourced lattice spacing/time normalization",
    },
    {
        "id": "imported_hbar_c_lattice_spacing_template",
        "description": "external hbar/c/a template for physical units",
        "internal_artifact": False,
        "positive_nonzero_reference": True,
        "profile_independent": True,
        "scale_invariant_under_rho_rescaling": True,
        "dimensionful_action_unit_exported": True,
        "independent_physical_calibration": True,
        "blocker": "passes only by importing empirical dimensional constants and lattice spacing",
    },
)

REQUIRED_GATES = (
    "internal_artifact",
    "positive_nonzero_reference",
    "profile_independent",
    "scale_invariant_under_rho_rescaling",
    "dimensionful_action_unit_exported",
    "independent_physical_calibration",
)


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def dirichlet_energy(profile: tuple[int, ...]) -> float:
    return 0.5 * sum((profile[(i + 1) % Z12_SIZE] - profile[i]) ** 2 for i in range(Z12_SIZE))


def binary_energy_rows() -> list[dict[str, Any]]:
    rows = []
    for profile in itertools.product((0, 1), repeat=Z12_SIZE):
        energy = dirichlet_energy(profile)
        boundary_edges = int(round(2 * energy))
        rows.append({
            "profile_bits": "".join(str(bit) for bit in profile),
            "dirichlet_energy": energy,
            "boundary_edges": boundary_edges,
            "constant_profile": boundary_edges == 0,
            "nonzero_internal_reference": boundary_edges > 0,
        })
    return rows


def energy_spectrum(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    values = sorted({row["boundary_edges"] for row in rows})
    return [{
        "boundary_edges": value,
        "dirichlet_energy": value / 2,
        "profile_count": sum(1 for row in rows if row["boundary_edges"] == value),
    } for value in values]


def gate_rows() -> list[dict[str, Any]]:
    rows = []
    for candidate in UNIT_CANDIDATES:
        for gate in REQUIRED_GATES:
            passed = bool(candidate[gate])
            rows.append({
                "candidate": candidate["id"],
                "required_gate": gate,
                "gate_passed": passed,
                "blocking_if_failed": not passed,
                "detail": "passed" if passed else candidate["blocker"],
            })
    return rows


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for candidate in UNIT_CANDIDATES:
        subset = [row for row in gates if row["candidate"] == candidate["id"]]
        out.append({
            "candidate": candidate["id"],
            "passed_gates": sum(1 for row in subset if row["gate_passed"]),
            "failed_gates": sum(1 for row in subset if not row["gate_passed"]),
            "accepted_nonimported_dimension_action_unit_source": all(row["gate_passed"] for row in subset),
        })
    return out


def build_payload() -> dict[str, Any]:
    p3080 = read_json(P3080)
    greps = content_grep()
    binary_rows = binary_energy_rows()
    spectrum = energy_spectrum(binary_rows)
    gates = gate_rows()
    aggs = aggregate(gates)
    accepted = [row for row in aggs if row["accepted_nonimported_dimension_action_unit_source"]]
    nonconstant = [row for row in binary_rows if not row["constant_profile"]]
    nonzero_energies = [row["dirichlet_energy"] for row in nonconstant]
    proof_obligations = [
        {"obligation": "read_p3080_interface_atom", "satisfied": True, "detail": "P3080 selected dimension/action-unit source as exactly one missing interface atom"},
        {"obligation": "compute_exact_binary_dirichlet_spectrum", "satisfied": True, "detail": "all 2^12 binary profiles are evaluated exactly"},
        {"obligation": "construct_candidate_unit_source_matrix", "satisfied": True, "detail": "6 candidate source classes x 6 required gates = 36 rows"},
        {"obligation": "separate_internal_dimensionless_witnesses_from_imported_units", "satisfied": True, "detail": "internal positive references are dimensionless; imported hbar/c/a is not internal"},
        {"obligation": "export_nonimported_dimension_action_unit_source", "satisfied": False, "detail": "0 candidates pass all six gates"},
    ]
    return {
        "status": "P3081_DIRICHLET_DIMENSION_ACTION_UNIT_SOURCE_BOUNDED_NO_GO",
        "input_hashes": {"P3080": hashlib.sha256(P3080.read_bytes()).hexdigest() if P3080.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": greps,
            "unit_source_audit_object": {
                "object": "Z12DirichletDimensionActionUnitSourceAudit",
                "source_reused": "P3080 missing dimension/action-unit interface atom for E_D(rho)",
                "required_gates": list(REQUIRED_GATES),
                "candidate_sources": [candidate["id"] for candidate in UNIT_CANDIDATES],
                "acceptance_predicate": "internal positive profile-independent scale-stable source that exports dimensionful action/energy/time units with independent calibration",
            },
            "binary_dirichlet_energy_rows": binary_rows,
            "binary_dirichlet_energy_spectrum": spectrum,
            "candidate_gate_rows": gates,
            "candidate_aggregate_certificate": aggs,
        },
        "finite_certificate": {
            "content_grep_lanes": len(greps),
            "content_grep_hits": sum(row["hit_count"] for row in greps),
            "p3080_accepted_standard_physics_interface_objects": p3080.get("finite_certificate", {}).get("accepted_standard_physics_interface_objects"),
            "binary_profile_rows": len(binary_rows),
            "constant_profile_rows": sum(1 for row in binary_rows if row["constant_profile"]),
            "nonconstant_profile_rows": len(nonconstant),
            "energy_spectrum_rows": len(spectrum),
            "minimum_nonzero_dirichlet_energy": min(nonzero_energies),
            "maximum_dirichlet_energy": max(row["dirichlet_energy"] for row in binary_rows),
            "unit_source_candidates": len(UNIT_CANDIDATES),
            "required_gates": len(REQUIRED_GATES),
            "candidate_gate_rows": len(gates),
            "accepted_nonimported_dimension_action_unit_sources": len(accepted),
            "proof_obligations": len(proof_obligations),
            "satisfied_proof_obligations": sum(1 for row in proof_obligations if row["satisfied"]),
        },
        "proof_obligations": proof_obligations,
        "decision": {
            "bounded_result": "P3081 constructs the requested dimension/action-unit source audit for the internal Dirichlet energy scalar.  The exact 2^12 binary-profile spectrum gives real positive internal reference numbers, including a minimum nonzero E_D=1.0 and maximum E_D=6.0, but these are dimensionless profile-energy witnesses.  Internal counts, alpha_geo, entropy-bit references, and spectral gaps do not export action/energy/time dimensions; hbar/c/lattice-spacing templates pass only by external import.  Therefore no non-imported dimension/action-unit source is exported.",
            "negative_export_flags": {key: False for key in ["dimension_action_unit_source_exported", "energy_unit_source_exported", "time_unit_source_exported", "action_quantum_exported", "physical_hamiltonian_exported", "observed_light_exported", "gauge_photon_sector_exported", "spacetime_eom_exported", "empirical_physics_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"binary_dirichlet_spectrum_computed": True, "positive_internal_reference_witnesses_found": True, "unit_source_gate_matrix_executed": True},
            "next_honest_step": "Pivot to exactly one different standard-physics interface atom rather than replaying units: construct a bounded continuum-limit functor obstruction/witness audit for the Z12 Dirichlet/Laplacian branch, testing whether any current artifact supplies a non-imported lattice-spacing/refinement map and error-controlled Z12-to-continuum passage.  Keep Lorentz/gauge/observed-physics promotion blocked unless that single functor atom is actually sourced.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3081/S2031 Dirichlet dimension/action-unit source audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- P3080 accepted standard-physics interface objects: `{c['p3080_accepted_standard_physics_interface_objects']}`",
        f"- binary profile rows: `{c['binary_profile_rows']}`",
        f"- constant profile rows: `{c['constant_profile_rows']}`",
        f"- nonconstant profile rows: `{c['nonconstant_profile_rows']}`",
        f"- energy spectrum rows: `{c['energy_spectrum_rows']}`",
        f"- minimum nonzero Dirichlet energy: `{c['minimum_nonzero_dirichlet_energy']}`",
        f"- maximum Dirichlet energy: `{c['maximum_dirichlet_energy']}`",
        f"- unit source candidates: `{c['unit_source_candidates']}`",
        f"- candidate gate rows: `{c['candidate_gate_rows']}`",
        f"- accepted non-imported dimension/action-unit sources: `{c['accepted_nonimported_dimension_action_unit_sources']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3081/S2031 Dirichlet dimension/action-unit source audit", "## P3081/S2031 Dirichlet dimension/action-unit source audit\n\n`P3081/S2031` attacks exactly one P3080 interface atom: a non-imported dimension/action-unit source for `E_D(rho)`.  It computes all `2^12 = 4096` binary Z12 profiles, their exact Dirichlet-energy spectrum, and a `6 x 6 = 36` candidate unit-source gate matrix.  The binary spectrum supplies real positive internal dimensionless witnesses (`min E_D=1`, `max E_D=6`), but no current internal candidate exports action/energy/time dimensions with independent calibration; the `hbar/c/a` template succeeds only by import.  No physical Hamiltonian, spacetime EOM, observed light, gauge photons, empirical physics, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3081/S2031 Dirichlet units remain unsourced", "## P3081/S2031 Dirichlet units remain unsourced\n\n`P3081/S2031` strengthens the Dirichlet/Laplacian audit by computing an exact finite energy spectrum for binary Z12 profiles.  The result gives internal positive dimensionless reference values, but it still does not source an action unit, energy unit, time unit, physical Hamiltonian normalization, or physical `L_total` term.\n")
    append_once(AGENTS, "Current Dirichlet dimension/action-unit source guardrail (P3081/S2031, 2026-06-24)", "## Current Dirichlet dimension/action-unit source guardrail (P3081/S2031, 2026-06-24)\n\n- P3081 follows the P3080 recommendation and audits one missing interface atom: a non-imported dimension/action-unit source for the internal Dirichlet energy scalar `E_D(rho)`.\n- The finite audit computes all `4096` binary Z12 profiles, `7` energy-spectrum rows, and `36` candidate unit-source gate rows; positive internal dimensionless witnesses exist, but `0` candidates export a non-imported action/energy/time unit.\n- Do not promote Dirichlet energy gaps, Z12 counts, `alpha_geo`, entropy-bit references, spectral gaps, or imported `hbar/c/a` templates to physical Hamiltonian, spacetime EOM, observed light, gauge photons, empirical physics, `QW-2191` discharge, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one different interface atom: a bounded continuum-limit functor obstruction/witness audit for the Z12 Dirichlet/Laplacian branch, unless a genuinely new unit-source theorem is introduced.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
