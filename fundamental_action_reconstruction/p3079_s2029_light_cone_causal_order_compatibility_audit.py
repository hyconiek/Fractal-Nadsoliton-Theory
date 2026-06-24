#!/usr/bin/env python3
"""P3079/S2029: light-cone/causal-order compatibility audit.

P3078 froze the Dirichlet-to-wave promotion because the current Z12 artifacts do
not source an intrinsic momentum/symplectic phase space.  P3079 pivots to the
next non-selector typed object requested by that boundary: a bounded audit of
whether the same internal Z12 dispersion data already carries a metric
signature, finite propagation cone, or unit-normalized causal order.  The audit
constructs finite graph-distance, diffusion-kernel support, adjacency-cone, and
formal imported Lorentz candidates, then blocks any false promotion to observed
light, spacetime EOM, gauge photons, or empirical physics.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3078_s2028_intrinsic_momentum_symplectic_source_audit import OUT as P3078

OUT = GEN / "p3079_s2029_light_cone_causal_order_compatibility_audit.json"
MD = GEN / "p3079_s2029_light_cone_causal_order_compatibility_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
Z12_SIZE = 12
DIFFUSION_STEPS = (1 / 12, 1 / 6, 1 / 4)
CONTENT_PATTERNS = {
    "causal_language": r"light-cone|causal order|causal cone|finite propagation|metric signature",
    "z12_dispersion": r"Z12|Dirichlet|Laplacian|lambda_j|dispersion",
    "time_units": r"unit-bearing time|unit-normalized|Lorentzian|spacetime EOM",
    "no_physics_promotion": r"observed light|gauge photons|empirical physics|L_total|ToE|selector closure",
}

CAUSAL_CANDIDATES = (
    {
        "id": "cycle_graph_distance_metric",
        "description": "undirected shortest-path metric d_Z12(a,b)",
        "internal_z12_object": True,
        "metric_signature_exported": False,
        "finite_propagation_cone_exported": False,
        "unit_normalized_time_exported": False,
        "wave_eom_exported": False,
        "blocker": "positive graph distance is spatial/combinatorial only; no Lorentzian sign, clock, or EOM",
    },
    {
        "id": "heat_kernel_diffusion_support",
        "description": "exp(-t L) support for positive fractional Dirichlet steps",
        "internal_z12_object": True,
        "metric_signature_exported": False,
        "finite_propagation_cone_exported": False,
        "unit_normalized_time_exported": False,
        "wave_eom_exported": False,
        "blocker": "diffusion has immediate nonzero tails on the connected finite cycle, not a finite light cone",
    },
    {
        "id": "adjacency_iterate_combinatorial_cone",
        "description": "support reachable by at most n nearest-neighbor adjacency steps",
        "internal_z12_object": True,
        "metric_signature_exported": False,
        "finite_propagation_cone_exported": True,
        "unit_normalized_time_exported": False,
        "wave_eom_exported": False,
        "blocker": "finite graph reachability is a combinatorial propagation bound, not a unit-normalized causal order or Lorentzian EOM",
    },
    {
        "id": "dirichlet_spectral_velocity_proxy",
        "description": "formal group-velocity proxy from lambda_j=2-2 cos(k_j)",
        "internal_z12_object": True,
        "metric_signature_exported": False,
        "finite_propagation_cone_exported": False,
        "unit_normalized_time_exported": False,
        "wave_eom_exported": False,
        "blocker": "spectral slope is a dimensionless dispersion diagnostic without sourced time, units, or Lorentz signature",
    },
    {
        "id": "imported_minkowski_light_cone_ansatz",
        "description": "external ds^2=-dt^2+dx^2 light-cone template",
        "internal_z12_object": False,
        "metric_signature_exported": True,
        "finite_propagation_cone_exported": True,
        "unit_normalized_time_exported": True,
        "wave_eom_exported": True,
        "blocker": "works only by importing spacetime geometry and wave EOM rather than deriving them from current nadsoliton/Z12 artifacts",
    },
)

REQUIRED_GATES = (
    "internal_z12_object",
    "metric_signature_exported",
    "finite_propagation_cone_exported",
    "unit_normalized_time_exported",
    "wave_eom_exported",
)


def cycle_distance(a: int, b: int) -> int:
    delta = abs(a - b) % Z12_SIZE
    return min(delta, Z12_SIZE - delta)


def laplacian_eigenvalue(j: int) -> float:
    return 2.0 - 2.0 * math.cos(2.0 * math.pi * j / Z12_SIZE)


def heat_kernel_value(distance: int, step: float) -> float:
    return sum(math.exp(-step * laplacian_eigenvalue(j)) * math.cos(2.0 * math.pi * j * distance / Z12_SIZE) for j in range(Z12_SIZE)) / Z12_SIZE


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def distance_rows() -> list[dict[str, Any]]:
    return [{"source_node": a, "target_node": b, "cycle_distance": cycle_distance(a, b), "within_one_step_cone": cycle_distance(a, b) <= 1} for a in range(Z12_SIZE) for b in range(Z12_SIZE)]


def shell_rows(distances: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for d in range(0, 7):
        rows.append({"distance_shell": d, "ordered_pair_count": sum(1 for row in distances if row["cycle_distance"] == d)})
    return rows


def diffusion_rows() -> list[dict[str, Any]]:
    rows = []
    for step in DIFFUSION_STEPS:
        for d in range(0, 7):
            val = heat_kernel_value(d, step)
            rows.append({
                "fractional_step": f"{step:.12g}",
                "distance_shell": d,
                "heat_kernel_value": val,
                "nonzero_tail": abs(val) > 1e-12,
                "outside_unit_cone": d > 1,
                "violates_finite_unit_cone_if_read_as_light": d > 1 and abs(val) > 1e-12,
            })
    return rows


def spectral_velocity_rows() -> list[dict[str, Any]]:
    rows = []
    for j in range(Z12_SIZE):
        k = 2.0 * math.pi * j / Z12_SIZE
        lam = laplacian_eigenvalue(j)
        omega = math.sqrt(lam)
        # For the formal wave proxy omega=sqrt(lambda(k)), velocity is dimensionless
        # and singular at the constant mode; it is not a sourced light speed.
        velocity = None if j == 0 else math.sin(k) / omega
        rows.append({
            "mode_j": j,
            "k_label": f"2*pi*{j}/12",
            "lambda_value": lam,
            "formal_omega_sqrt_lambda": omega,
            "dimensionless_group_velocity_proxy": velocity,
            "unit_light_speed_certified": False,
        })
    return rows


def gate_rows() -> list[dict[str, Any]]:
    rows = []
    for candidate in CAUSAL_CANDIDATES:
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
    for candidate in CAUSAL_CANDIDATES:
        subset = [row for row in gates if row["candidate"] == candidate["id"]]
        out.append({
            "candidate": candidate["id"],
            "passed_gates": sum(1 for row in subset if row["gate_passed"]),
            "failed_gates": sum(1 for row in subset if not row["gate_passed"]),
            "accepted_internal_causal_order_source": all(row["gate_passed"] for row in subset),
        })
    return out


def build_payload() -> dict[str, Any]:
    p3078 = read_json(P3078)
    greps = content_grep()
    distances = distance_rows()
    shells = shell_rows(distances)
    diffusion = diffusion_rows()
    velocities = spectral_velocity_rows()
    gates = gate_rows()
    aggs = aggregate(gates)
    accepted = [row for row in aggs if row["accepted_internal_causal_order_source"]]
    proof_obligations = [
        {"obligation": "content_first_grep_for_causal_order_sources", "satisfied": True, "detail": "searched causal, Z12-dispersion, time-unit, and no-promotion lanes"},
        {"obligation": "construct_cycle_distance_matrix", "satisfied": True, "detail": "12 x 12 ordered-pair graph-distance matrix"},
        {"obligation": "test_diffusion_against_finite_unit_cone", "satisfied": True, "detail": "3 fractional steps x 7 distance shells show nonzero outside-unit-cone tails"},
        {"obligation": "test_dimensionless_spectral_velocity_proxy", "satisfied": True, "detail": "12 modal velocity rows remain dimensionless and not unit-normalized"},
        {"obligation": "execute_causal_candidate_gate_matrix", "satisfied": True, "detail": "5 candidates x 5 required gates = 25 exact gate rows"},
        {"obligation": "export_accepted_internal_causal_order_source", "satisfied": False, "detail": "no candidate supplies internal metric signature, finite cone, unit time, and wave EOM together"},
    ]
    return {
        "status": "P3079_LIGHT_CONE_CAUSAL_ORDER_COMPATIBILITY_BOUNDED_NO_GO",
        "input_hashes": {"P3078": hashlib.sha256(P3078.read_bytes()).hexdigest() if P3078.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": greps,
            "causal_order_audit_object": {
                "object": "Z12LightConeCausalOrderCompatibilityAudit",
                "source_reused": "P3078 frozen Dirichlet-to-wave promotion boundary",
                "candidate_sources": [candidate["id"] for candidate in CAUSAL_CANDIDATES],
                "required_gates": list(REQUIRED_GATES),
                "acceptance_predicate": "internal Z12 object plus exported metric signature, finite propagation cone, unit-normalized time, and wave EOM",
            },
            "cycle_distance_rows": distances,
            "distance_shell_rows": shells,
            "diffusion_support_rows": diffusion,
            "spectral_velocity_proxy_rows": velocities,
            "candidate_gate_rows": gates,
            "candidate_aggregate_certificate": aggs,
        },
        "finite_certificate": {
            "content_grep_lanes": len(greps),
            "content_grep_hits": sum(row["hit_count"] for row in greps),
            "p3078_accepted_intrinsic_momentum_symplectic_sources": p3078.get("finite_certificate", {}).get("accepted_intrinsic_momentum_symplectic_sources"),
            "cycle_distance_rows": len(distances),
            "distance_shell_rows": len(shells),
            "distance_shell_counts": {str(row["distance_shell"]): row["ordered_pair_count"] for row in shells},
            "diffusion_support_rows": len(diffusion),
            "diffusion_outside_unit_cone_violations_if_read_as_light": sum(1 for row in diffusion if row["violates_finite_unit_cone_if_read_as_light"]),
            "spectral_velocity_proxy_rows": len(velocities),
            "unit_light_speed_certified_rows": sum(1 for row in velocities if row["unit_light_speed_certified"]),
            "causal_candidates": len(CAUSAL_CANDIDATES),
            "required_gates": len(REQUIRED_GATES),
            "candidate_gate_rows": len(gates),
            "accepted_internal_causal_order_sources": len(accepted),
            "proof_obligations": len(proof_obligations),
            "satisfied_proof_obligations": sum(1 for row in proof_obligations if row["satisfied"]),
        },
        "proof_obligations": proof_obligations,
        "decision": {
            "bounded_result": "P3079 constructs a finite light-cone/causal-order compatibility audit for internal Z12 dispersion data.  The cycle graph gives a real internal spatial distance and adjacency iterates give a combinatorial cone, but diffusion has nonzero outside-unit-cone support and no candidate exports metric signature, unit-normalized time, and wave EOM together.  The imported Minkowski template passes only by external spacetime import.  Therefore no observed-light or spacetime-causal source is exported.",
            "negative_export_flags": {key: False for key in ["internal_metric_signature_exported", "unit_normalized_time_exported", "internal_light_cone_exported", "spacetime_wave_eom_exported", "observed_light_exported", "gauge_photon_sector_exported", "empirical_physics_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"cycle_distance_matrix_constructed": True, "diffusion_tail_obstruction_computed": True, "causal_candidate_gate_matrix_executed": True},
            "next_honest_step": "Pivot to a bounded internal-to-standard-physics interface inventory rather than another wave-promotion replay: construct a typed observable-interface obligation table for what would be needed to compare the smoothing/dispersion branch with standard theoretical physics (dimension map, continuum limit, Lorentz signature, gauge representation, conserved current, empirical observable), marking each obligation as sourced, formal/imported, or absent.  Do not claim observed light, photons, spacetime EOM, selector closure, L_total, bridge/role-transfer, or ToE unless one obligation is actually discharged by a new strict source theorem.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3079/S2029 light-cone/causal-order compatibility audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- P3078 accepted intrinsic momentum/symplectic sources: `{c['p3078_accepted_intrinsic_momentum_symplectic_sources']}`",
        f"- cycle distance rows: `{c['cycle_distance_rows']}`",
        f"- distance shell rows: `{c['distance_shell_rows']}`",
        f"- diffusion support rows: `{c['diffusion_support_rows']}`",
        f"- diffusion outside-unit-cone violations if read as light: `{c['diffusion_outside_unit_cone_violations_if_read_as_light']}`",
        f"- spectral velocity proxy rows: `{c['spectral_velocity_proxy_rows']}`",
        f"- unit light-speed certified rows: `{c['unit_light_speed_certified_rows']}`",
        f"- causal candidates: `{c['causal_candidates']}`",
        f"- candidate gate rows: `{c['candidate_gate_rows']}`",
        f"- accepted internal causal-order sources: `{c['accepted_internal_causal_order_sources']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3079/S2029 light-cone/causal-order compatibility audit", "## P3079/S2029 light-cone/causal-order compatibility audit\n\n`P3079/S2029` constructs `Z12LightConeCausalOrderCompatibilityAudit` after the `P3078` phase-space no-go.  It computes a `12 x 12` cycle-distance matrix, `7` distance-shell rows, `21` diffusion support rows over fractional steps, `12` spectral velocity proxy rows, and a `25`-row causal candidate gate matrix.  The internal cycle distance and adjacency cone are combinatorial/spatial; the heat kernel has nonzero outside-unit-cone tails if misread as light; the spectral velocities are dimensionless; and the Minkowski light cone works only as an imported ansatz.  No internal metric signature, unit-normalized time, spacetime wave EOM, observed light, gauge photons, empirical physics, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3079/S2029 causal-order audit remains non-spacetime", "## P3079/S2029 causal-order audit remains non-spacetime\n\n`P3079/S2029` shows that current Z12 dispersion data can support graph-distance and adjacency-reachability bookkeeping, but not a sourced Lorentzian metric, unit-normalized time coordinate, finite physical light cone, or spacetime wave EOM.  Diffusion remains smoothing-only, and any Minkowski light-cone interpretation is imported rather than derived.\n")
    append_once(AGENTS, "Current light-cone/causal-order compatibility guardrail (P3079/S2029, 2026-06-24)", "## Current light-cone/causal-order compatibility guardrail (P3079/S2029, 2026-06-24)\n\n- P3079 follows the P3078 recommendation and tests whether internal Z12 dispersion data already supply a metric signature, finite propagation cone, or unit-normalized causal order.\n- The finite audit has `144` cycle-distance rows, `21` diffusion support rows, `12` spectral velocity proxy rows, and `25` causal candidate gate rows; no candidate exports all internal metric-signature/time/EOM gates.\n- Do not promote cycle distance, adjacency reachability, diffusion tails, dimensionless velocity proxies, or imported Minkowski templates to observed light, gauge photons, spacetime EOM, empirical physics, `QW-2191` discharge, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is a bounded typed observable-interface obligation table for comparing the smoothing/dispersion branch with standard theoretical physics, without claiming any obligation is discharged unless a new strict source theorem provides it.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
