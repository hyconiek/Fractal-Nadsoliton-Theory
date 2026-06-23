#!/usr/bin/env python3
"""P3040/S1990: retarded path-anisotropy source obstruction.

This is the next single-premise audit after P3039.  P3038 used a candidate
c-retarded split
    Delta_i = cos(omega*(1+rho*M_i)) - cos(omega*(1-rho*M_i))
where M_i is a memory/viscosity trace.  P3040 asks whether current strict data
source the missing path-anisotropy theorem: a nonpremise law fixing the retarded
parallel/perpendicular split and a positive rho, rather than inserting them as
candidate slots.

The result is a bounded obstruction.  The split is finite and nonzero, and its
sign torsor is explicit, but rho/sign/path choices remain candidate gauges,
target-tuned receivers, or chart-order receivers rather than strict sourced path
geometry.
"""
from __future__ import annotations

import hashlib, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3023_s1973_kernel_dissipation_time_order_candidate_obstruction import N, k_strict
from p3038_s1988_viscous_retarded_chiral_selector_candidate_obstruction import OUT as P3038, RETARD_OMEGA, RETARD_RHO, memory_viscosity_trace, chiral_projector
from p3039_s1989_chiral_projector_sign_localizer_obstruction import OUT as P3039

OUT = GEN / "p3040_s1990_retarded_path_anisotropy_source_obstruction.json"
MD = GEN / "p3040_s1990_retarded_path_anisotropy_source_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
TOL = 1e-12
RHO_GRID = [0.0, 0.05, 0.10, 0.25, 0.50, 1.00]
PATH_SIGNS = [+1, -1]


def labels() -> list[int]:
    return list(range(N))


def kernel_vector() -> list[float]:
    return [k_strict(i + 1) for i in labels()]


def cyclic_gradient(values: list[float]) -> list[float]:
    return [values[(i + 1) % N] - values[(i - 1) % N] for i in labels()]


def cyclic_curvature(values: list[float]) -> list[float]:
    return [values[(i + 1) % N] - 2.0 * values[i] + values[(i - 1) % N] for i in labels()]


def retarded_split(trace: list[float], rho: float, path_sign: int = +1, omega: float = RETARD_OMEGA) -> list[float]:
    signed_rho = path_sign * rho
    return [math.cos(omega * (1.0 + signed_rho * value)) - math.cos(omega * (1.0 - signed_rho * value)) for value in trace]


def branch_readout_score(delta: list[float]) -> float:
    # Keep the P3038 readout slot fixed only as a receiver; P3039 blocks its source status.
    return sum(d * c for d, c in zip(delta, chiral_projector()))


def candidate_grid(trace: list[float]) -> list[dict[str, Any]]:
    rows = []
    for rho in RHO_GRID:
        for sign in PATH_SIGNS:
            delta = retarded_split(trace, rho, sign)
            rows.append({
                "rho": rho,
                "path_sign": sign,
                "delta_l1": round(sum(abs(x) for x in delta), 15),
                "readout_score": round(branch_readout_score(delta), 15),
                "finite_nonzero": any(abs(x) > TOL for x in delta),
            })
    return rows


def best_grid_rows(grid: list[dict[str, Any]]) -> list[dict[str, Any]]:
    max_abs = max(abs(row["readout_score"]) for row in grid)
    return [row for row in grid if abs(abs(row["readout_score"]) - max_abs) <= TOL]


def build_matrix() -> dict[str, Any]:
    k = kernel_vector()
    trace = memory_viscosity_trace(k)
    base_delta = retarded_split(trace, RETARD_RHO, +1)
    flipped_delta = retarded_split(trace, RETARD_RHO, -1)
    sign_flips_delta = max(abs(a + b) for a, b in zip(base_delta, flipped_delta)) < 1e-12
    grid = candidate_grid(trace)
    best = best_grid_rows(grid)
    grad = cyclic_gradient(k)
    curv = cyclic_curvature(trace)
    receivers = [
        {
            "receiver": "bare_retarded_path_split_torsor",
            "formula": "Delta_{rho,sigma}(i)=cos(omega*(1+sigma*rho*M_i))-cos(omega*(1-sigma*rho*M_i))",
            "finite_nonzero": any(abs(x) > TOL for x in base_delta),
            "path_sign_torsor_flips_delta": sign_flips_delta,
            "accepted_as_path_anisotropy_source": False,
            "failure": "the finite split is real, but sigma and rho are inserted path-gauge slots rather than sourced geometry",
        },
        {
            "receiver": "rho_grid_max_branch_readout_receiver",
            "grid_rows": grid,
            "best_rows": best,
            "finite_nonzero": any(row["finite_nonzero"] for row in grid),
            "accepted_as_path_anisotropy_source": False,
            "failure": "maximizing the P3038 branch readout tunes rho/sign to the target receiver; it is not a target-independent path theorem",
        },
        {
            "receiver": "kernel_cyclic_gradient_path_receiver",
            "gradient_l1": round(sum(abs(x) for x in grad), 15),
            "finite_nonzero": any(abs(x) > TOL for x in grad),
            "accepted_as_path_anisotropy_source": False,
            "failure": "cyclic gradients are computable but require a chosen chart order and do not name retarded path geometry",
        },
        {
            "receiver": "memory_curvature_path_receiver",
            "curvature_l1": round(sum(abs(x) for x in curv), 15),
            "finite_nonzero": any(abs(x) > TOL for x in curv),
            "accepted_as_path_anisotropy_source": False,
            "failure": "memory curvature is a receiver-side shape diagnostic, not an exported law fixing parallel/perpendicular paths or rho",
        },
    ]
    obligations = [
        {"obligation": "exact_p3038_path_anisotropy_premise_targeted", "satisfied": True, "detail": "only the retarded path-anisotropy source premise is audited"},
        {"obligation": "finite_retarded_split_torsor_constructed", "satisfied": True, "detail": "rho grid and path-sign torsor are explicit"},
        {"obligation": "nonzero_candidate_split_verified", "satisfied": any(abs(x) > TOL for x in base_delta), "detail": "P3038 rho gives a nonzero split"},
        {"obligation": "path_sign_flip_verified", "satisfied": sign_flips_delta, "detail": "sigma -> -sigma sends Delta to -Delta"},
        {"obligation": "nonpremise_positive_rho_source", "satisfied": False, "detail": "rho remains a parameter/gauge; grid choice is target tuning"},
        {"obligation": "strict_parallel_perpendicular_path_geometry", "satisfied": False, "detail": "no strict source names the two retarded path sectors"},
        {"obligation": "chart_independent_path_localizer", "satisfied": False, "detail": "gradient/curvature receivers depend on the finite cyclic chart/order"},
        {"obligation": "coupling_back_to_p3038_as_source", "satisfied": False, "detail": "without sourced rho and path sectors, P3038 remains a candidate operator"},
    ]
    return {
        "object": "RetardedPathAnisotropySource_ObstructionMatrix",
        "target_formula": "Delta_i = cos(omega*(1+rho*M_i))-cos(omega*(1-rho*M_i)) from P3038",
        "parameters": {"retard_omega": RETARD_OMEGA, "p3038_rho": RETARD_RHO, "rho_grid": RHO_GRID, "path_signs": PATH_SIGNS},
        "torsor": {"path_sign_flip_verified": sign_flips_delta, "rho_grid_rows": len(grid), "best_grid_rows": len(best)},
        "receiver_rows": receivers,
        "proof_obligations": obligations,
        "finite_certificate": {
            "receiver_rows": len(receivers),
            "finite_nonzero_rows": sum(1 for row in receivers if row["finite_nonzero"]),
            "accepted_path_anisotropy_source_rows": sum(1 for row in receivers if row["accepted_as_path_anisotropy_source"]),
            "rho_grid_rows": len(grid),
            "best_grid_rows": len(best),
            "path_sign_flip_verified": sign_flips_delta,
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "sourced_retardation_path_anisotropy_exported": all(row["satisfied"] for row in obligations) and all(row["accepted_as_path_anisotropy_source"] for row in receivers),
        },
    }


def build_payload() -> dict[str, Any]:
    read_json(P3038)
    read_json(P3039)
    matrix = build_matrix()
    return {
        "status": "P3040_RETARDED_PATH_ANISOTROPY_SOURCE_OBSTRUCTION_NO_SOURCE_EXPORT",
        "input_hashes": {
            "P3038": hashlib.sha256(P3038.read_bytes()).hexdigest() if P3038.exists() else None,
            "P3039": hashlib.sha256(P3039.read_bytes()).hexdigest() if P3039.exists() else None,
        },
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "The P3038 retarded split is finite and nonzero, and the path-sign torsor is explicit: reversing the sign sends Delta to -Delta.  But current receivers only insert rho/sign, tune them to the branch readout, or use chart-order gradient/curvature diagnostics.  No strict path geometry exports a positive rho, parallel/perpendicular retarded sectors, or a chart-independent path localizer.",
            "negative_export_flags": {k: False for k in ["sourced_retardation_path_anisotropy_exported", "p3038_promoted_to_selector_source", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay rho-grid tuning, cyclic gradients, memory curvature, or inserted parallel/perpendicular labels as a sourced retarded geometry.  The remaining direct P3038 premise is a physical unit/readout coupling theorem for the integrated density; otherwise preserve P3038-P3040 as a finite branch-separating but unsourced selector-candidate certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3040/S1990 retarded path-anisotropy source obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- receiver rows: `{c['receiver_rows']}`",
        f"- finite nonzero rows: `{c['finite_nonzero_rows']}`",
        f"- accepted path-anisotropy source rows: `{c['accepted_path_anisotropy_source_rows']}`",
        f"- rho grid rows: `{c['rho_grid_rows']}`",
        f"- best grid rows: `{c['best_grid_rows']}`",
        f"- path sign flip verified: `{c['path_sign_flip_verified']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- sourced retardation path anisotropy exported: `{c['sourced_retardation_path_anisotropy_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3040/S1990 retarded path-anisotropy source obstruction", "## P3040/S1990 retarded path-anisotropy source obstruction\n\n`P3040/S1990` attacks exactly one remaining P3038 premise: a sourced retarded path-anisotropy theorem for `Delta_i = cos(omega*(1+rho*M_i))-cos(omega*(1-rho*M_i))`.  The finite split is nonzero and the path-sign torsor is explicit (`sigma -> -sigma` sends `Delta` to `-Delta`), but the tested rows either insert `rho/sign`, tune a rho grid to the P3038 readout, or use chart-order gradient/curvature receivers.  No strict source exports positive `rho`, parallel/perpendicular retarded path sectors, or a chart-independent path localizer.  Thus P3040 does not promote P3038 to selector closure and does not discharge `QW-2191`, observed physics, `L_total`, bridge/role transfer, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3040/S1990 retarded path anisotropy `L_total` guard", "## P3040/S1990 retarded path anisotropy `L_total` guard\n\n`P3040/S1990` adds no physical `L_total` term.  It verifies a nonzero finite retarded split and path-sign torsor, but no strict path geometry sources `rho`, parallel/perpendicular sectors, or chart-independent path localization; the P3038 density remains a candidate receiver.\n")
    append_once(AGENTS, "Current retarded path-anisotropy source guardrail (P3040/S1990, 2026-06-23)", "## Current retarded path-anisotropy source guardrail (P3040/S1990, 2026-06-23)\n\n- P3040 attacks exactly one P3038 missing source premise: a sourced retarded path-anisotropy theorem for the `Delta_i` split.\n- The split is finite and nonzero, and the path-sign torsor flips `Delta`, but `rho`/sign insertion, rho-grid readout tuning, cyclic gradients, and memory curvature are receiver-side diagnostics rather than strict path geometry.\n- Do not promote rho-grid tuning, inserted parallel/perpendicular labels, cyclic-gradient path receivers, memory-curvature receivers, or finite retarded splitting alone to `QW-2191` discharge, selector closure, observed physics, unit-bearing action/EOM, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move may attack the remaining direct P3038 premise: physical unit/readout coupling theorem for the integrated density, but only with a genuinely new strict source law.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
