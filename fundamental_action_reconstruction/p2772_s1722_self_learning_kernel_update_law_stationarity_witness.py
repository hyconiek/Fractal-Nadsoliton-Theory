#!/usr/bin/env python3
"""P2772/S1722: self-learning kernel update-law stationarity witness.

P2771 tested one explicit geometric self-coupling candidate and recommended either
an ontologically sourced alternative C_geo or an explicit self-learning update
law.  This script takes the second route in bounded form: define a candidate
energy/loss from the same finite geometric self-coupling residual and use its
finite-difference gradient as a self-learning kernel-parameter update law.

The test is intentionally narrow.  If the current kernel tuple were already a
stationary learned fixed point for this candidate law, the gradient of the
normalized C_geo scalar-eigenclosure loss would vanish.  The finite calculation
finds nonzero gradients and nonzero one-step updates for both current kernel
tuples, so this candidate does not export a self-learning kernel theorem.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
P2771 = GEN / "p2771_s1721_geometric_self_coupling_operator_witness.json"
OUT = GEN / "p2772_s1722_self_learning_kernel_update_law_stationarity_witness.json"
MD = GEN / "p2772_s1722_self_learning_kernel_update_law_stationarity_witness.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 13
SHELLS = list(range((N - 1) // 2 + 1))
PARAMETER_SETS = {
    "K_legacy_ont": {
        "parameter_order": ["alpha_geo", "omega", "phi", "beta_tors"],
        "values": {"alpha_geo": 4.0 * math.log(2.0), "omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 0.01},
        "positive_parameters": {"alpha_geo", "beta_tors"},
    },
    "K_strict_gate": {
        "parameter_order": ["omega", "phi", "beta", "eta"],
        "values": {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 1.8},
        "positive_parameters": {"beta", "eta"},
    },
}
NEGATIVE_EXPORT_FLAGS = [
    "self_learning_kernel_update_theorem_exported",
    "learned_stationary_kernel_fixed_point_exported",
    "geometric_self_coupling_kernel_theorem_exported",
    "kernel_fully_expresses_nadsoliton_characteristics",
    "legacy_strict_bridge_closed",
    "role_bearing_ltotal_promoted",
    "selector_closure_exported",
    "toe_closure_exported",
]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def radial_distance(x: int, n: int = N) -> int:
    y = x % n
    return min(y, n - y)


def kernel_value(name: str, params: dict[str, float], d: int) -> float:
    if name == "K_legacy_ont":
        return params["alpha_geo"] * math.cos(params["omega"] * d + params["phi"]) / (1.0 + params["beta_tors"] * d)
    if name == "K_strict_gate":
        if params["beta"] <= 0.0 or params["eta"] <= 0.0:
            return float("nan")
        return math.cos(params["omega"] * d + params["phi"]) / (1.0 + params["beta"] * (d ** params["eta"]))
    raise ValueError(name)


def shell_samples(name: str, params: dict[str, float]) -> dict[int, float]:
    return {d: kernel_value(name, params, d) for d in SHELLS}


def cyclic_samples(name: str, params: dict[str, float]) -> list[float]:
    shells = shell_samples(name, params)
    return [shells[radial_distance(x)] for x in range(N)]


def c_geo(samples: list[float]) -> list[float]:
    return [sum(samples[x] * samples[(d - x) % N] for x in range(N)) for d in range(N)]


def radialize(values: list[float]) -> dict[int, float]:
    rows: dict[int, list[float]] = {d: [] for d in SHELLS}
    for x, value in enumerate(values):
        rows[radial_distance(x)].append(value)
    return {d: sum(vals) / len(vals) for d, vals in rows.items()}


def least_squares_scale(source: dict[int, float], target: dict[int, float]) -> float:
    numerator = sum(source[d] * target[d] for d in SHELLS)
    denominator = sum(source[d] * source[d] for d in SHELLS)
    return numerator / denominator


def normalized_c_geo_loss(name: str, params: dict[str, float]) -> float:
    if any(not math.isfinite(v) for v in params.values()):
        return float("inf")
    samples = shell_samples(name, params)
    if any(not math.isfinite(v) for v in samples.values()):
        return float("inf")
    coupled = radialize(c_geo(cyclic_samples(name, params)))
    gamma = least_squares_scale(samples, coupled)
    residual_sq = sum((coupled[d] - gamma * samples[d]) ** 2 for d in SHELLS)
    target_sq = sum(coupled[d] ** 2 for d in SHELLS)
    if target_sq == 0.0:
        return residual_sq
    return 0.5 * residual_sq / target_sq


def finite_difference_gradient(name: str, params: dict[str, float], order: list[str]) -> dict[str, float]:
    gradient: dict[str, float] = {}
    for key in order:
        base = params[key]
        step = 1e-5 * max(1.0, abs(base))
        p_plus = dict(params)
        p_minus = dict(params)
        p_plus[key] = base + step
        p_minus[key] = base - step
        if key in PARAMETER_SETS[name]["positive_parameters"] and p_minus[key] <= 0.0:
            p_minus[key] = base
            gradient[key] = (normalized_c_geo_loss(name, p_plus) - normalized_c_geo_loss(name, p_minus)) / (p_plus[key] - p_minus[key])
        else:
            gradient[key] = (normalized_c_geo_loss(name, p_plus) - normalized_c_geo_loss(name, p_minus)) / (2.0 * step)
    return gradient


def gradient_norm(gradient: dict[str, float]) -> float:
    return math.sqrt(sum(v * v for v in gradient.values()))


def one_step_update(params: dict[str, float], gradient: dict[str, float], learning_rate: float) -> dict[str, float]:
    return {key: params[key] - learning_rate * gradient[key] for key in params}


def kernel_update_witness(name: str) -> dict[str, Any]:
    spec = PARAMETER_SETS[name]
    params = dict(spec["values"])
    order = list(spec["parameter_order"])
    loss0 = normalized_c_geo_loss(name, params)
    grad = finite_difference_gradient(name, params, order)
    norm = gradient_norm(grad)
    learning_rate = 1e-3
    updated = one_step_update(params, grad, learning_rate)
    loss1 = normalized_c_geo_loss(name, updated)
    return {
        "kernel": name,
        "parameter_order": order,
        "initial_parameters": params,
        "loss_definition": "L_geo(theta)=0.5*||C_geo_N[K_theta]-gamma(theta)K_theta||^2/||C_geo_N[K_theta]||^2 on Z/13Z radial shells",
        "initial_loss": loss0,
        "finite_difference_gradient": grad,
        "gradient_norm": norm,
        "stationary_tolerance": 1e-9,
        "passes_stationary_tolerance": norm < 1e-9,
        "learning_rate": learning_rate,
        "one_step_updated_parameters": updated,
        "one_step_loss": loss1,
        "one_step_delta_loss": loss1 - loss0,
        "nonzero_update_parameters": [key for key in order if abs(updated[key] - params[key]) > 1e-12],
    }


def all_witnesses() -> dict[str, Any]:
    rows = [kernel_update_witness(name) for name in PARAMETER_SETS]
    return {
        "candidate_update_law": "theta_{t+1}=theta_t - lr * grad_theta L_geo(theta_t)",
        "finite_geometry": "Z/13Z radial shells inherited from P2771 C_geo_N",
        "row_count": len(rows),
        "rows": rows,
        "all_current_tuples_stationary": all(row["passes_stationary_tolerance"] for row in rows),
        "failed_stationary_kernels": [row["kernel"] for row in rows if not row["passes_stationary_tolerance"]],
        "min_gradient_norm": min(row["gradient_norm"] for row in rows),
        "max_gradient_norm": max(row["gradient_norm"] for row in rows),
    }


def acceptance_matrix(witness: dict[str, Any], p2771: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2771_c_geo_no_closure_boundary_present": p2771.get("status") == "P2771_FINITE_GEOMETRIC_SELF_COUPLING_OPERATOR_WITNESS_BOUNDED_NO_GO_NO_CLOSURE",
        "explicit_update_law_supplied": True,
        "all_current_kernel_tuples_stationary": witness["all_current_tuples_stationary"],
        "candidate_loss_ontologically_sourced": False,
        "candidate_update_coupled_to_physical_Ltotal": False,
        "candidate_exports_convergence_theorem": False,
    }
    return {
        "facts": facts,
        "accepted_as_self_learning_kernel_theorem": False,
        "accepted_as_learned_stationary_fixed_point": False,
        "accepted_as_geometric_self_coupling_theorem": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "The explicit gradient update law is computationally well-defined, but both current kernel tuples have nonzero gradients for the finite C_geo residual loss.  The loss is also not yet ontologically sourced or coupled to a physical L_total, so this is a bounded nonstationarity witness, not a self-learning kernel theorem.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    witness = payload["self_learning_update_witness"]
    lines = [
        "# P2772/S1722 self-learning kernel update-law stationarity witness",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Candidate update law",
        witness["candidate_update_law"],
        "",
        "## Result",
        f"- all_current_tuples_stationary={witness['all_current_tuples_stationary']}",
        f"- failed_stationary_kernels={witness['failed_stationary_kernels']}",
        f"- min_gradient_norm={witness['min_gradient_norm']}",
        f"- max_gradient_norm={witness['max_gradient_norm']}",
        "",
        "## Decision",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2771 = read_json(P2771)
    witness = all_witnesses()
    acceptance = acceptance_matrix(witness, p2771)
    payload = {
        "status": "P2772_SELF_LEARNING_KERNEL_UPDATE_LAW_STATIONARITY_WITNESS_BOUNDED_NO_GO_NO_CLOSURE",
        "input_hashes": {"P2771": sha(P2771)},
        "input_statuses": {"P2771": p2771.get("status")},
        "audited_object": "finite gradient self-learning update law for the P2771 C_geo residual loss",
        "self_learning_update_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not promote the current kernels to self-learning nadsoliton closure from this finite gradient update law.  The next honest move must supply either an ontologically sourced learning functional whose stationary equations are derived from the nadsoliton geometry, or a theorem that the current kernel tuple is a stable fixed point of such a sourced update.  Without that source/fixed-point theorem, preserve the P2697-P2772 no-full-expression/no-closure certificate and avoid L_total, bridge, role-transfer, selector, or ToE promotion.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2772/S1722 self-learning kernel update-law stationarity witness", "## P2772/S1722 self-learning kernel update-law stationarity witness\n\n`P2772/S1722` follows the other P2770/P2771 branch by introducing an explicit finite self-learning candidate update law `theta_{t+1}=theta_t-lr*grad L_geo(theta_t)`, where `L_geo` is the normalized P2771 `C_geo_N` scalar-eigenclosure residual on `Z/13Z`.  Finite-difference gradients at both current kernel tuples are nonzero, so neither tuple is a stationary learned fixed point of this candidate law.  The loss is also not yet ontologically sourced or coupled to a physical `L_total`; no self-learning kernel theorem, geometric self-coupling theorem, bridge closure, role transfer, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2772/S1722 update-law Ltotal guard", "## P2772/S1722 update-law Ltotal guard\n\n`P2772/S1722` adds no variational source term.  It defines one finite gradient update candidate from the P2771 self-coupling residual and finds nonzero gradients for both kernel tuples; therefore the candidate cannot promote the kernels to role-bearing `L_total` or self-learning nadsoliton dynamics.\n")
    append_once(AGENTS, "Current self-learning kernel update-law stationarity witness guardrail (P2772/S1722, 2026-06-15)", "## Current self-learning kernel update-law stationarity witness guardrail (P2772/S1722, 2026-06-15)\n\n- P2772 introduces the explicit P2771-recommended update candidate `theta_{t+1}=theta_t-lr*grad L_geo(theta_t)`, using the normalized finite `C_geo_N` scalar-eigenclosure residual as `L_geo`.\n- Both `K_legacy_ont` and `K_strict_gate` have nonzero finite-difference gradients, so neither current tuple is a stationary learned fixed point of this candidate; the loss is not ontologically sourced or coupled to physical `L_total`.\n- Do not promote this bounded update law to self-learning kernel closure, geometric self-coupling closure, full kernel expression, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure.  A next admissible move must supply an ontologically sourced learning functional or a sourced fixed-point theorem for the current tuple.\n")
    return payload


if __name__ == "__main__":
    main()
