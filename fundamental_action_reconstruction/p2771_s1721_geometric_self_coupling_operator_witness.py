#!/usr/bin/env python3
"""P2771/S1721: finite geometric self-coupling operator witness.

P2770 recommended introducing one explicit typed object before making any stronger
claim that the kernel fully expresses a geometrically self-coupled nadsoliton.
This script introduces a minimal selector-free candidate object:

    C_geo_N[K](d) = sum_x K(r_N(x)) K(r_N(d-x))

on the finite cyclic geometry Z/NZ with radial distance r_N(y)=min(y,N-y).
It then tests whether the current kernel samples are closed under this geometric
self-coupling up to a single scalar normalization, i.e. whether C_geo_N[K] is a
nonzero scalar multiple of K on all radial shells.  Failure is a bounded no-go
for this candidate only; it is not a no-go for all possible self-coupling laws.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
P2770 = GEN / "p2770_s1720_kernel_characteristic_expressivity_audit.json"
OUT = GEN / "p2771_s1721_geometric_self_coupling_operator_witness.json"
MD = GEN / "p2771_s1721_geometric_self_coupling_operator_witness.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 13
SHELLS = list(range((N - 1) // 2 + 1))
KERNEL_PARAMETERS = {
    "K_legacy_ont": {"alpha_geo": 4.0 * math.log(2.0), "omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 0.01},
    "K_strict_gate": {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 1.8},
}
NEGATIVE_EXPORT_FLAGS = [
    "geometric_self_coupling_kernel_theorem_exported",
    "kernel_self_coupling_fixed_point_exported",
    "kernel_fully_expresses_nadsoliton_characteristics",
    "legacy_strict_bridge_closed",
    "physical_role_transfer_started",
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


def kernel_value(name: str, d: int) -> float:
    p = KERNEL_PARAMETERS[name]
    if name == "K_legacy_ont":
        return p["alpha_geo"] * math.cos(p["omega"] * d + p["phi"]) / (1.0 + p["beta_tors"] * d)
    if name == "K_strict_gate":
        return math.cos(p["omega"] * d + p["phi"]) / (1.0 + p["beta"] * (d ** p["eta"]))
    raise ValueError(name)


def shell_samples(name: str) -> dict[int, float]:
    return {d: kernel_value(name, d) for d in SHELLS}


def cyclic_samples(name: str) -> list[float]:
    shells = shell_samples(name)
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
    if abs(denominator) < 1e-15:
        raise ValueError("zero source norm")
    return numerator / denominator


def relative_l2_residual(source: dict[int, float], target: dict[int, float], scale: float) -> float:
    residual = math.sqrt(sum((target[d] - scale * source[d]) ** 2 for d in SHELLS))
    norm = math.sqrt(sum(target[d] ** 2 for d in SHELLS))
    return residual / norm if norm else residual


def self_coupling_witness(name: str) -> dict[str, Any]:
    source = shell_samples(name)
    coupled = radialize(c_geo(cyclic_samples(name)))
    scale = least_squares_scale(source, coupled)
    residual_by_shell = {d: coupled[d] - scale * source[d] for d in SHELLS}
    max_abs_residual = max(abs(v) for v in residual_by_shell.values())
    rel_l2 = relative_l2_residual(source, coupled, scale)
    return {
        "kernel": name,
        "finite_geometry": "Z/13Z radial shells with r(x)=min(x,13-x)",
        "shell_order": SHELLS,
        "source_shell_values": source,
        "c_geo_shell_values": coupled,
        "least_squares_scalar": scale,
        "residual_by_shell": residual_by_shell,
        "max_abs_residual": max_abs_residual,
        "relative_l2_residual": rel_l2,
        "passes_scalar_eigenclosure_tolerance_1e_minus_9": rel_l2 < 1e-9,
    }


def all_witnesses() -> dict[str, Any]:
    rows = [self_coupling_witness(name) for name in KERNEL_PARAMETERS]
    return {
        "operator": "C_geo_N[K](d)=sum_x K(r_N(x))*K(r_N(d-x)) on N=13 cyclic radial geometry",
        "normalization_test": "C_geo_N[K] must equal gamma*K on all radial shells for one scalar gamma",
        "row_count": len(rows),
        "rows": rows,
        "all_kernels_pass_scalar_eigenclosure": all(row["passes_scalar_eigenclosure_tolerance_1e_minus_9"] for row in rows),
        "failed_kernels": [row["kernel"] for row in rows if not row["passes_scalar_eigenclosure_tolerance_1e_minus_9"]],
        "max_relative_l2_residual": max(row["relative_l2_residual"] for row in rows),
    }


def acceptance_matrix(witness: dict[str, Any], p2770: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2770_no_full_expression_boundary_present": p2770.get("status") == "P2770_KERNEL_CHARACTERISTIC_EXPRESSIVITY_AUDIT_NO_FULL_EXPRESSION_NO_CLOSURE",
        "explicit_c_geo_candidate_supplied": True,
        "all_kernel_samples_scalar_eigenclosed": witness["all_kernels_pass_scalar_eigenclosure"],
        "residuals_zero_to_tolerance": witness["max_relative_l2_residual"] < 1e-9,
        "candidate_is_unique_or_ontologically_sourced": False,
        "candidate_supplies_learning_update_law": False,
    }
    return {
        "facts": facts,
        "accepted_as_geometric_self_coupling_theorem": False,
        "accepted_as_kernel_full_expression_theorem": False,
        "accepted_as_learning_update_law": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "The explicit finite C_geo candidate is not a scalar eigenclosure for either current kernel, and it is not yet uniquely or ontologically sourced.  It therefore gives a bounded no-go for this candidate self-coupling law, not a full geometric self-coupling theorem.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    witness = payload["geometric_self_coupling_witness"]
    lines = [
        "# P2771/S1721 finite geometric self-coupling operator witness",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Operator",
        witness["operator"],
        "",
        "## Result",
        f"- all_kernels_pass_scalar_eigenclosure={witness['all_kernels_pass_scalar_eigenclosure']}",
        f"- failed_kernels={witness['failed_kernels']}",
        f"- max_relative_l2_residual={witness['max_relative_l2_residual']}",
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
    p2770 = read_json(P2770)
    witness = all_witnesses()
    acceptance = acceptance_matrix(witness, p2770)
    payload = {
        "status": "P2771_FINITE_GEOMETRIC_SELF_COUPLING_OPERATOR_WITNESS_BOUNDED_NO_GO_NO_CLOSURE",
        "input_hashes": {"P2770": sha(P2770)},
        "input_statuses": {"P2770": p2770.get("status")},
        "audited_object": "selector-free finite cyclic radial geometric self-coupling operator C_geo_N[K]",
        "geometric_self_coupling_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not promote the current kernels to geometrically self-coupled nadsoliton closure from this C_geo candidate.  The next honest move is either to supply a strictly sourced geometric self-coupling operator with a different, justified geometry/normalization and rerun the scalar-eigenclosure residual table, or to pivot to the other P2770 branch: an explicit self-learning kernel-parameter update law with provenance and a bounded convergence/consistency witness.  Without one of those, preserve the P2697-P2771 no-full-expression/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2771/S1721 finite geometric self-coupling operator witness", "## P2771/S1721 finite geometric self-coupling operator witness\n\n`P2771/S1721` follows the P2770 recommendation by introducing one explicit typed candidate, the selector-free finite cyclic radial operator `C_geo_N[K](d)=sum_x K(r_N(x))*K(r_N(d-x))` on `Z/13Z`.  The bounded scalar-eigenclosure test asks whether `C_geo_N[K] = gamma K` on all radial shells.  Both `K_legacy_ont` and `K_strict_gate` fail with nonzero residuals, so this candidate does not export a geometric self-coupling theorem, kernel full-expression theorem, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2771/S1721 C_geo Ltotal guard", "## P2771/S1721 C_geo Ltotal guard\n\n`P2771/S1721` adds no variational source term.  It tests one explicit finite geometric self-coupling candidate and finds nonzero scalar-eigenclosure residuals for both current kernels; therefore this candidate cannot promote the kernel formulas to role-bearing `L_total` or a full self-coupled nadsoliton dynamics.\n")
    append_once(AGENTS, "Current finite geometric self-coupling operator witness guardrail (P2771/S1721, 2026-06-15)", "## Current finite geometric self-coupling operator witness guardrail (P2771/S1721, 2026-06-15)\n\n- P2771 introduces the explicit P2770-requested candidate object `C_geo_N[K](d)=sum_x K(r_N(x))*K(r_N(d-x))` on selector-free `Z/13Z` radial geometry and tests scalar eigenclosure `C_geo_N[K]=gamma K`.\n- Both `K_legacy_ont` and `K_strict_gate` fail the finite scalar-eigenclosure residual test; the candidate is also not uniquely/ontologically sourced and supplies no learning update law.\n- Do not promote this bounded candidate to geometric self-coupling closure, full kernel expression, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure.  A next admissible move must either supply a strictly sourced alternative `C_geo` with justified geometry/normalization and rerun residuals, or pivot to an explicit self-learning kernel-parameter update law with bounded convergence/provenance witness.\n")
    return payload


if __name__ == "__main__":
    main()
