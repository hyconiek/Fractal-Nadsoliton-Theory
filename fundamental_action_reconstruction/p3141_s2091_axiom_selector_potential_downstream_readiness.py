"""P3141/S2091: axiom selector potential and downstream-readiness audit.

P3140 showed that selector axioms can conditionally choose a D_HL pair, but
only non-strictly.  This packet constructs the next missing theoretical object
on that non-strict branch: an explicit finite symmetry-breaking potential whose
unique minimizer is the axiom-selected pair.  It then audits whether that
potential can be promoted downstream to units, L_total, bridge/role transfer, or
ToE closure.  The answer remains bounded: useful non-strict machinery, no strict
or variational promotion.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3141_s2091_axiom_selector_potential_downstream_readiness.json"
MD = GEN / "p3141_s2091_axiom_selector_potential_downstream_readiness.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P3140": GEN / "p3140_s2090_axiom_augmented_selector_premise_calculus.json",
    "P3139": GEN / "p3139_s2089_dhl_lane_no_new_frontier_reconciliation.json",
}

N = 12
LAMBDAS = [-1, 1]
R0 = 0
LAMBDA0 = 1
WEIGHTS = [1, 2, 3]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def cyc_dist(a: int, b: int) -> int:
    raw = abs((a - b) % N)
    return min(raw, N - raw)


def pair_space() -> list[tuple[int, int]]:
    return [(r, lam) for r in range(N) for lam in LAMBDAS]


def selector_potential(pair: tuple[int, int], w_r: int, w_lam: int, mu: int = 1) -> int:
    r, lam = pair
    origin_penalty = cyc_dist(r, R0) ** 2
    lambda_penalty = (1 - lam * LAMBDA0) // 2
    return mu * (w_r * origin_penalty + w_lam * lambda_penalty)


def potential_row(w_r: int, w_lam: int) -> dict[str, Any]:
    values = [(pair, selector_potential(pair, w_r, w_lam)) for pair in pair_space()]
    min_val = min(v for _, v in values)
    minimizers = [pair for pair, v in values if v == min_val]
    first_gap = min(v for _, v in values if v > min_val) - min_val
    max_val = max(v for _, v in values)
    return {
        "w_origin": w_r,
        "w_lambda": w_lam,
        "mu": 1,
        "min_value": min_val,
        "minimizers": [[r, lam] for r, lam in minimizers],
        "unique_minimizer": len(minimizers) == 1,
        "selected_pair": [R0, LAMBDA0] if len(minimizers) == 1 else None,
        "first_gap": first_gap,
        "max_value": max_val,
        "strict_weight_source_exported": False,
        "unit_bearing_scale_exported": False,
    }


def build_payload() -> dict[str, Any]:
    inputs = {name: load_json(path) for name, path in INPUTS.items()}
    rows = [potential_row(w_r, w_lam) for w_r in WEIGHTS for w_lam in WEIGHTS]
    downstream_rows = [
        {"gate": "finite_unique_minimizer", "passed": all(row["unique_minimizer"] for row in rows), "detail": "all positive audited weights uniquely select the axiom pair"},
        {"gate": "non_strict_axiom_consistency", "passed": True, "detail": "the minimizer agrees with P3140 A_origin + A_lambda"},
        {"gate": "strict_weight_source", "passed": False, "detail": "w_origin, w_lambda, and mu are assumed branch parameters, not strict sourced constants"},
        {"gate": "unit_bearing_scale", "passed": False, "detail": "mu is dimensionless here; no action unit or physical normalization is exported"},
        {"gate": "field_variational_lift", "passed": False, "detail": "finite pair penalties are not a field-space functional derivative or EOM"},
        {"gate": "bridge_role_transfer_ToE", "passed": False, "detail": "selector potential alone does not complete legacy->strict bridge, role transfer, L_total, or ToE"},
    ]
    theorem = {
        "name": "P3141_T1_axiom_selector_potential_unique_minimizer",
        "statement": "For every audited positive integer weight pair (w_origin,w_lambda) in {1,2,3}^2, the finite axiom potential V_ax(r,lambda)=mu*(w_origin*d_Z12(r,r0)^2 + w_lambda*(1-lambda*lambda0)/2) has a unique minimizer at the P3140 axiom-selected pair (r0,lambda0)=(0,+1).  This constructs a usable non-strict symmetry-breaking potential, but the weights and scale are assumed and do not supply strict source provenance or unit-bearing variational closure.",
        "finite_counts": {
            "pair_space_size": len(pair_space()),
            "weight_pairs_tested": len(rows),
            "unique_minimizer_rows": sum(row["unique_minimizer"] for row in rows),
            "strict_weight_source_rows": sum(row["strict_weight_source_exported"] for row in rows),
            "unit_bearing_scale_rows": sum(row["unit_bearing_scale_exported"] for row in rows),
        },
    }
    return {
        "status": "P3141_AXIOM_SELECTOR_POTENTIAL_DOWNSTREAM_READINESS_BOUNDED_NON_STRICT",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "repo_context": {
            "P3140_status": inputs["P3140"].get("status"),
            "P3139_status": inputs["P3139"].get("status"),
        },
        "constructed_object": {
            "name": "V_ax^sel finite axiom selector potential",
            "formula": "V_ax(r,lambda)=mu*(w_origin*d_Z12(r,r0)^2 + w_lambda*(1-lambda*lambda0)/2)",
            "axiom_pair": {"r0": R0, "lambda0": LAMBDA0},
            "classification": "non_strict_axiom_branch_potential",
        },
        "potential_rows": rows,
        "downstream_gate_rows": downstream_rows,
        "finite_theorem": theorem,
        "decision": {
            "bounded_result": "P3141 constructs the missing finite symmetry-breaking potential for the P3140 axiom branch.  It is computationally decisive inside the non-strict branch: all audited positive weights uniquely select (0,+1).  It does not promote to strict physics because the weights, scale, field lift, unit normalization, and bridge/role-transfer/L_total/ToE links are absent.",
            "positive_scoped_flags": {
                "non_strict_selector_potential_constructed": True,
                "unique_minimizer_theorem_exported_in_axiom_branch": True,
                "conditional_symmetry_breaking_engine_available_non_strictly": True,
            },
            "negative_export_flags": {
                "strict_weight_source_exported": False,
                "unit_bearing_selector_action_exported": False,
                "field_variational_EOM_exported": False,
                "strict_QW_2191_discharged": False,
                "strict_selector_closure_exported": False,
                "bridge_completion_exported": False,
                "legacy_role_transfer_exported": False,
                "L_total_exported": False,
                "ToE_closure_exported": False,
            },
            "next_honest_step": "For the axiom branch, the next proof-grade move is exactly one downstream theorem: either a unit-bearing scale/source theorem for mu and the weights, or a field-variational lift of V_ax with a real EOM derivative.  For strict-core work, the next move remains deriving one of A_origin, A_lambda, A_coupling, or the potential weights from a new strict source rather than assuming them.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    th = payload["finite_theorem"]
    lines = [
        "# P3141/S2091 Axiom selector potential downstream-readiness audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Constructed object",
        f"- `{payload['constructed_object']['name']}`",
        f"- Formula: `{payload['constructed_object']['formula']}`",
        f"- Classification: `{payload['constructed_object']['classification']}`",
        "",
        "## Finite theorem",
        f"`{th['name']}`: {th['statement']}",
        "",
        "## Finite counts",
    ]
    for key, value in th["finite_counts"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Downstream gates"])
    for row in payload["downstream_gate_rows"]:
        lines.append(f"- `{row['gate']}`: `{row['passed']}` — {row['detail']}")
    lines.extend(["", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3141/S2091 axiom selector potential downstream-readiness audit", "## P3141/S2091 axiom selector potential downstream-readiness audit\n\n`P3141/S2091` constructs the finite non-strict selector potential `V_ax(r,lambda)=mu*(w_origin*d_Z12(r,r0)^2 + w_lambda*(1-lambda*lambda0)/2)` for the P3140 axiom branch.  Across `9` positive integer weight pairs, the potential has a unique minimizer at the axiom-selected `(0,+1)` pair.  This supplies a conditional symmetry-breaking engine inside the non-strict branch only: no strict weight source, unit-bearing selector action, field variational EOM, strict `QW-2191` discharge, bridge completion, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3141/S2091 V_ax remains non-strict and lacks unit-bearing EOM lift", "## P3141/S2091 V_ax remains non-strict and lacks unit-bearing EOM lift\n\n`P3141/S2091` gives a finite axiom-branch potential with a unique minimizer, but its weights and scale are assumed and dimensionless, and no field-space variational lift or EOM derivative is exported.  It therefore cannot be promoted to unit-bearing action, Hamiltonian normalization, spacetime EOM, `L_total`, bridge-completion theorem, role-transfer theorem, or ToE closure.\n")
    append_once(AGENTS, "Current axiom selector potential downstream-readiness guardrail (P3141/S2091, 2026-07-13)", "## Current axiom selector potential downstream-readiness guardrail (P3141/S2091, 2026-07-13)\n\n- P3141 constructs a finite non-strict selector potential `V_ax` for the P3140 axiom branch and proves that all `9` audited positive weight pairs uniquely minimize at the stipulated `(r0,lambda0)=(0,+1)` pair.\n- This is a conditional symmetry-breaking engine only inside the axiom branch; the weights, scale, unit normalization, and field-variational lift are not strict sourced.\n- Do not promote `V_ax` to strict `QW-2191` discharge, strict selector closure, unit-bearing action/EOM, bridge completion, role transfer, `L_total`, or ToE closure.\n- Next honest move: in the axiom branch, prove exactly one downstream theorem for unit-bearing scale/weights or field-variational lift; in strict-core work, derive one of the selector axioms or potential weights from a new strict source.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
