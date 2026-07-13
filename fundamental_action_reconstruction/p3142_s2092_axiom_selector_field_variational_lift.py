"""P3142/S2092: non-strict axiom selector field-variational lift.

P3141 built a finite axiom-branch selector potential but left the next
axiomatic downstream object open: a field/local-chart variational lift with an
actual derivative and Hessian.  This packet constructs that object honestly as a
non-strict axiom-branch toy functional, computes its stationary equations on a
local chart, and audits why it still cannot be promoted to strict-source,
unit-bearing L_total, bridge, role-transfer, or ToE closure.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3142_s2092_axiom_selector_field_variational_lift.json"
MD = GEN / "p3142_s2092_axiom_selector_field_variational_lift.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P3141": GEN / "p3141_s2091_axiom_selector_potential_downstream_readiness.json",
    "P3140": GEN / "p3140_s2090_axiom_augmented_selector_premise_calculus.json",
}

THETA0 = 0
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


def chart_symbolics() -> dict[str, Any]:
    theta, s, mu, w_theta, w_s, kappa = sp.symbols("theta s mu w_theta w_s kappa", positive=True, real=True)
    # Local chart lift of P3141: theta is a local coordinate at selected r0;
    # s is a continuous polarity field constrained by a double-well around ±1.
    functional = mu * (w_theta * theta**2 + w_s * (s - LAMBDA0) ** 2 / 2 + kappa * (s**2 - 1) ** 2)
    gradient = [sp.diff(functional, theta), sp.diff(functional, s)]
    hessian = sp.hessian(functional, (theta, s))
    at_selected = {theta: THETA0, s: LAMBDA0}
    gradient_at_selected = [sp.simplify(g.subs(at_selected)) for g in gradient]
    hessian_at_selected = sp.simplify(hessian.subs(at_selected))
    eigenvals = {str(k): str(v) for k, v in hessian_at_selected.eigenvals().items()}
    return {
        "variables": ["theta", "s"],
        "formula": str(functional),
        "gradient": [str(g) for g in gradient],
        "gradient_at_selected": [str(g) for g in gradient_at_selected],
        "hessian_at_selected": str(hessian_at_selected),
        "hessian_eigenvalues_at_selected": eigenvals,
        "positive_definite_for_positive_parameters": True,
    }


def parameter_row(w_theta: int, w_s: int, kappa: int) -> dict[str, Any]:
    theta, s = sp.symbols("theta s", real=True)
    mu = 1
    functional = mu * (w_theta * theta**2 + w_s * (s - 1) ** 2 / 2 + kappa * (s**2 - 1) ** 2)
    grad = [sp.diff(functional, theta), sp.diff(functional, s)]
    hessian = sp.hessian(functional, (theta, s)).subs({theta: 0, s: 1})
    eigs = [int(ev) for ev in hessian.eigenvals().keys()]
    return {
        "w_theta": w_theta,
        "w_s": w_s,
        "kappa": kappa,
        "gradient_at_selected": [str(g.subs({theta: 0, s: 1})) for g in grad],
        "hessian_at_selected": [[int(hessian[i, j]) for j in range(2)] for i in range(2)],
        "positive_eigenvalues": eigs,
        "strict_weight_source_exported": False,
        "unit_bearing_measure_exported": False,
        "global_Z12_chart_exported": False,
    }


def build_payload() -> dict[str, Any]:
    inputs = {name: load_json(path) for name, path in INPUTS.items()}
    rows = [parameter_row(a, b, c) for a in WEIGHTS for b in WEIGHTS for c in WEIGHTS]
    gates = [
        {"gate": "local_variational_derivative", "passed": True, "detail": "explicit dV/dtheta and dV/ds are computed"},
        {"gate": "selected_point_stationary", "passed": all(row["gradient_at_selected"] == ["0", "0"] for row in rows), "detail": "all audited positive parameter triples are stationary at the axiom-selected branch"},
        {"gate": "local_positive_hessian", "passed": all(min(row["positive_eigenvalues"]) > 0 for row in rows), "detail": "the selected point is a strict local minimum in the local chart"},
        {"gate": "strict_source_weights", "passed": False, "detail": "mu, w_theta, w_s, and kappa are assumed axiom-branch parameters"},
        {"gate": "global_Z12_field_chart", "passed": False, "detail": "the construction uses a chosen local chart at r0 and does not derive the global Z12 origin"},
        {"gate": "unit_bearing_L_total_ToE", "passed": False, "detail": "no action measure, physical units, spacetime coupling, bridge role-transfer, or ToE closure is exported"},
    ]
    theorem = {
        "name": "P3142_T1_axiom_selector_local_variational_lift",
        "statement": "The local-chart functional V_lift(theta,s)=mu*(w_theta*theta^2 + w_s*(s-1)^2/2 + kappa*(s^2-1)^2) has zero gradient and positive Hessian at the P3140/P3141 axiom-selected branch (theta,s)=(0,+1) for every audited positive integer triple in {1,2,3}^3.  This is a non-strict field-variational lift of the axiom selector potential, not a strict source theorem or unit-bearing action.",
        "finite_counts": {
            "parameter_triples_tested": len(rows),
            "stationary_rows": sum(row["gradient_at_selected"] == ["0", "0"] for row in rows),
            "positive_hessian_rows": sum(min(row["positive_eigenvalues"]) > 0 for row in rows),
            "strict_weight_source_rows": sum(row["strict_weight_source_exported"] for row in rows),
            "unit_bearing_measure_rows": sum(row["unit_bearing_measure_exported"] for row in rows),
        },
    }
    return {
        "status": "P3142_AXIOM_SELECTOR_FIELD_VARIATIONAL_LIFT_BOUNDED_NON_STRICT",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "repo_context": {"P3141_status": inputs["P3141"].get("status"), "P3140_status": inputs["P3140"].get("status")},
        "constructed_object": {
            "name": "V_lift^ax local-chart axiom selector functional",
            "classification": "non_strict_axiom_branch_local_variational_lift",
            "depends_on_assumptions": ["A_origin", "A_lambda", "A_coupling", "positive assumed weights", "chosen local chart at r0"],
        },
        "symbolic_variation": chart_symbolics(),
        "parameter_rows": rows,
        "gate_rows": gates,
        "finite_theorem": theorem,
        "decision": {
            "bounded_result": "P3142 constructs the requested axiomatic downstream field-variational object and verifies its derivative/Hessian computationally.  It gives a real conditional symmetry-breaking engine in the non-strict branch, but the engine still imports the selector axioms, weights, local chart, scale, and unit measure.",
            "negative_export_flags": {
                "strict_QW_2191_discharged": False,
                "strict_selector_closure_exported": False,
                "global_Z12_origin_derived": False,
                "unit_bearing_action_exported": False,
                "spacetime_EOM_exported": False,
                "bridge_completion_exported": False,
                "legacy_role_transfer_exported": False,
                "L_total_exported": False,
                "ToE_closure_exported": False,
            },
            "theoretical_perspective": "A selector is symmetry-breaking data: it chooses one representative/orientation from an otherwise symmetric orbit.  In this axiom branch it can start a conditional engine by making one channel locally stable, but strict ToE probability does not become high until the selector is sourced non-premise and then coupled to units, variational dynamics, bridge completion, and role transfer.",
            "next_honest_step": "Do one narrow source audit for the weights/scale of V_lift: construct exactly one candidate strict source for mu, w_theta, w_s, or kappa and test whether it is import-free, unit-bearing, and compatible with the global Z12 quotient.  If no such source is supplied, preserve the P3140-P3142 non-strict axiom-branch boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    th = payload["finite_theorem"]
    lines = [
        "# P3142/S2092 Axiom selector field-variational lift audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Constructed object",
        f"- `{payload['constructed_object']['name']}`",
        f"- Classification: `{payload['constructed_object']['classification']}`",
        f"- Assumptions: `{', '.join(payload['constructed_object']['depends_on_assumptions'])}`",
        "",
        "## Symbolic variation",
        f"- Formula: `{payload['symbolic_variation']['formula']}`",
        f"- Gradient: `{payload['symbolic_variation']['gradient']}`",
        f"- Hessian at selected branch: `{payload['symbolic_variation']['hessian_at_selected']}`",
        "",
        "## Finite theorem",
        f"`{th['name']}`: {th['statement']}",
        "",
        "## Finite counts",
    ]
    for key, value in th["finite_counts"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Gates"])
    for row in payload["gate_rows"]:
        lines.append(f"- `{row['gate']}`: `{row['passed']}` — {row['detail']}")
    lines.extend(["", "## Decision", payload["decision"]["bounded_result"], "", "## Theoretical perspective", payload["decision"]["theoretical_perspective"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3142/S2092 axiom selector field-variational lift audit", "## P3142/S2092 axiom selector field-variational lift audit\n\n`P3142/S2092` constructs the non-strict local-chart lift `V_lift(theta,s)=mu*(w_theta*theta^2 + w_s*(s-1)^2/2 + kappa*(s^2-1)^2)` for the P3140/P3141 axiom selector branch.  Across `27` positive parameter triples, the selected branch `(theta,s)=(0,+1)` has zero gradient and positive Hessian.  This is a real conditional variational engine only in the axiom branch: no strict source for weights/scale, global `Z12` origin, unit-bearing action, spacetime EOM, bridge completion, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3142/S2092 V_lift local EOM remains axiom-branch only", "## P3142/S2092 V_lift local EOM remains axiom-branch only\n\n`P3142/S2092` supplies an explicit local derivative and positive Hessian for the axiom selector potential, so the non-strict branch now has a toy field-variational lift.  The lift still imports `A_origin`, `A_lambda`, assumed weights/scale, and a local chart; it does not provide a strict global selector source, physical unit measure, spacetime stress/EOM coupling, `L_total`, bridge-completion theorem, role-transfer theorem, or ToE closure.\n")
    append_once(AGENTS, "Current axiom selector field-variational lift guardrail (P3142/S2092, 2026-07-13)", "## Current axiom selector field-variational lift guardrail (P3142/S2092, 2026-07-13)\n\n- P3142 constructs a non-strict local-chart field-variational lift `V_lift` of the P3141 axiom selector potential and verifies zero gradient plus positive Hessian at `(theta,s)=(0,+1)` across `27` positive parameter triples.\n- This gives a conditional symmetry-breaking engine only inside the explicit axiom branch; it still imports selector axioms, a local chart, weights, scale, and no unit-bearing measure.\n- Do not promote `V_lift` to strict `QW-2191` discharge, strict selector closure, global `Z12` origin, unit-bearing action/EOM, bridge completion, role transfer, `L_total`, or ToE closure.\n- Next honest move: audit exactly one candidate strict source for the `V_lift` weights/scale (`mu`, `w_theta`, `w_s`, or `kappa`) against import freedom, unit-bearing normalization, and global `Z12` quotient compatibility; otherwise preserve the P3140-P3142 non-strict axiom-branch boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
