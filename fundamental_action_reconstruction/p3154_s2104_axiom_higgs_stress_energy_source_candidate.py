"""P3154/S2104: symbolic stress-energy candidate for the axiom branch.

P3153 identified the next missing GR/EH object: one stress-energy source
candidate built from the P3149 matter/Higgs local Lagrangian convention.  This
packet constructs the homogeneous Higgs-scalar stress-energy tensor on a flat
FRW receiver, derives the exact conservation residual, and audits whether the
current repo exports the state/VEV/unit data needed to use it as a strict EH
source.  The result is a real symbolic T_mu_nu candidate, but not a strict
source closure.
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
OUT = GEN / "p3154_s2104_axiom_higgs_stress_energy_source_candidate.json"
MD = GEN / "p3154_s2104_axiom_higgs_stress_energy_source_candidate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P3153": GEN / "p3153_s2103_axiom_gr_eh_nonproxy_coupling_audit.json",
    "P3149": GEN / "p3149_s2099_brs_ltotal_interface_invariance_audit.json",
    "P3148": GEN / "p3148_s2098_sm_representation_registry_completion_audit.json",
    "P3146": GEN / "p3146_s2096_axiom_unit_action_measure_bridge.json",
    "P3152": GEN / "p3152_s2102_y_sm_charge_unit_normalization_obstruction.json",
}

t = sp.symbols("t", positive=True)
p, k, v0, V0, mu2, lam = sp.symbols("p k v0 V0 mu2 lam", nonzero=True)
h = sp.Function("h")(t)
a = t**p
Hubble = sp.diff(a, t) / a
V = -mu2 * h**2 / 2 + lam * h**4 / 4


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def scalar_stress_energy() -> dict[str, sp.Expr]:
    hdot = sp.diff(h, t)
    rho = sp.simplify(hdot**2 / 2 + V)
    pressure = sp.simplify(hdot**2 / 2 - V)
    conservation = sp.simplify(sp.diff(rho, t) + 3 * Hubble * (rho + pressure))
    kg = sp.simplify(sp.diff(h, t, 2) + 3 * Hubble * hdot + sp.diff(V, h))
    factorized = sp.factor(conservation - hdot * kg)
    return {
        "rho": rho,
        "pressure": pressure,
        "T_00": rho,
        "T_ii_covariant_each": sp.simplify(a**2 * pressure),
        "conservation_residual": conservation,
        "kg_operator": kg,
        "conservation_minus_hdot_times_kg": factorized,
    }


def candidate_state_rows() -> list[dict[str, Any]]:
    data = scalar_stress_energy()
    rows = []
    substitutions = [
        ("zero_higgs_state", {h: 0, sp.diff(h, t): 0, sp.diff(h, t, 2): 0}, "trivial zero stress-energy if V(0)=0; not a matter/GR source"),
        ("constant_imported_vev_v0", {h: v0, sp.diff(h, t): 0, sp.diff(h, t, 2): 0}, "constant condensate is conserved but v0 and potential scale are imported"),
        ("rolling_log_profile_free_potential", {h: k * sp.log(t), sp.diff(h, t): k / t, sp.diff(h, t, 2): -k / t**2, mu2: 0, lam: 0}, "dimensionless rolling profile has nonzero conservation residual except p=1/3"),
    ]
    for name, subs, note in substitutions:
        rho = sp.simplify(data["rho"].subs(subs))
        pressure = sp.simplify(data["pressure"].subs(subs))
        residual = sp.simplify(data["conservation_residual"].subs(subs))
        kg_residual = sp.simplify(data["kg_operator"].subs(subs))
        rows.append({
            "state_candidate": name,
            "rho": str(rho),
            "pressure": str(pressure),
            "conservation_residual": str(residual),
            "kg_residual": str(kg_residual),
            "conserved_without_extra_condition": bool(sp.simplify(residual) == 0),
            "nonzero_stress_energy": bool(sp.simplify(rho) != 0 or sp.simplify(pressure) != 0),
            "strict_state_or_vev_source_exported": False,
            "unit_dimension_source_exported": False,
            "note": note,
        })
    return rows


def source_gate_rows() -> list[dict[str, Any]]:
    return [
        {"gate": "local_P3149_Higgs_lagrangian_form", "satisfied": True, "detail": "local Higgs kinetic/potential form can define a symbolic T_mu_nu candidate"},
        {"gate": "noncircular_Higgs_state_or_VEV", "satisfied": False, "detail": "no strict state, condensate, or VEV source is exported"},
        {"gate": "dimensionful_unit_and_Newton_coupling", "satisfied": False, "detail": "P3146 units are axiomatic and do not fix G_N or Higgs potential units strictly"},
        {"gate": "metric_variation_bundle", "satisfied": False, "detail": "candidate is FRW-reduced; no full metric bundle variation theorem is exported here"},
        {"gate": "covariant_conservation_on_sourced_solution", "satisfied": False, "detail": "conservation requires the KG equation/state dynamics; not sourced noncircularly"},
    ]


def build_payload() -> dict[str, Any]:
    tensor = scalar_stress_energy()
    states = candidate_state_rows()
    gates = source_gate_rows()
    accepted = [r for r in states if r["conserved_without_extra_condition"] and r["nonzero_stress_energy"] and r["strict_state_or_vev_source_exported"] and r["unit_dimension_source_exported"]]
    counts = {
        "symbolic_tensor_components": 2,
        "conservation_identity_residual_zero": int(sp.simplify(tensor["conservation_minus_hdot_times_kg"]) == 0),
        "state_candidate_rows": len(states),
        "conserved_state_rows_without_extra_condition": sum(r["conserved_without_extra_condition"] for r in states),
        "nonzero_stress_energy_rows": sum(r["nonzero_stress_energy"] for r in states),
        "source_gate_rows": len(gates),
        "source_gates_satisfied": sum(r["satisfied"] for r in gates),
        "accepted_strict_stress_energy_sources": len(accepted),
    }
    return {
        "status": "P3154_AXIOM_HIGGS_STRESS_ENERGY_SYMBOLIC_CANDIDATE_NO_STRICT_SOURCE",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "constructed_object": {
            "name": "T_H^ax homogeneous Higgs stress-energy candidate",
            "classification": "symbolic_stress_energy_candidate_with_source_gap",
            "scope": "P3149 local Higgs scalar form on flat FRW receiver; conservation and coupling-source audit",
        },
        "symbolic_stress_energy": {key: str(value) for key, value in tensor.items()},
        "state_candidate_rows": states,
        "source_gate_rows": gates,
        "finite_theorem": {
            "name": "P3154_T1_higgs_stress_energy_candidate_and_source_obstruction",
            "statement": "The P3149 Higgs local Lagrangian form yields a valid symbolic homogeneous stress-energy candidate on flat FRW: rho=1/2 hdot^2+V(h), p=1/2 hdot^2-V(h), and covariant conservation reduces exactly to hdot times the Higgs/KG equation.  This is a real GR/EH source candidate object, but current artifacts do not export a noncircular Higgs state/VEV, dimensionful potential/Newton units, or full metric-variation bundle; therefore it cannot close EH coupling or L_total.",
            "finite_counts": counts,
        },
        "decision": {
            "bounded_result": "P3154 advances P3153 by constructing the actual symbolic T_mu_nu candidate, but the source gates remain open: the candidate needs a strict state/VEV and unit/coupling theorem before it can source nonflat EH rows.",
            "why_not_strict": "The zero state is conserved but trivial, the constant VEV row imports v0/potential scale, and the rolling profile is not conserved without extra dynamics.  None supplies noncircular strict state provenance or units.",
            "next_honest_step": "Construct P3155 as one VEV/state-source audit: test whether any current strict scalar invariant fixes a nonzero Higgs condensate value and potential scale without importing SM phenomenology.  If no such source is found, freeze the axiom-branch EH stress-energy route as conditional and pivot back to strict unit/action or selector-source intake.",
            "negative_export_flags": {
                "strict_higgs_state_or_vev_exported": False,
                "strict_stress_energy_source_exported": False,
                "newton_constant_source_exported": False,
                "einstein_hilbert_coupling_exported": False,
                "unit_bearing_L_total_exported": False,
                "strict_SM_generation_exported": False,
                "strict_GR_generation_exported": False,
                "selector_closure_exported": False,
                "bridge_or_role_transfer_exported": False,
                "ToE_closure_exported": False,
            },
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    th = payload["finite_theorem"]
    lines = ["# P3154/S2104 axiom Higgs stress-energy source candidate", "", f"Status: `{payload['status']}`", "", "## Constructed object", f"- `{payload['constructed_object']['name']}`", f"- Classification: `{payload['constructed_object']['classification']}`", f"- Scope: `{payload['constructed_object']['scope']}`", "", "## Finite theorem", f"`{th['name']}`: {th['statement']}", "", "## Finite counts"]
    for k, v in th["finite_counts"].items():
        lines.append(f"- `{k}`: `{v}`")
    lines.extend(["", "## Symbolic stress-energy"])
    for k, v in payload["symbolic_stress_energy"].items():
        lines.append(f"- `{k}`: `{v}`")
    lines.extend(["", "## State candidate rows"])
    for row in payload["state_candidate_rows"]:
        lines.append(f"- `{row['state_candidate']}`: rho `{row['rho']}`, p `{row['pressure']}`, conservation `{row['conservation_residual']}`, conserved `{row['conserved_without_extra_condition']}`, nonzero `{row['nonzero_stress_energy']}`; {row['note']}")
    lines.extend(["", "## Source gates"])
    for row in payload["source_gate_rows"]:
        lines.append(f"- `{row['gate']}`: `{row['satisfied']}` — {row['detail']}")
    lines.extend(["", "## Decision", payload["decision"]["bounded_result"], "", "## Why this is not strict", payload["decision"]["why_not_strict"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3154/S2104 axiom Higgs stress-energy source candidate", "## P3154/S2104 axiom Higgs stress-energy source candidate\n\n`P3154/S2104` constructs `T_H^ax`, a symbolic homogeneous Higgs stress-energy candidate from the P3149 local matter/Higgs Lagrangian form on a flat FRW receiver.  The audit derives `rho=hdot^2/2+V(h)`, `p=hdot^2/2-V(h)`, and verifies that covariant conservation reduces exactly to `hdot` times the Higgs/KG equation.  This is a real candidate object, but no noncircular Higgs state/VEV, dimensionful potential/Newton unit, or full metric-variation bundle is exported; no EH coupling closure, unit-bearing `L_total`, strict SM/GR generation, bridge/role transfer, selector closure, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3154/S2104 Higgs Tmunu remains conditional", "## P3154/S2104 Higgs Tmunu remains conditional\n\n`P3154/S2104` turns the local P3149 Higgs receiver into a symbolic stress-energy tensor candidate.  The tensor is conserved only on a sourced Higgs/KG state, and current artifacts do not export the required nonzero condensate/state, potential scale, Newton coupling, or full metric variation theorem.\n")
    append_once(AGENTS, "Current axiom Higgs stress-energy candidate guardrail (P3154/S2104, 2026-07-13)", "## Current axiom Higgs stress-energy candidate guardrail (P3154/S2104, 2026-07-13)\n\n- P3154 constructs `T_H^ax`, a symbolic homogeneous Higgs stress-energy candidate from the P3149 local Higgs/matter Lagrangian interface.\n- The exact conservation identity holds only as `nabla_mu T^{mu0}=hdot*(KG_h)`, so conservation requires a Higgs/KG state equation rather than following from the local receiver alone.\n- Current artifacts do not export a noncircular Higgs state/VEV, dimensionful potential/Newton unit, or full metric-variation bundle.\n- Do not promote `T_H^ax` to EH coupling closure, unit-bearing `L_total`, strict SM/GR generation, selector closure, bridge/role transfer, or ToE.\n- Next honest move: P3155 should audit exactly one VEV/state-source candidate from current strict scalar invariants; otherwise freeze the axiom-branch EH stress-energy route as conditional.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
