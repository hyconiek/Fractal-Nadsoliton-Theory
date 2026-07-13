"""P3157/S2107: Omega_dim mass-unit torsor audit.

P3156 showed that alpha_geo-to-Higgs coupling formulas become dimensionally
valid only after a mass-squared scale is supplied.  This packet constructs the
missing mass-unit object explicitly as Omega_M, a positive mass-unit torsor, and
tests whether coupling it to alpha_geo fixes a canonical representative.  It
does not: all positive rescalings remain equally compatible unless a new strict
unit-source law breaks the torsor.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3157_s2107_omega_dim_mass_unit_torsor_audit.json"
MD = GEN / "p3157_s2107_omega_dim_mass_unit_torsor_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P3156": GEN / "p3156_s2106_alpha_geo_to_higgs_coupling_formula_audit.json",
    "P3155": GEN / "p3155_s2105_higgs_vev_state_source_audit.json",
    "P3116": GEN / "p3116_s2066_k_dim_dimension_source_functor_audit.json",
    "P3118": GEN / "p3118_s2068_r_dim_action_length_time_relation_audit.json",
    "P3119": GEN / "p3119_s2069_xi_lt_axis_source_object_audit.json",
}

S_VALUE = "alpha_geo"
SCALES = [Fraction(1, 2), Fraction(1), Fraction(2), Fraction(3)]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def fstr(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def torsor_rows() -> list[dict[str, Any]]:
    rows = []
    for c in SCALES:
        rows.append({
            "scale_c": fstr(c),
            "Omega_M_representative": f"{fstr(c)}*M_*",
            "mu2_formula": f"{S_VALUE}*({fstr(c)}*M_*)^2",
            "lambda_formula": S_VALUE,
            "v2_formula": f"({fstr(c)}*M_*)^2",
            "dimensionally_valid": True,
            "same_dimensionless_ratio": True,
            "canonical_representative_selected": c == 1,
            "strict_source_law_for_c": False,
        })
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [
        {"gate": "mass_dimension_carrier_constructed", "satisfied": True, "detail": "Omega_M supplies a formal mass-dimension carrier"},
        {"gate": "alpha_geo_coupling_dimensionally_valid", "satisfied": True, "detail": "mu2=alpha_geo*Omega_M^2 and lambda=alpha_geo are dimensionally valid"},
        {"gate": "positive_torsor_representative_unique", "satisfied": False, "detail": "Omega_M -> c Omega_M leaves dimensionless compatibility intact"},
        {"gate": "strict_nadsoliton_source_law_for_mass_scale", "satisfied": False, "detail": "no current artifact selects c or M_* nonconventionally"},
        {"gate": "action_length_time_coupling", "satisfied": False, "detail": "P3116/P3118/P3119 leave K_dim/R_dim/Xi_LT unit chain unexported"},
    ]


def build_payload() -> dict[str, Any]:
    rows = torsor_rows()
    gates = gate_rows()
    accepted = [r for r in rows if r["canonical_representative_selected"] and r["strict_source_law_for_c"]]
    counts = {
        "torsor_scale_rows": len(rows),
        "dimensionally_valid_rows": sum(r["dimensionally_valid"] for r in rows),
        "canonical_rows_by_gauge_label": sum(r["canonical_representative_selected"] for r in rows),
        "strict_source_selected_rows": len(accepted),
        "gate_rows": len(gates),
        "gates_satisfied": sum(r["satisfied"] for r in gates),
    }
    return {
        "status": "P3157_OMEGA_DIM_MASS_UNIT_TORSOR_FORMAL_CARRIER_NO_STRICT_SOURCE",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "constructed_object": {
            "name": "Omega_M positive mass-unit torsor",
            "classification": "formal_dimension_carrier_with_unfixed_positive_scale",
            "scope": "supplies the missing mass dimension for P3156 but tests canonical source selection under Omega_M -> c Omega_M",
        },
        "torsor_rows": rows,
        "gate_rows": gates,
        "finite_theorem": {
            "name": "P3157_T1_mass_unit_torsor_source_obstruction",
            "statement": "Adding a formal mass-unit torsor Omega_M makes the alpha_geo Higgs-coupling schema dimensionally valid: mu2=alpha_geo*Omega_M^2 and lambda=alpha_geo.  However, the positive rescaling Omega_M -> c Omega_M changes mu2 and v while preserving all dimensionless compatibility.  Current artifacts do not export a strict source law selecting c or coupling Omega_M to the action/length/time unit chain, so Omega_M is a useful missing-object schema but not a strict unit source.",
            "finite_counts": counts,
        },
        "decision": {
            "bounded_result": "P3157 constructs the missing mass-unit torsor explicitly and shows why it does not close the Higgs/SM/EH lane: formal dimension can be supplied, but canonical scale selection is still absent.",
            "next_honest_step": "The next honest proof-grade move is not another Higgs formula.  Construct exactly one strict source law for the positive mass/action unit torsor, preferably by returning to the P3116 Omega_dim/K_dim obligation with a new nonconventional source object; otherwise preserve the no-strict-unit boundary.",
            "negative_export_flags": {
                "strict_mass_unit_source_exported": False,
                "K_dim_functor_exported": False,
                "strict_higgs_couplings_exported": False,
                "strict_higgs_vev_exported": False,
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
    lines = ["# P3157/S2107 Omega_dim mass-unit torsor audit", "", f"Status: `{payload['status']}`", "", "## Constructed object", f"- `{payload['constructed_object']['name']}`", f"- Classification: `{payload['constructed_object']['classification']}`", f"- Scope: `{payload['constructed_object']['scope']}`", "", "## Finite theorem", f"`{th['name']}`: {th['statement']}", "", "## Finite counts"]
    for k, v in th["finite_counts"].items():
        lines.append(f"- `{k}`: `{v}`")
    lines.extend(["", "## Torsor rows"])
    for row in payload["torsor_rows"]:
        lines.append(f"- `c={row['scale_c']}`: Omega `{row['Omega_M_representative']}`, mu2 `{row['mu2_formula']}`, v2 `{row['v2_formula']}`, strict source `{row['strict_source_law_for_c']}`")
    lines.extend(["", "## Gates"])
    for row in payload["gate_rows"]:
        lines.append(f"- `{row['gate']}`: `{row['satisfied']}` — {row['detail']}")
    lines.extend(["", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3157/S2107 Omega_dim mass-unit torsor audit", "## P3157/S2107 Omega_dim mass-unit torsor audit\n\n`P3157/S2107` constructs `Omega_M`, a formal positive mass-unit torsor that supplies the missing mass dimension for the P3156 `alpha_geo` Higgs-coupling schema.  The schema `mu2=alpha_geo*Omega_M^2`, `lambda=alpha_geo` is dimensionally valid, but the orbit `Omega_M -> c Omega_M` leaves dimensionless compatibility intact while changing `mu2` and `v`.  No current strict source law selects `c` or couples `Omega_M` to the P3116/P3118/P3119 unit chain, so no strict mass unit, K_dim functor, Higgs coupling/VEV, EH coupling, `L_total`, SM/GR generation, selector closure, bridge/role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3157/S2107 Omega_M is formal not sourced", "## P3157/S2107 Omega_M is formal not sourced\n\n`P3157/S2107` shows that a formal mass unit would repair the dimensional defect in P3156, but without a strict source law for the positive torsor representative it is only a gauge/postulate.  The Higgs/SM/EH lane remains blocked by the same unit-source problem isolated by P3116.\n")
    append_once(AGENTS, "Current Omega_dim mass-unit torsor audit guardrail (P3157/S2107, 2026-07-13)", "## Current Omega_dim mass-unit torsor audit guardrail (P3157/S2107, 2026-07-13)\n\n- P3157 constructs `Omega_M`, a formal positive mass-unit torsor for the P3156 alpha-to-Higgs schema.\n- `mu2=alpha_geo*Omega_M^2` and `lambda=alpha_geo` are dimensionally valid, but `Omega_M -> c Omega_M` leaves compatibility intact and no strict source law selects `c`.\n- Do not promote `Omega_M` to a strict mass unit, K_dim functor, Higgs coupling/VEV source, EH coupling, unit-bearing `L_total`, strict SM/GR generation, selector closure, bridge/role transfer, or ToE.\n- Next honest move: construct exactly one strict source law for the positive mass/action unit torsor via the P3116 `Omega_dim/K_dim` obligation, or preserve the no-strict-unit boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
