"""P3156/S2106: alpha_geo-to-Higgs-coupling formula audit.

P3155 froze the generic scalar inventory and allowed continuation only for one
new formula mapping a strict scalar invariant to both Higgs couplings with units.
This packet tests exactly one schema: Phi_alpha^H uses the strict positive
scalar S=alpha_geo=4 ln 2 as the dimensionless seed for (mu2, lambda).  The
audit separates ratio selection from dimensionful mass-source import.
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
OUT = GEN / "p3156_s2106_alpha_geo_to_higgs_coupling_formula_audit.json"
MD = GEN / "p3156_s2106_alpha_geo_to_higgs_coupling_formula_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P3155": GEN / "p3155_s2105_higgs_vev_state_source_audit.json",
    "P3154": GEN / "p3154_s2104_axiom_higgs_stress_energy_source_candidate.json",
    "P3109": GEN / "p3109_s2059_alpha_geo_scale_orbit_quotient_source_law_audit.json",
    "P3116": GEN / "p3116_s2066_k_dim_dimension_source_functor_audit.json",
    "ALPHA": GEN / "alpha_geo_strict_derived_v1.json",
}

S, M2 = sp.symbols("S M2", positive=True)
# Dimension vector convention in four-dimensional natural-unit scalar theory:
# [mass_power]. lambda is dimensionless, mu2 has mass dimension 2, v has mass 1.
DIMENSION = {"dimensionless": 0, "mass2": 2}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def formula_rows() -> list[dict[str, Any]]:
    rows = [
        {
            "row": "raw_dimensionless_alpha",
            "mu2_formula": S,
            "lambda_formula": S,
            "mu2_dimension": DIMENSION["dimensionless"],
            "lambda_dimension": DIMENSION["dimensionless"],
            "imports_mass_scale": False,
            "comment": "uses alpha_geo directly for both couplings",
        },
        {
            "row": "alpha_ratio_with_imported_mass_scale",
            "mu2_formula": M2 * S,
            "lambda_formula": S,
            "mu2_dimension": DIMENSION["mass2"],
            "lambda_dimension": DIMENSION["dimensionless"],
            "imports_mass_scale": True,
            "comment": "dimensionally valid only after importing an external mass^2 unit M2",
        },
        {
            "row": "alpha_squared_lambda_with_imported_mass_scale",
            "mu2_formula": M2 * S,
            "lambda_formula": S**2,
            "mu2_dimension": DIMENSION["mass2"],
            "lambda_dimension": DIMENSION["dimensionless"],
            "imports_mass_scale": True,
            "comment": "changes dimensionless ratio but still imports M2",
        },
    ]
    audited = []
    for row in rows:
        mu2 = row["mu2_formula"]
        lam = row["lambda_formula"]
        v2 = sp.simplify(mu2 / lam)
        dimensionally_valid = row["mu2_dimension"] == 2 and row["lambda_dimension"] == 0
        strict_unit_source = dimensionally_valid and not row["imports_mass_scale"]
        audited.append({
            **{k: (str(v) if isinstance(v, sp.Expr) else v) for k, v in row.items()},
            "v2_formula": str(v2),
            "v2_dimension": row["mu2_dimension"] - row["lambda_dimension"],
            "dimensionally_valid_higgs_couplings": dimensionally_valid,
            "strict_mass_unit_source_exported": strict_unit_source,
            "accepted_strict_formula": dimensionally_valid and strict_unit_source,
        })
    return audited


def build_payload() -> dict[str, Any]:
    rows = formula_rows()
    accepted = [r for r in rows if r["accepted_strict_formula"]]
    counts = {
        "formula_schema_count": 1,
        "formula_variant_rows": len(rows),
        "dimensionally_valid_rows": sum(r["dimensionally_valid_higgs_couplings"] for r in rows),
        "rows_importing_mass_scale": sum(r["imports_mass_scale"] for r in rows),
        "accepted_strict_formula_rows": len(accepted),
    }
    return {
        "status": "P3156_ALPHA_GEO_TO_HIGGS_COUPLING_FORMULA_AUDIT_NO_STRICT_UNIT_SOURCE",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "constructed_object": {
            "name": "Phi_alpha^H alpha_geo-to-Higgs-coupling formula schema",
            "classification": "single_formula_schema_dimension_and_source_audit",
            "scope": "tests whether S=alpha_geo can source both mu2 and lambda with units",
        },
        "formula_rows": rows,
        "finite_theorem": {
            "name": "P3156_T1_alpha_geo_higgs_coupling_unit_obstruction",
            "statement": "The single audited schema Phi_alpha^H can make dimensionless Higgs-coupling ratios from S=alpha_geo, but the Higgs mass parameter mu2 has mass dimension two.  Raw alpha_geo formulas are dimensionally invalid for mu2; dimensionally valid variants require an imported mass^2 scale M2.  Current artifacts therefore do not export a strict scalar-to-Higgs-coupling source for both mu2 and lambda with units.",
            "finite_counts": counts,
        },
        "decision": {
            "bounded_result": "P3156 satisfies the P3155 continuation rule by auditing exactly one new formula schema.  It fails as a strict source because alpha_geo is dimensionless and the only dimensionally valid variants import M2.",
            "next_honest_step": "Do not continue the axiom-branch Higgs/SM/EH lane without a genuinely new strict mass-unit source.  The next proof-grade move should pivot to the independent Omega_dim/K_dim unit-source obligation from P3116, or to a selector-source intake if a new non-premise source object is supplied.",
            "negative_export_flags": {
                "strict_higgs_couplings_exported": False,
                "strict_mass_unit_exported": False,
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
    lines = ["# P3156/S2106 alpha_geo to Higgs-coupling formula audit", "", f"Status: `{payload['status']}`", "", "## Constructed object", f"- `{payload['constructed_object']['name']}`", f"- Classification: `{payload['constructed_object']['classification']}`", f"- Scope: `{payload['constructed_object']['scope']}`", "", "## Finite theorem", f"`{th['name']}`: {th['statement']}", "", "## Finite counts"]
    for k, v in th["finite_counts"].items():
        lines.append(f"- `{k}`: `{v}`")
    lines.extend(["", "## Formula rows"])
    for row in payload["formula_rows"]:
        lines.append(f"- `{row['row']}`: mu2 `{row['mu2_formula']}` dim `{row['mu2_dimension']}`, lambda `{row['lambda_formula']}` dim `{row['lambda_dimension']}`, v2 `{row['v2_formula']}`, valid `{row['dimensionally_valid_higgs_couplings']}`, imports M2 `{row['imports_mass_scale']}`, accepted `{row['accepted_strict_formula']}`")
    lines.extend(["", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3156/S2106 alpha_geo to Higgs-coupling formula audit", "## P3156/S2106 alpha_geo to Higgs-coupling formula audit\n\n`P3156/S2106` audits exactly one new formula schema, `Phi_alpha^H`, mapping the strict dimensionless scalar `alpha_geo` to Higgs couplings.  Raw alpha formulas can set dimensionless ratios but cannot source `mu2`, which has mass dimension two; dimensionally valid variants require importing a mass scale `M2`.  Thus no strict Higgs coupling source, VEV source, EH coupling, unit-bearing `L_total`, strict SM/GR generation, selector closure, bridge/role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3156/S2106 alpha_geo Higgs formula remains unit-blocked", "## P3156/S2106 alpha_geo Higgs formula remains unit-blocked\n\n`P3156/S2106` shows that `alpha_geo` can only provide dimensionless Higgs-coupling ratios unless a mass-squared unit is imported.  The axiom-branch Higgs/SM/EH lane therefore remains blocked by the independent strict mass/action-unit source problem.\n")
    append_once(AGENTS, "Current alpha_geo to Higgs-coupling formula audit guardrail (P3156/S2106, 2026-07-13)", "## Current alpha_geo to Higgs-coupling formula audit guardrail (P3156/S2106, 2026-07-13)\n\n- P3156 audits exactly one new formula schema `Phi_alpha^H` from `alpha_geo` to Higgs couplings.\n- Raw `alpha_geo` formulas are dimensionally invalid for `mu2`; dimensionally valid variants import an external mass-squared scale `M2`.\n- Do not continue the axiom-branch Higgs/SM/EH lane without a genuinely new strict mass-unit source.\n- Do not promote P3156 to strict Higgs couplings, VEV, EH coupling, unit-bearing `L_total`, strict SM/GR generation, selector closure, bridge/role transfer, or ToE.\n- Next honest move: pivot to the independent `Omega_dim/K_dim` unit-source obligation from P3116, or to selector-source intake only if a new non-premise source object is supplied.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
