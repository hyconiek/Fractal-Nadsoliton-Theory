"""P3150/S2100: source-selection audit for R_SM^ax hypercharge data.

P3149 left the honest next question: can the axiom-branch SM registry be
source-selected rather than merely installed?  This packet attacks exactly one
part of that problem, the hypercharge vector.  It derives the one-family SM
hypercharge *ray* from local Yukawa invariance plus anomaly constraints, then
records the remaining strict-source gaps: representation content, normalization,
and nadsoliton provenance are still not selected.
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
OUT = GEN / "p3150_s2100_rsm_hypercharge_ray_source_selection_audit.json"
MD = GEN / "p3150_s2100_rsm_hypercharge_ray_source_selection_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P3149": GEN / "p3149_s2099_brs_ltotal_interface_invariance_audit.json",
    "P3148": GEN / "p3148_s2098_sm_representation_registry_completion_audit.json",
    "P3147": GEN / "p3147_s2097_axiom_lie_smgr_fit_readiness_matrix.json",
    "P1962": GEN / "p1962_s912_strict_matter_higgs_brst_extension_registry_audit.json",
}

FIELDS = ["q", "u", "d", "l", "e", "h"]
STANDARD_RAY = {"q": sp.Rational(1, 6), "u": sp.Rational(-2, 3), "d": sp.Rational(1, 3), "l": sp.Rational(-1, 2), "e": sp.Integer(1), "h": sp.Rational(1, 2)}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def linear_system() -> tuple[list[sp.Symbol], list[dict[str, Any]], sp.Matrix]:
    q, u, d, l, e, h = sp.symbols("q u d l e h")
    vars_ = [q, u, d, l, e, h]
    equations = [
        ("up_yukawa", q + h + u, "Y(Q)+Y(H)+Y(u_c)=0"),
        ("down_yukawa", q - h + d, "Y(Q)-Y(H)+Y(d_c)=0"),
        ("lepton_yukawa", l - h + e, "Y(L)-Y(H)+Y(e_c)=0"),
        ("SU2_SU2_U1", 3 * q + l, "3Y(Q)+Y(L)=0"),
        ("gravity_gravity_U1", 6 * q + 3 * u + 3 * d + 2 * l + e, "sum dim(R)Y(R)=0"),
    ]
    rows = []
    matrix_rows = []
    for name, expr, text in equations:
        coeffs = [sp.expand(expr).coeff(v) for v in vars_]
        matrix_rows.append(coeffs)
        rows.append({"equation": name, "formula": text, "coefficients_q_u_d_l_e_h": [str(c) for c in coeffs]})
    return vars_, rows, sp.Matrix(matrix_rows)


def derive_ray() -> dict[str, Any]:
    vars_, rows, mat = linear_system()
    nullspace = mat.nullspace()
    generator = nullspace[0]
    q, u, d, l, e, h = vars_
    scale = sp.symbols("s")
    solution = {str(v): sp.simplify(scale * generator[i] / generator[0]) for i, v in enumerate(vars_)}
    # Normalize to h=1/2 to compare with P3148 convention.
    h_value = solution["h"]
    norm_factor = sp.Rational(1, 2) / h_value
    normalized = {k: sp.simplify(v * norm_factor) for k, v in solution.items()}
    su3 = 2 * q + u + d
    cubic = 6 * q**3 + 3 * u**3 + 3 * d**3 + 2 * l**3 + e**3
    substit = {q: solution["q"], u: solution["u"], d: solution["d"], l: solution["l"], e: solution["e"], h: solution["h"]}
    su3_on_ray = sp.simplify(su3.subs(substit))
    cubic_on_ray = sp.factor(cubic.subs(substit))
    matches_standard = all(sp.simplify(normalized[k] - STANDARD_RAY[k]) == 0 for k in FIELDS)
    return {
        "equation_rows": rows,
        "matrix_rank": int(mat.rank()),
        "unknown_count": len(vars_),
        "nullity": len(nullspace),
        "ray_generator_q_u_d_l_e_h": [str(x) for x in generator],
        "general_solution_with_q_scale": {k: str(v) for k, v in solution.items()},
        "normalized_h_equals_one_half": {k: str(v) for k, v in normalized.items()},
        "matches_P3148_standard_hypercharges": matches_standard,
        "redundant_SU3_SU3_U1_on_ray": str(su3_on_ray),
        "cubic_U1_anomaly_on_ray": str(cubic_on_ray),
        "redundant_checks_vanish": su3_on_ray == 0 and cubic_on_ray == 0,
    }


def candidate_source_rows(ray: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "candidate": "Y_SM^ray consistency constraints",
            "selects_hypercharge_ray": ray["matches_P3148_standard_hypercharges"] and ray["redundant_checks_vanish"],
            "selects_absolute_normalization": False,
            "selects_representation_content": False,
            "strict_nadsoliton_source_law": False,
            "noncircular_current_artifact": False,
            "verdict": "conditional ray witness only: uses SM-like field/Yukawa/anomaly obligations",
        },
        {
            "candidate": "P3146 unit/action axioms",
            "selects_hypercharge_ray": False,
            "selects_absolute_normalization": False,
            "selects_representation_content": False,
            "strict_nadsoliton_source_law": False,
            "noncircular_current_artifact": True,
            "verdict": "dimension scales do not select gauge charges",
        },
        {
            "candidate": "P1960/P1961 local Lie/BRST algebra",
            "selects_hypercharge_ray": False,
            "selects_absolute_normalization": False,
            "selects_representation_content": False,
            "strict_nadsoliton_source_law": False,
            "noncircular_current_artifact": True,
            "verdict": "local algebra supplies consistency rules, not a unique matter registry",
        },
        {
            "candidate": "P3148 installed R_SM^ax registry",
            "selects_hypercharge_ray": True,
            "selects_absolute_normalization": True,
            "selects_representation_content": True,
            "strict_nadsoliton_source_law": False,
            "noncircular_current_artifact": False,
            "verdict": "self-referential as a source: it is the installed target registry",
        },
    ]


def build_payload() -> dict[str, Any]:
    ray = derive_ray()
    candidates = candidate_source_rows(ray)
    accepted_strict = [r for r in candidates if r["selects_hypercharge_ray"] and r["selects_absolute_normalization"] and r["selects_representation_content"] and r["strict_nadsoliton_source_law"] and r["noncircular_current_artifact"]]
    counts = {
        "linear_equation_rows": len(ray["equation_rows"]),
        "unknowns": ray["unknown_count"],
        "matrix_rank": ray["matrix_rank"],
        "nullity": ray["nullity"],
        "ray_witnesses_matching_P3148": int(ray["matches_P3148_standard_hypercharges"]),
        "redundant_anomaly_checks_vanishing": int(ray["redundant_checks_vanish"]),
        "candidate_source_rows": len(candidates),
        "strict_accepted_source_rows": len(accepted_strict),
    }
    return {
        "status": "P3150_RSM_HYPERCHARGE_RAY_SOURCE_SELECTION_CONDITIONAL_RAY_NO_STRICT_SOURCE",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "constructed_object": {
            "name": "Y_SM^ray hypercharge consistency ray",
            "classification": "conditional_source_selection_witness_for_hypercharge_ratios_only",
            "scope": "derives the P3148 hypercharge ratios from Yukawa plus anomaly equations, not the absolute normalization or representation content",
        },
        "hypercharge_ray_derivation": ray,
        "candidate_source_rows": candidates,
        "finite_theorem": {
            "name": "P3150_T1_hypercharge_ray_derivation_and_source_obstruction",
            "statement": "The one-family Yukawa-invariance plus anomaly equations form a 5 x 6 linear system of rank 5, so they select a one-dimensional hypercharge ray.  Normalizing the Higgs charge to 1/2 gives exactly the P3148 assignments q=1/6, u=-2/3, d=1/3, l=-1/2, e=1, h=1/2; the SU3^2 U1 and U1^3 anomaly checks vanish on the same ray.  This is a genuine conditional source-selection witness for hypercharge ratios, but it is not a strict nadsoliton source: it assumes the SM-like field content/Yukawa obligations, leaves absolute normalization free until a unit/charge convention is imposed, and does not select the representation content itself.",
            "finite_counts": counts,
        },
        "decision": {
            "bounded_result": "P3150 improves the P3148/P3149 axiom branch: the hypercharge ratios are not arbitrary once the one-family field pattern, Yukawa terms, and anomaly constraints are assumed; they are the unique ray of the finite system.",
            "why_not_strict": "The witness is conditional on an installed SM-like field pattern and Yukawa/anomaly obligations.  No current strict object selects that field pattern, the Higgs normalization, the charge unit, or the registry from nadsoliton data.",
            "next_honest_step": "Construct P3151 as a finite source-selection audit for the representation content itself: test whether any strict object selects the five chiral multiplet pattern plus Higgs dimensions (3,2), (bar3,1), (bar3,1), (1,2), (1,1), (1,2) without importing SM family data.  If no source exists, preserve R_SM^ax/Y_SM^ray as conditional phenomenology and pivot to unit-charge normalization or GR/EH nonproxy coupling.",
            "negative_export_flags": {
                "strict_registry_source_exported": False,
                "absolute_charge_normalization_exported": False,
                "representation_content_source_exported": False,
                "unit_bearing_L_total_exported": False,
                "global_BV_BRST_exported": False,
                "strict_SM_generation_exported": False,
                "strict_GR_generation_exported": False,
                "ToE_closure_exported": False,
            },
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    th = payload["finite_theorem"]
    ray = payload["hypercharge_ray_derivation"]
    lines = [
        "# P3150/S2100 R_SM hypercharge-ray source-selection audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Constructed object",
        f"- `{payload['constructed_object']['name']}`",
        f"- Classification: `{payload['constructed_object']['classification']}`",
        f"- Scope: `{payload['constructed_object']['scope']}`",
        "",
        "## Finite theorem",
        f"`{th['name']}`: {th['statement']}",
        "",
        "## Finite counts",
    ]
    for key, value in th["finite_counts"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Normalized ray"])
    for key, value in ray["normalized_h_equals_one_half"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Candidate source rows"])
    for row in payload["candidate_source_rows"]:
        lines.append(f"- `{row['candidate']}`: ray `{row['selects_hypercharge_ray']}`, normalization `{row['selects_absolute_normalization']}`, representation content `{row['selects_representation_content']}`, strict source `{row['strict_nadsoliton_source_law']}`; {row['verdict']}")
    lines.extend(["", "## Decision", payload["decision"]["bounded_result"], "", "## Why this is not strict", payload["decision"]["why_not_strict"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3150/S2100 R_SM hypercharge-ray source-selection audit", "## P3150/S2100 R_SM hypercharge-ray source-selection audit\n\n`P3150/S2100` constructs `Y_SM^ray`, a conditional hypercharge source-selection witness for `R_SM^ax`.  The Yukawa/anomaly equation system has `5` rows, `6` unknowns, rank `5`, and nullity `1`; after normalizing `Y(H)=1/2`, it reproduces the P3148 hypercharges exactly and the redundant `SU3^2 U1` plus `U1^3` anomaly checks vanish.  This selects only a conditional hypercharge ray: it does not source the SM field content, absolute charge unit, unit-bearing `L_total`, global BV/BRST, GR/EH coupling, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3150/S2100 hypercharge ray is conditional", "## P3150/S2100 hypercharge ray is conditional\n\n`P3150/S2100` proves that the installed one-family field/Yukawa/anomaly pattern fixes the hypercharge ratios up to scale.  It does not derive that field pattern from the strict kernel, fix the absolute charge normalization from a physical unit theorem, or install the result in a unit-bearing Lagrangian/EOM.\n")
    append_once(AGENTS, "Current R_SM hypercharge-ray source-selection guardrail (P3150/S2100, 2026-07-13)", "## Current R_SM hypercharge-ray source-selection guardrail (P3150/S2100, 2026-07-13)\n\n- P3150 constructs `Y_SM^ray`, a conditional source-selection witness for the hypercharge ratios in `R_SM^ax`.\n- The exact Yukawa/anomaly system has `5` equations, `6` unknowns, rank `5`, and nullity `1`; with `Y(H)=1/2` it reproduces the P3148 hypercharges, and redundant `SU3^2 U1` plus `U1^3` anomaly checks vanish.\n- This is real algebraic progress but not strict source closure: it assumes the SM-like field/Yukawa/anomaly pattern, leaves absolute charge normalization conventional, and does not source representation content from the nadsoliton.\n- Do not promote `Y_SM^ray` to strict SM generation, unit-bearing `L_total`, global BV/BRST, GR/EH closure, bridge/role transfer, or ToE.\n- Next honest move: audit whether any strict object selects the representation content itself — the five chiral multiplets plus Higgs dimensions — without importing SM family data; otherwise preserve this as conditional phenomenology and pivot to charge-unit normalization or GR/EH nonproxy coupling.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
