"""P3152/S2102: charge-unit normalization obstruction for Y_SM^ray.

P3150 selected the hypercharge ratios only up to scale and P3151 showed the
representation content is still installed rather than strict-sourced.  This
packet attacks the remaining finite normalization question: can current strict
artifacts fix the convention Y(H)=1/2 without importing the SM electric-charge
unit?  The answer is a bounded no-go: all homogeneous Yukawa/anomaly equations
and local BRST compatibility gates are invariant under a nonzero rational scale
torsor Y -> cY.
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
OUT = GEN / "p3152_s2102_y_sm_charge_unit_normalization_obstruction.json"
MD = GEN / "p3152_s2102_y_sm_charge_unit_normalization_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P3151": GEN / "p3151_s2101_rsm_representation_content_source_selection_audit.json",
    "P3150": GEN / "p3150_s2100_rsm_hypercharge_ray_source_selection_audit.json",
    "P3149": GEN / "p3149_s2099_brs_ltotal_interface_invariance_audit.json",
    "P3146": GEN / "p3146_s2096_axiom_unit_action_measure_bridge.json",
    "P3116": GEN / "p3116_s2066_k_dim_dimension_source_functor_audit.json",
}

BASE_Y = {
    "q": Fraction(1, 6),
    "u": Fraction(-2, 3),
    "d": Fraction(1, 3),
    "l": Fraction(-1, 2),
    "e": Fraction(1, 1),
    "h": Fraction(1, 2),
}
SCALES = [Fraction(-3), Fraction(-2), Fraction(-1), Fraction(-1, 2), Fraction(1, 2), Fraction(1), Fraction(2), Fraction(3)]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def fstr(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def scaled_y(scale: Fraction) -> dict[str, Fraction]:
    return {k: scale * v for k, v in BASE_Y.items()}


def equation_residuals(y: dict[str, Fraction]) -> dict[str, Fraction]:
    q, u, d, l, e, h = (y[k] for k in ["q", "u", "d", "l", "e", "h"])
    return {
        "up_yukawa": q + h + u,
        "down_yukawa": q - h + d,
        "lepton_yukawa": l - h + e,
        "SU2_SU2_U1": 3 * q + l,
        "gravity_gravity_U1": 6 * q + 3 * u + 3 * d + 2 * l + e,
        "SU3_SU3_U1_redundant": 2 * q + u + d,
        "U1_U1_U1_redundant": 6 * q**3 + 3 * u**3 + 3 * d**3 + 2 * l**3 + e**3,
    }


def electric_charges(y: dict[str, Fraction]) -> dict[str, Fraction]:
    # Left-handed Weyl convention: Q doublet components have T3=+/-1/2;
    # singlets use T3=0.  These values demonstrate scale dependence of charge
    # units rather than selecting a physical convention.
    return {
        "Q_up_component": Fraction(1, 2) + y["q"],
        "Q_down_component": Fraction(-1, 2) + y["q"],
        "u_c": y["u"],
        "d_c": y["d"],
        "L_neutrino_component": Fraction(1, 2) + y["l"],
        "L_electron_component": Fraction(-1, 2) + y["l"],
        "e_c": y["e"],
        "H_up_component": Fraction(1, 2) + y["h"],
        "H_down_component": Fraction(-1, 2) + y["h"],
    }


def gcd_fraction(values: list[Fraction]) -> Fraction:
    den_lcm = 1
    for v in values:
        den_lcm = den_lcm * v.denominator // __import__("math").gcd(den_lcm, v.denominator)
    nums = [abs(v.numerator * (den_lcm // v.denominator)) for v in values if v]
    if not nums:
        return Fraction(0)
    g = nums[0]
    for n in nums[1:]:
        g = __import__("math").gcd(g, n)
    return Fraction(g, den_lcm)


def scale_orbit_rows() -> list[dict[str, Any]]:
    rows = []
    for c in SCALES:
        y = scaled_y(c)
        residuals = equation_residuals(y)
        charges = electric_charges(y)
        rows.append({
            "scale_c": fstr(c),
            "Y_H": fstr(y["h"]),
            "all_yukawa_and_anomaly_residuals_zero": all(v == 0 for v in residuals.values()),
            "residuals": {k: fstr(v) for k, v in residuals.items()},
            "electric_charge_values": {k: fstr(v) for k, v in charges.items()},
            "hypercharge_unit_gcd": fstr(gcd_fraction(list(y.values()))),
            "electric_charge_gcd": fstr(gcd_fraction(list(charges.values()))),
        })
    return rows


def build_payload() -> dict[str, Any]:
    rows = scale_orbit_rows()
    invariant_count = sum(row["all_yukawa_and_anomaly_residuals_zero"] for row in rows)
    candidates = [
        {"candidate": "Y_SM^ray homogeneous Yukawa/anomaly system", "fixes_scale": False, "strict_source_law": False, "noncircular": True, "verdict": "equations select only a projective ray; every sampled nonzero rational scale remains valid"},
        {"candidate": "P3148/P3149 installed SM electric-charge convention", "fixes_scale": True, "strict_source_law": False, "noncircular": False, "verdict": "sets Y(H)=1/2 by convention/installed registry, circular for strict source selection"},
        {"candidate": "P3146 axiom unit/action bridge", "fixes_scale": False, "strict_source_law": False, "noncircular": True, "verdict": "length/time/action postulates do not define the U(1) charge unit"},
        {"candidate": "P3116 dimension-source functor inventory", "fixes_scale": False, "strict_source_law": False, "noncircular": True, "verdict": "current physical-unit functors remain unsourced and do not export a charge valuation"},
    ]
    accepted = [r for r in candidates if r["fixes_scale"] and r["strict_source_law"] and r["noncircular"]]
    counts = {
        "sampled_nonzero_rational_scales": len(rows),
        "invariant_scales": invariant_count,
        "distinct_Y_H_values": len({row["Y_H"] for row in rows}),
        "candidate_source_rows": len(candidates),
        "strict_accepted_source_rows": len(accepted),
    }
    return {
        "status": "P3152_YSM_CHARGE_UNIT_NORMALIZATION_SCALE_TORSOR_NO_STRICT_SOURCE",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "constructed_object": {"name": "U_Y^tors charge-unit normalization torsor", "classification": "bounded_scale_torsor_obstruction", "scope": "nonzero rational rescalings of the P3150 hypercharge ray, audited against Yukawa/anomaly residuals and local electric-charge values"},
        "scale_orbit_rows": rows,
        "candidate_source_rows": candidates,
        "finite_theorem": {"name": "P3152_T1_hypercharge_charge_unit_scale_torsor_obstruction", "statement": "The P3150 hypercharge equations and redundant anomaly checks are homogeneous in Y.  On the sampled nonzero rational scale orbit, every row has zero Yukawa/anomaly residuals while Y(H), hypercharge gcd, and electric-charge values vary.  Therefore current algebraic compatibility fixes only a projective hypercharge ray, not the physical charge unit or the convention Y(H)=1/2.", "finite_counts": counts},
        "decision": {"bounded_result": "P3152 constructs the missing charge-unit object U_Y^tors and proves a finite scale-torsor obstruction for strict normalization.", "why_not_strict": "The only current row that fixes Y(H)=1/2 imports the installed SM charge convention.  Existing unit/action and dimension-source audits do not export a noncircular U(1) charge valuation from strict nadsoliton data.", "next_honest_step": "Pivot away from SM charge normalization unless a new strict charge-valuation object is introduced.  The next proof-grade route is P3153: a GR/EH nonproxy coupling audit for the axiom branch, testing one explicit metric/EH source interface without claiming L_total, SM generation, selector closure, bridge completion, role transfer, or ToE.", "negative_export_flags": {"absolute_charge_normalization_exported": False, "strict_charge_valuation_exported": False, "strict_SM_generation_exported": False, "unit_bearing_L_total_exported": False, "global_BV_BRST_exported": False, "strict_GR_generation_exported": False, "selector_closure_exported": False, "bridge_or_role_transfer_exported": False, "ToE_closure_exported": False}},
    }


def write_markdown(payload: dict[str, Any]) -> None:
    th = payload["finite_theorem"]
    lines = ["# P3152/S2102 Y_SM charge-unit normalization obstruction", "", f"Status: `{payload['status']}`", "", "## Constructed object", f"- `{payload['constructed_object']['name']}`", f"- Classification: `{payload['constructed_object']['classification']}`", f"- Scope: `{payload['constructed_object']['scope']}`", "", "## Finite theorem", f"`{th['name']}`: {th['statement']}", "", "## Finite counts"]
    for k, v in th["finite_counts"].items():
        lines.append(f"- `{k}`: `{v}`")
    lines.extend(["", "## Scale orbit rows"])
    for row in payload["scale_orbit_rows"]:
        lines.append(f"- `c={row['scale_c']}`: `Y(H)={row['Y_H']}`, residuals zero `{row['all_yukawa_and_anomaly_residuals_zero']}`, hypercharge gcd `{row['hypercharge_unit_gcd']}`, electric-charge gcd `{row['electric_charge_gcd']}`")
    lines.extend(["", "## Candidate source rows"])
    for row in payload["candidate_source_rows"]:
        lines.append(f"- `{row['candidate']}`: fixes scale `{row['fixes_scale']}`, strict source `{row['strict_source_law']}`, noncircular `{row['noncircular']}`; {row['verdict']}")
    lines.extend(["", "## Decision", payload["decision"]["bounded_result"], "", "## Why this is not strict", payload["decision"]["why_not_strict"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3152/S2102 Y_SM charge-unit normalization obstruction", "## P3152/S2102 Y_SM charge-unit normalization obstruction\n\n`P3152/S2102` constructs `U_Y^tors`, the nonzero-scale torsor for the P3150 hypercharge ray.  The sampled rational orbit has zero Yukawa/anomaly residuals at every scale while `Y(H)`, hypercharge gcd, and electric-charge values vary.  Thus current algebraic compatibility fixes only a projective ray and does not export a strict charge unit, strict SM generation, unit-bearing `L_total`, global BV/BRST, GR/EH closure, selector closure, bridge/role transfer, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3152/S2102 charge unit remains unsourced", "## P3152/S2102 charge unit remains unsourced\n\n`P3152/S2102` proves that the local matter/Higgs algebra remains invariant under nonzero hypercharge rescaling.  Fixing `Y(H)=1/2` still imports the SM electric-charge convention; a unit-bearing Lagrangian/EOM would require a noncircular strict U(1) charge valuation or a different physical-unit source theorem.\n")
    append_once(AGENTS, "Current Y_SM charge-unit normalization guardrail (P3152/S2102, 2026-07-13)", "## Current Y_SM charge-unit normalization guardrail (P3152/S2102, 2026-07-13)\n\n- P3152 constructs `U_Y^tors`, the charge-unit normalization torsor for the P3150 hypercharge ray.\n- All sampled nonzero rational rescalings preserve the Yukawa/anomaly residual equations, while `Y(H)`, the hypercharge gcd, and electric-charge values vary.\n- This closes the current axiom-branch SM source-selection stack as conditional: representation content, hypercharge ray, and charge unit are not strict-sourced by current artifacts.\n- Do not promote `U_Y^tors`, `Y_SM^ray`, `R_shape^scan`, or `R_SM^ax` to strict SM generation, unit-bearing `L_total`, global BV/BRST, GR/EH closure, selector closure, bridge/role transfer, or ToE.\n- Next honest move: pivot to a GR/EH nonproxy coupling audit for the axiom branch, unless a genuinely new strict U(1) charge-valuation object is introduced.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
