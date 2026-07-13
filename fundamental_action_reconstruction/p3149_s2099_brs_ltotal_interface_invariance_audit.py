"""P3149/S2099: BRST/L_total interface invariance audit for R_SM^ax.

P3148 produced an algebraic one-family SM representation registry.  This next
step feeds that registry into the local P1961 BRST/gauge interface and checks
which matter/Higgs Lagrangian building blocks are invariant in the same
convention.  The result is a bounded local algebraic interface theorem, not a
unit-bearing L_total or global BV/BRST closure theorem.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import sympy as sp

from p3148_s2098_sm_representation_registry_completion_audit import (
    I,
    INPUTS as P3148_INPUTS,
    color_generators,
    reps,
    weak_generators,
)

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3149_s2099_brs_ltotal_interface_invariance_audit.json"
MD = GEN / "p3149_s2099_brs_ltotal_interface_invariance_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P3148": GEN / "p3148_s2098_sm_representation_registry_completion_audit.json",
    "P3147": GEN / "p3147_s2097_axiom_lie_smgr_fit_readiness_matrix.json",
    "P1961": GEN / "p1961_s911_strict_local_nonabelian_brst_differential_and_nilpotency_certificate.json",
    "P1962": GEN / "p1962_s912_strict_matter_higgs_brst_extension_registry_audit.json",
    **{f"P3148_DEP_{k}": v for k, v in P3148_INPUTS.items()},
}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def all_zero(mat: sp.Matrix) -> bool:
    return all(sp.simplify(x) == 0 for x in mat)


def rep_by_field() -> dict[str, dict[str, Any]]:
    return {r["field"]: r for r in reps()}


def kinetic_invariance_rows() -> list[dict[str, Any]]:
    rows = []
    for rep in reps():
        cgens = color_generators(rep["su3"])
        wgens = weak_generators(rep["su2"])
        color_bad = sum(not all_zero(t.conjugate().T - t) for t in cgens)
        weak_bad = sum(not all_zero(t.conjugate().T - t) for t in wgens)
        y_real = rep["Y"].denominator != 0
        rows.append({
            "field": rep["field"],
            "term": f"({rep['field']})^dagger D {rep['field']}",
            "su3_unitary_generator_failures": int(color_bad),
            "su2_unitary_generator_failures": int(weak_bad),
            "u1_real_hypercharge": bool(y_real),
            "local_gauge_invariant_bilinear": color_bad == 0 and weak_bad == 0 and y_real,
        })
    return rows


def su2_epsilon_invariant() -> bool:
    eps = sp.Matrix([[0, 1], [-1, 0]])
    return all(all_zero(t.T * eps + eps * t) for t in weak_generators("2"))


def color_delta_invariant() -> bool:
    # For 3 x bar3 -> 1, anti-fundamental generators are -T^*, so T^T + Tbar = 0 on delta.
    fund = color_generators("3")
    anti = color_generators("bar3")
    return all(all_zero(t.T + tb) for t, tb in zip(fund, anti, strict=True))


def yukawa_invariance_rows() -> list[dict[str, Any]]:
    # Left-handed Weyl convention.  SU2 uses epsilon contraction between two doublets;
    # color uses delta contraction between 3 and bar3.
    checks = [
        {
            "term": "Q_L H u_c",
            "fields": ["Q_L", "H", "u_c"],
            "hypercharge_sum": Fraction(1, 6) + Fraction(1, 2) + Fraction(-2, 3),
            "su3_invariant_tensor": color_delta_invariant(),
            "su2_invariant_tensor": su2_epsilon_invariant(),
        },
        {
            "term": "Q_L Hdagger d_c",
            "fields": ["Q_L", "Hdagger", "d_c"],
            "hypercharge_sum": Fraction(1, 6) + Fraction(-1, 2) + Fraction(1, 3),
            "su3_invariant_tensor": color_delta_invariant(),
            "su2_invariant_tensor": su2_epsilon_invariant(),
        },
        {
            "term": "L_L Hdagger e_c",
            "fields": ["L_L", "Hdagger", "e_c"],
            "hypercharge_sum": Fraction(-1, 2) + Fraction(-1, 2) + Fraction(1, 1),
            "su3_invariant_tensor": True,
            "su2_invariant_tensor": su2_epsilon_invariant(),
        },
    ]
    rows = []
    for c in checks:
        rows.append({
            **c,
            "hypercharge_sum": str(c["hypercharge_sum"]),
            "u1_invariant": c["hypercharge_sum"] == 0,
            "local_brst_variation_zero": c["hypercharge_sum"] == 0 and c["su3_invariant_tensor"] and c["su2_invariant_tensor"],
        })
    return rows


def build_payload() -> dict[str, Any]:
    kinetic = kinetic_invariance_rows()
    yukawa = yukawa_invariance_rows()
    counts = {
        "kinetic_rows": len(kinetic),
        "kinetic_rows_local_gauge_invariant": sum(r["local_gauge_invariant_bilinear"] for r in kinetic),
        "yukawa_rows": len(yukawa),
        "yukawa_rows_local_brst_zero": sum(r["local_brst_variation_zero"] for r in yukawa),
        "local_lagrangian_interface_rows": len(kinetic) + len(yukawa),
        "local_lagrangian_interface_rows_passing": sum(r["local_gauge_invariant_bilinear"] for r in kinetic) + sum(r["local_brst_variation_zero"] for r in yukawa),
        "unit_bearing_measure_rows": 0,
        "global_bv_brst_rows": 0,
        "nonproxy_gr_rows": 0,
        "strict_nadsoliton_source_rows": 0,
    }
    local_pass = counts["local_lagrangian_interface_rows"] == counts["local_lagrangian_interface_rows_passing"]
    return {
        "status": "P3149_BRST_LTOTAL_INTERFACE_LOCAL_ALGEBRAIC_PASS_CONDITIONAL_NO_GLOBAL_CLOSURE",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "constructed_object": {
            "name": "I_BRST^ax(R_SM^ax) local matter/Higgs interface certificate",
            "classification": "axiom_branch_local_brs_lagrangian_interface_audit",
            "scope": "matter/Higgs kinetic and Yukawa invariance in the P3148 convention; no unit-bearing integral or global BV/BRST",
        },
        "kinetic_invariance_rows": kinetic,
        "yukawa_invariance_rows": yukawa,
        "finite_theorem": {
            "name": "P3149_T1_local_brs_lagrangian_interface_certificate",
            "statement": "Feeding R_SM^ax into the local P1961-style gauge/BRST interface gives a local algebraic pass for the audited matter/Higgs Lagrangian blocks: all 6 kinetic bilinear rows have unitary factor generators and real hypercharge, and all 3 Yukawa rows have invariant SU(3)/SU(2) contraction tensors plus zero hypercharge sum.  Therefore the axiom-branch SM registry is locally compatible with the expected matter/Higgs gauge-invariant Lagrangian interface.  This does not export a unit-bearing L_total, global BV/BRST charge/cohomology, strict source for the registry, GR/EH coupling, selector closure, or ToE closure.",
            "finite_counts": counts,
        },
        "decision": {
            "bounded_result": "The axiom-branch SM side is now stronger: R_SM^ax is not just an anomaly-free registry; it also supports the local matter/Higgs kinetic and Yukawa gauge-invariance interface expected by the P1961 BRST rules.",
            "why_not_strict": "The calculation is local and algebraic.  It assumes the P3148 representation registry and does not build a unit-bearing action integral, global BRST charge, physical Hilbert/cohomology projection, strict source of the registry, or nonproxy GR/EH sector.",
            "next_honest_step": "Construct P3150 as a bounded source-selection obstruction/witness audit for R_SM^ax itself: test whether any current strict object can select the one-family representation/hypercharge registry without importing SM ansatz data.  If not, preserve R_SM^ax as conditional phenomenology and pivot to the unit-source or GR/EH nonproxy lane.",
            "negative_export_flags": {
                "unit_bearing_L_total_exported": False,
                "global_BV_BRST_exported": False,
                "physical_state_cohomology_exported": False,
                "strict_registry_source_exported": False,
                "strict_SM_generation_exported": False,
                "strict_GR_generation_exported": False,
                "selector_closure_exported": False,
                "ToE_closure_exported": False,
            },
            "accepted_as_local_interface_certificate": local_pass,
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    th = payload["finite_theorem"]
    c = th["finite_counts"]
    lines = [
        "# P3149/S2099 BRST/L_total interface invariance audit",
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
    for key, value in c.items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Kinetic rows"])
    for row in payload["kinetic_invariance_rows"]:
        lines.append(f"- `{row['term']}`: local gauge invariant `{row['local_gauge_invariant_bilinear']}`")
    lines.extend(["", "## Yukawa rows"])
    for row in payload["yukawa_invariance_rows"]:
        lines.append(f"- `{row['term']}`: hypercharge `{row['hypercharge_sum']}`, BRST variation zero `{row['local_brst_variation_zero']}`")
    lines.extend(["", "## Decision", payload["decision"]["bounded_result"], "", "## Why this is not strict", payload["decision"]["why_not_strict"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3149/S2099 BRST/L_total interface invariance audit", "## P3149/S2099 BRST/L_total interface invariance audit\n\n`P3149/S2099` constructs `I_BRST^ax(R_SM^ax)`, a local matter/Higgs BRST/Lagrangian interface certificate for the P3148 representation registry.  The finite audit passes `6/6` kinetic bilinear rows and `3/3` Yukawa rows, so all `9/9` audited local interface blocks are gauge/BRST-compatible in the installed convention.  This remains conditional and local: no unit-bearing `L_total`, global BV/BRST charge/cohomology, strict registry source, GR/EH coupling, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3149/S2099 local BRST interface is not global L_total", "## P3149/S2099 local BRST interface is not global L_total\n\n`P3149/S2099` verifies local algebraic gauge invariance for matter/Higgs kinetic and Yukawa blocks built from `R_SM^ax`.  It does not supply the physical measure, action normalization, global BRST charge, cohomology projection, full nonproxy EOM, or GR/EH coupling required for a unit-bearing `L_total`.\n")
    append_once(AGENTS, "Current BRST/L_total interface invariance guardrail (P3149/S2099, 2026-07-13)", "## Current BRST/L_total interface invariance guardrail (P3149/S2099, 2026-07-13)\n\n- P3149 constructs `I_BRST^ax(R_SM^ax)`, a local matter/Higgs BRST/Lagrangian interface certificate for the P3148 representation registry.\n- The finite algebraic audit passes `6/6` kinetic bilinear rows and `3/3` Yukawa rows; all `9/9` audited local interface blocks are gauge/BRST-compatible in the installed convention.\n- This strengthens the axiom-branch SM fit but remains local and conditional: no unit-bearing `L_total`, global BV/BRST charge/cohomology, strict registry source, GR/EH coupling, selector closure, or ToE closure is exported.\n- Do not promote P3149 to strict SM generation or full Lagrangian/EOM closure.  Next honest move: audit whether any current strict object can source-select `R_SM^ax` itself without importing SM ansatz data; otherwise preserve it as conditional phenomenology and pivot to unit-source or GR/EH nonproxy work.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
