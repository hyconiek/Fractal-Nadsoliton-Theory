#!/usr/bin/env python3
"""
QW-2074: Strict no-fit derivation round for missing SM+GR parameters.

Principle:
- no optimization, no retune, no threshold search,
- explicit dependency labeling to avoid false "first-principles" claims.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2074_strict_nofit_missing_parameter_derivations.json"
OUT_MD = ROOT / "RAPORT_QW2074_STRICT_NOFIT_MISSING_PARAMETER_DERIVATIONS.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def get_ref(groups: Dict, param_id: str):
    for _, items in groups.items():
        for item in items:
            if item["id"] == param_id:
                return item["value"]
    raise KeyError(f"Missing parameter in registry: {param_id}")


def solve_mw_from_alpha_gf_mz(alpha: float, g_f: float, m_z: float, delta_r: float = 0.0) -> float:
    # On-shell relation:
    # M_W^2 (1 - M_W^2/M_Z^2) = (pi alpha)/(sqrt(2) G_F) * 1/(1-delta_r)
    a = (math.pi * alpha) / (math.sqrt(2.0) * g_f) * (1.0 / (1.0 - delta_r))
    mz2 = m_z * m_z
    disc = 1.0 - 4.0 * a / mz2
    if disc < 0.0:
        raise ValueError("No real MW solution in on-shell equation.")
    x = 0.5 * mz2 * (1.0 + math.sqrt(disc))
    return float(math.sqrt(max(x, 0.0)))


def main() -> None:
    reg = load_json("report_qw2068_sm_gr_parameter_registry.json")
    groups = reg["groups"]

    alpha_inv = float(get_ref(groups, "alpha_em_inv_mz"))
    alpha = 1.0 / alpha_inv
    g_f = float(get_ref(groups, "g_f"))
    m_z = float(get_ref(groups, "m_z"))

    # No-fit anchor-dependent physical relations.
    v_from_gf = (math.sqrt(2.0) * g_f) ** (-0.5)
    m_w_from_tree_onshell = solve_mw_from_alpha_gf_mz(alpha=alpha, g_f=g_f, m_z=m_z, delta_r=0.0)

    # SI-definitional constants (not dynamic derivations).
    c_light = 299792458.0
    h_planck_exact = 6.62607015e-34
    hbar_from_h = h_planck_exact / (2.0 * math.pi)

    updates: List[Dict] = [
        {
            "id": "v_higgs",
            "predicted_value": float(v_from_gf),
            "status": "derived_nofit_anchor_dependent",
            "strict_level": "physical_relation_anchor_dependent",
            "method": "v=(sqrt(2)G_F)^(-1/2)",
            "notes": "Depends on G_F anchor; not standalone first-principles closure.",
        },
        {
            "id": "m_w",
            "predicted_value": float(m_w_from_tree_onshell),
            "status": "derived_nofit_anchor_dependent",
            "strict_level": "physical_relation_anchor_dependent",
            "method": "on_shell_relation(alpha,G_F,M_Z,delta_r=0)",
            "notes": "Tree-level/on-shell baseline; full EW radiative precision still separate.",
        },
        {
            "id": "c_light",
            "predicted_value": float(c_light),
            "status": "definition_constant",
            "strict_level": "si_definition",
            "method": "SI exact definition",
            "notes": "Metrological definition constant, not dynamic field derivation.",
        },
        {
            "id": "hbar",
            "predicted_value": float(hbar_from_h),
            "status": "definition_constant",
            "strict_level": "si_definition",
            "method": "hbar=h/(2pi), h exact SI",
            "notes": "Metrological definition chain, not dynamic field derivation.",
        },
    ]

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {"registry": "report_qw2068_sm_gr_parameter_registry.json"},
        "updates": updates,
        "count_updates": len(updates),
        "verdict": "STRICT_NOFIT_MISSING_PARAMETER_DERIVATION_ROUND1",
        "required_next_step": "INTEGRATE_UPDATES_IN_QW2069_WITH_EPISTEMIC_LABELS",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2074: STRICT NO-FIT MISSING PARAMETER DERIVATIONS",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- updates: {len(updates)}",
        "",
        "## Key Values",
        f"- v_higgs_from_G_F: {v_from_gf:.6f} GeV",
        f"- m_w_on_shell_tree: {m_w_from_tree_onshell:.6f} GeV",
        f"- c_light (SI): {c_light:.1f} m/s",
        f"- hbar_from_h: {hbar_from_h:.12e} J*s",
        "",
        "## Epistemic Rule",
        "- These updates are explicitly labeled as anchor-dependent or definitional.",
        "- They are not promoted to strict first-principles closure by themselves.",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2074] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2074] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2074] verdict={out['verdict']} updates={len(updates)}")


if __name__ == "__main__":
    main()
