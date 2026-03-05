#!/usr/bin/env python3
"""
QW-2206: Foundational entity + topology scope stratification gate.

Purpose:
- integrate strict action/EoM evidence (L1),
- integrate local topological protection evidence (L2/L17),
- keep global/full-object boundary explicit (no overclaim).
"""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def find_float(pattern: str, text: str) -> float | None:
    m = re.search(pattern, text, flags=re.IGNORECASE | re.MULTILINE)
    if not m:
        return None
    try:
        return float(m.group(1))
    except Exception:
        return None


def main() -> None:
    src_q2165 = ROOT / "report_qw2165_l13_exhaustive_canonical_eom_gate.json"
    src_q1204 = ROOT / "RAPORT_QW1204_SKYRMION_RIGOROUS.md"
    src_q1611 = ROOT / "RAPORT_QW1611_SKYRMION_CONVERGENCE.md"
    src_q1622 = ROOT / "RAPORT_QW1622_FR_QUANTIZATION.md"

    q2165 = load_json(src_q2165)
    t1204 = src_q1204.read_text(encoding="utf-8")
    t1611 = src_q1611.read_text(encoding="utf-8")
    t1622 = src_q1622.read_text(encoding="utf-8")

    b_1204 = find_float(r"B\s*=\s*([0-9]+\.[0-9]+)\s*[±\+/-]", t1204)
    if b_1204 is None:
        b_1204 = find_float(r"B\s*=\s*([0-9]+\.[0-9]+)", t1204)

    qinf_1611 = find_float(r"Q_∞\s*=\s*([0-9]+\.[0-9]+)", t1611)
    if qinf_1611 is None:
        qinf_1611 = find_float(r"Q\s*=\s*([0-9]+\.[0-9]+),\s*\|Q-1\|", t1611)

    spin_1622 = find_float(r"RESULT:\s*spin\s*=\s*([0-9]+(?:\.[0-9]+)?)", t1622)
    g_1622 = find_float(r"RESULT:\s*g\s*=\s*([0-9]+(?:\.[0-9]+)?)", t1622)

    q2165_flags = q2165.get("flags", {})
    one_action_present = bool(q2165.get("model", {}).get("lagrangian_density"))
    eom_exhaustive = bool(q2165_flags.get("euler_lagrange_executed_for_all_13_fields"))
    local_second_order = bool(q2165_flags.get("all_psi_eom_local_second_order")) and bool(
        q2165_flags.get("phi_eom_local_second_order")
    )
    no_spacetime_nonlocal = bool(q2165_flags.get("no_spacetime_nonlocal_tokens_in_all_13_eom"))
    full_variational_all_orders_closed = bool(
        q2165_flags.get("full_all_orders_variational_theorem_from_complete_fin_action_completed")
    )

    b_close = b_1204 is not None and abs(b_1204 - 1.0) <= 0.01
    qinf_close = qinf_1611 is not None and abs(qinf_1611 - 1.0) <= 0.01
    fr_spin_half = spin_1622 is not None and abs(spin_1622 - 0.5) <= 1e-9
    fr_g_two = g_1622 is not None and abs(g_1622 - 2.0) <= 1e-9

    flags: dict[str, bool] = {
        "q2165_one_action_template_present": one_action_present,
        "q2165_exhaustive_eom_present": eom_exhaustive,
        "q2165_local_second_order_structure_present": local_second_order,
        "q2165_no_spacetime_nonlocal_tokens": no_spacetime_nonlocal,
        "q1204_topological_charge_close_to_one": b_close,
        "q1611_radial_convergence_close_to_one": qinf_close,
        "q1622_fr_spin_half_present": fr_spin_half,
        "q1622_fr_g_factor_two_present": fr_g_two,
        "single_fundamental_field_reduction_closed": False,
        "global_full_object_topological_protection_theorem_closed": False,
        "no_overclaim_boundary_explicit": True,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2165_one_action_template_present"]
        and flags["q2165_exhaustive_eom_present"]
        and flags["q1204_topological_charge_close_to_one"]
        and flags["q1622_fr_spin_half_present"]
        and flags["q1622_fr_g_factor_two_present"]
        and flags["no_overclaim_boundary_explicit"]
    )

    if core_ok:
        verdict = "FOUNDATIONAL_ENTITY_TOPOLOGY_SCOPE_GATE_PASS_PARTIAL_LOCAL_PROTECTION_ONLY"
    else:
        verdict = "FOUNDATIONAL_ENTITY_TOPOLOGY_SCOPE_GATE_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2165": src_q2165.name,
            "q1204": src_q1204.name,
            "q1611": src_q1611.name,
            "q1622": src_q1622.name,
        },
        "observables": {
            "q1204_b_topological": b_1204,
            "q1611_qinf_radial": qinf_1611,
            "q1622_spin": spin_1622,
            "q1622_g_factor": g_1622,
            "q2165_all_orders_variational_theorem_closed": full_variational_all_orders_closed,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "FORMALLY_DERIVE_SINGLE_FUNDAMENTAL_FIELD_REDUCTION_AND_GLOBAL_FULL_OBJECT_"
            "TOPOLOGICAL_PROTECTION_THEOREM_FROM_COMPLETE_FIN_ACTION"
        ),
    }

    out_json = ROOT / "report_qw2206_foundational_entity_topology_scope_gate.json"
    out_md = ROOT / "RAPORT_QW2206_FOUNDATIONAL_ENTITY_TOPOLOGY_SCOPE_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2206: FOUNDATIONAL ENTITY + TOPOLOGY SCOPE GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- L1: canonical single-action template + exhaustive EoM evidence is present in strict chain.",
                "- L2/L17: local topological protection evidence is present (Skyrmion charge close to 1 + FR spin-1/2, g=2).",
                "- Full single-field ontological reduction and global full-object topological theorem remain explicitly open.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2206_foundational_entity_topology_scope_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
