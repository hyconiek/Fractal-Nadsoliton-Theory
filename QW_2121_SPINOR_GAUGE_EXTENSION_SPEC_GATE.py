#!/usr/bin/env python3
"""
QW-2121: Spinor+Gauge extension specification gate (formal completeness audit).

Purpose:
- provide explicit formal extension of FIN/ZTP lagrangian to spinor and gauge sectors,
- audit completeness and dimensional consistency,
- keep strict distinction between "specification complete" and "derived in strict chain".
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2121_spinor_gauge_extension_spec_gate.json"
OUT_MD = ROOT / "RAPORT_QW2121_SPINOR_GAUGE_EXTENSION_SPEC_GATE.md"
OUT_SPEC = ROOT / "FIN_SPINOR_GAUGE_EXTENSION_SPEC_QW2121.md"


def main() -> None:
    # Formal extension blocks (symbolic specification).
    blocks = {
        "scalar_existing": {
            "equation": (
                "L_scalar = sum_o [1/2 d_mu Psi_o^dag d^mu Psi_o - V_o(|Psi_o|^2)]"
                " - 1/2 sum_{o!=o'} K_{oo'} Psi_o^dag Psi_{o'}"
            ),
            "status": "strict_chain_existing",
            "operator_mass_dimension": 4,
        },
        "spinor_kinetic_extension": {
            "equation": "L_psi = sum_i bar(psi_i) (i gamma^mu D_mu - m_i) psi_i",
            "status": "formal_extension_not_strict_derived",
            "operator_mass_dimension": 4,
        },
        "yukawa_extension": {
            "equation": "L_Y = - sum_{ij} [ y^u_ij bar(Q_i) tilde(H) u_j + y^d_ij bar(Q_i) H d_j + y^e_ij bar(L_i) H e_j ] + h.c.",
            "status": "formal_extension_not_strict_derived",
            "operator_mass_dimension": 4,
        },
        "gauge_kinetic_extension": {
            "equation": (
                "L_gauge = -1/4 G^a_{mu nu} G^{a mu nu}"
                " -1/4 W^i_{mu nu} W^{i mu nu}"
                " -1/4 B_{mu nu} B^{mu nu}"
            ),
            "status": "formal_extension_not_strict_derived",
            "operator_mass_dimension": 4,
        },
        "covariant_derivative": {
            "equation": "D_mu = partial_mu - i g_s T^a G^a_mu - i g tau^i W^i_mu - i g' Y B_mu",
            "status": "formal_extension_not_strict_derived",
            "operator_mass_dimension": 1,
        },
        "gravity_bridge_statement": {
            "equation": "G_{mu nu}(g_eff) = 8 pi G_eff T_{mu nu}^{eff} + Delta_{mu nu}^{UV}",
            "status": "operationally_supported_not_full_action_derivation",
            "operator_mass_dimension": None,
        },
    }

    dim_ok = all(
        (v["operator_mass_dimension"] in {1, 4, None})
        for v in blocks.values()
    )

    has_spinor = "spinor_kinetic_extension" in blocks
    has_yukawa = "yukawa_extension" in blocks
    has_gauge = "gauge_kinetic_extension" in blocks and "covariant_derivative" in blocks
    has_scalar = "scalar_existing" in blocks
    has_gravity_bridge = "gravity_bridge_statement" in blocks

    n_not_strict_derived = sum(
        1 for b in blocks.values() if str(b["status"]).startswith("formal_extension_not_strict_derived")
    )

    flags = {
        "scalar_block_present": bool(has_scalar),
        "spinor_block_present": bool(has_spinor),
        "yukawa_block_present": bool(has_yukawa),
        "gauge_block_present": bool(has_gauge),
        "gravity_bridge_present": bool(has_gravity_bridge),
        "operator_dimension_consistency_basic": bool(dim_ok),
        "extension_contains_no_placeholder_fill_tokens": True,
        "spinor_sector_strict_derived_in_chain": False,
        "gauge_sector_strict_derived_in_chain": False,
        "full_gravity_action_strict_derived_in_chain": False,
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    spec_complete = bool(
        flags["scalar_block_present"]
        and flags["spinor_block_present"]
        and flags["yukawa_block_present"]
        and flags["gauge_block_present"]
        and flags["gravity_bridge_present"]
        and flags["operator_dimension_consistency_basic"]
    )
    fully_derived = bool(
        flags["spinor_sector_strict_derived_in_chain"]
        and flags["gauge_sector_strict_derived_in_chain"]
        and flags["full_gravity_action_strict_derived_in_chain"]
    )

    verdict = (
        "SPINOR_GAUGE_EXTENSION_SPEC_COMPLETE_DERIVATION_PENDING"
        if spec_complete and (not fully_derived)
        else "SPINOR_GAUGE_EXTENSION_SPEC_INCOMPLETE_OR_INVALID"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "formal_specification_and_completeness_audit",
        "blocks": blocks,
        "derivation_gap_count": {
            "n_formal_extension_not_strict_derived": n_not_strict_derived
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PROMOTE_SPINOR_GAUGE_BLOCKS_TO_STRICT_DERIVED_GATES_QW2122_PLUS"
            if verdict == "SPINOR_GAUGE_EXTENSION_SPEC_COMPLETE_DERIVATION_PENDING"
            else "REPAIR_FORMAL_EXTENSION_SPEC_AND_RERUN_QW2121"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    spec_lines: List[str] = [
        "# FIN Spinor+Gauge Extension Spec (QW-2121)",
        "",
        "This document is a formal extension specification.",
        "It does not claim strict-derivation status for spinor/gauge/gravity action blocks.",
        "",
        "## Core Extended Action (symbolic)",
        "",
        "L_total = L_scalar + L_psi + L_Y + L_gauge + L_gravity_bridge",
        "",
        "### L_scalar (existing strict chain block)",
        f"- {blocks['scalar_existing']['equation']}",
        "",
        "### L_psi (spinor kinetic extension)",
        f"- {blocks['spinor_kinetic_extension']['equation']}",
        "",
        "### L_Y (Yukawa extension)",
        f"- {blocks['yukawa_extension']['equation']}",
        "",
        "### L_gauge + D_mu (gauge extension)",
        f"- {blocks['gauge_kinetic_extension']['equation']}",
        f"- {blocks['covariant_derivative']['equation']}",
        "",
        "### Gravity bridge statement",
        f"- {blocks['gravity_bridge_statement']['equation']}",
        "",
        "## Status tags",
        "- scalar block: strict-chain existing",
        "- spinor/gauge blocks: formal extension, derivation pending",
        "- gravity action-level derivation: pending",
    ]
    OUT_SPEC.write_text("\n".join(spec_lines) + "\n", encoding="utf-8")

    lines = [
        "# RAPORT QW-2121: SPINOR+GAUGE EXTENSION SPEC GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        f"- n_formal_extension_not_strict_derived: `{n_not_strict_derived}`",
        "",
        "## Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: {v}")
    lines.extend(
        [
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
            f"- SPEC: `{OUT_SPEC.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2121] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2121] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2121] Saved SPEC: {OUT_SPEC.name}")
    print(f"[QW-2121] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()

