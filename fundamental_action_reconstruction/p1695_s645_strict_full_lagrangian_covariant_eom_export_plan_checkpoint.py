#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1691 = GEN / "p1691_s641_strict_full_chain_lagrangian_to_qg_theorem_obligation_matrix.json"
IN1694 = GEN / "p1694_s644_strict_kernel_to_full_lagrangian_bidirectional_map_witness.json"
OUT = GEN / "p1695_s645_strict_full_lagrangian_covariant_eom_export_plan_checkpoint.json"


def main() -> None:
    p1691 = json.loads(IN1691.read_text(encoding="utf-8"))
    p1694 = json.loads(IN1694.read_text(encoding="utf-8"))

    full_L = p1691.get("full_lagrangian_density", {})

    covariant_eom_export_matrix = {
        "metric_g": {"lagrangian_block": "L_gravity + L_mix", "eom_status": "OPEN_COVARIANT_EXPORT_REQUIRED"},
        "gauge_Aa": {"lagrangian_block": "L_gauge + L_mix", "eom_status": "OPEN_COVARIANT_EXPORT_REQUIRED"},
        "higgs_H": {"lagrangian_block": "L_higgs + L_mix + L_yukawa", "eom_status": "OPEN_COVARIANT_EXPORT_REQUIRED"},
        "fermion_psif": {"lagrangian_block": "L_fermion + L_yukawa", "eom_status": "OPEN_COVARIANT_EXPORT_REQUIRED"},
        "scalar_phi": {"lagrangian_block": "L_scalar_phi + L_mix", "eom_status": "OPEN_COVARIANT_EXPORT_REQUIRED"},
    }

    theorem_bridge = {
        "counterterm_flow_closure": "KEEP_OPEN",
        "brst_cohomology_closure": "KEEP_OPEN",
        "cutkosky_full_sector_unitarity": "KEEP_OPEN",
        "background_family_independence": "KEEP_OPEN",
    }

    payload = {
        "checkpoint": "P1695_S645",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "forward_chain_anchor": p1694.get("full_chain"),
        "kernel": p1694.get("kernel"),
        "kernel_input": p1694.get("kernel_input"),
        "coefficient_anchor": p1694.get("forward_coefficient_map_symbolic"),
        "full_lagrangian_density": full_L,
        "covariant_eom_export_matrix": covariant_eom_export_matrix,
        "bidirectional_note": "Forward map exported; reverse global EOM-family to variational-origin theorem remains open.",
        "qg_theorem_obligations": theorem_bridge,
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "next_honest_step": "Wyeksportować pełne kowariantne EOM dla każdego bloku L_total i dołączyć wspólny theorem witness łączący counterterm-flow + BRST/Cutkosky + background-family.",
        "lay_summary": "Mamy pełny wzór teorii i plan dokładnego wyprowadzenia równań ruchu dla każdego sektora. Nadal brakuje końcowych, globalnych dowodów dla kwantowej grawitacji.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
