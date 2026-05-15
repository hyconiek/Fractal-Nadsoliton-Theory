#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1694 = GEN / "p1694_s644_strict_kernel_to_full_lagrangian_bidirectional_map_witness.json"
IN1662 = GEN / "p1662_s612_strict_full_lagrangian_explicit_density_summary.json"
IN1700 = GEN / "p1700_s650_strict_eom_variational_residual_identity_checkpoint.json"
IN1705 = GEN / "p1705_s655_strict_nonproxy_metric_spinor_variational_objects_export_checkpoint.json"
OUT = GEN / "p1706_s656_strict_full_chain_kernel_to_full_lagrangian_to_eom_dossier_checkpoint.json"


def main() -> None:
    p1694 = json.loads(IN1694.read_text(encoding="utf-8"))
    p1662 = json.loads(IN1662.read_text(encoding="utf-8"))
    p1700 = json.loads(IN1700.read_text(encoding="utf-8"))
    p1705 = json.loads(IN1705.read_text(encoding="utf-8"))

    payload = {
        "checkpoint": "P1706_S656",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "full_chain": {
            "kernel": p1694.get("kernel"),
            "kernel_input": p1694.get("kernel_input", {}),
            "coefficient_map": p1694.get("forward_coefficient_map_symbolic", {}),
            "full_lagrangian_explicit": p1662.get("full_lagrangian_density_explicit", {}),
            "reduced_eom_bundle_anchor": p1700.get("EOM_bundle", {}),
            "reverse_variational_identity_anchor": {
                "status": p1700.get("identity_status", "UNKNOWN"),
                "residuals": p1700.get("variational_identity_residuals", {}),
            },
            "nonproxy_variational_objects_anchor": {
                "metric": p1705.get("metric_variational_objects_export", {}),
                "spinor": p1705.get("spinor_variational_objects_export", {}),
            },
        },
        "bidirectional_readout": {
            "forward_kernel_to_coefficients": "EXPORTED",
            "forward_coefficients_to_full_lagrangian": "EXPORTED",
            "forward_lagrangian_to_eom": "EXPORTED_REDUCED_PLUS_NONPROXY_OBJECT_SCAFFOLD",
            "reverse_eom_to_variational_origin_local": p1700.get("identity_status", "UNKNOWN"),
            "reverse_global_nonproxy_theorem": "OPEN",
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": [
            "full_nonproxy_metric_tensor_eom_export",
            "full_nonproxy_spinor_connection_eom_export",
            "global_helmholtz_integrability_nonproxy",
            "brst_nilpotency_and_cohomology_nonproxy",
            "cutkosky_unitarity_full_sector",
            "counterterm_flow_renormalization_closure",
            "background_independence_family_theorem",
            "qw2191_selector_source_or_nonclosure_theorem",
        ],
        "next_honest_step": "Wyprowadzić jawny full nonproxy metric+spinor covariant EOM export (nie tylko obiekty kontraktowe), aby otworzyć formalny etap global_helmholtz_integrability_nonproxy i BRST nilpotency proof.",
        "lay_summary": "To pełny raport: od kernela strict, przez współczynniki i pełny lagranżian, aż po równania ruchu i kontrolę wsteczną. Nadal brakuje finalnych, globalnych dowodów nonproxy wymaganych do domknięcia teorii.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
