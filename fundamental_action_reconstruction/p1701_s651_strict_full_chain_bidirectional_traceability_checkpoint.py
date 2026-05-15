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
OUT = GEN / "p1701_s651_strict_full_chain_bidirectional_traceability_checkpoint.json"


def main() -> None:
    p1694 = json.loads(IN1694.read_text(encoding="utf-8"))
    p1662 = json.loads(IN1662.read_text(encoding="utf-8"))
    p1700 = json.loads(IN1700.read_text(encoding="utf-8"))

    full_density = p1662.get("full_lagrangian_density_explicit", {})
    ordered_keys = ["L_scalar", "L_gauge", "L_fermion", "L_higgs", "L_gravity", "L_mix"]
    full_lagrangian_explicit_ordered = {k: full_density.get(k, "MISSING") for k in ordered_keys}

    payload = {
        "checkpoint": "P1701_S651",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": {
            "forward": "K_strict -> coefficients -> full explicit L_total -> sector EOM bundle",
            "reverse": "EOM bundle -> EL residual identity (reduced proxy pass) -> global theorem-level inversion OPEN",
        },
        "kernel_anchor": {
            "kernel": p1694.get("kernel"),
            "kernel_input": p1694.get("kernel_input", {}),
        },
        "coefficient_anchor": p1694.get("forward_coefficient_map_symbolic", {}),
        "full_lagrangian_explicit": full_lagrangian_explicit_ordered,
        "eom_anchor_from_p1700": p1700.get("EOM_bundle", {}),
        "reverse_identity_anchor_from_p1700": {
            "variational_identity_residuals": p1700.get("variational_identity_residuals", {}),
            "variational_identity_zero_flags": p1700.get("variational_identity_zero_flags", {}),
            "identity_status": p1700.get("identity_status", "UNKNOWN"),
        },
        "bidirectional_matrix": {
            "kernel_to_coefficients": "EXPORTED",
            "coefficients_to_full_lagrangian": "EXPORTED",
            "full_lagrangian_to_eom_bundle": "EXPORTED_REDUCED_PROXY",
            "eom_to_variational_origin": "LOCAL_IDENTITY_PASS_REDUCED_PROXY",
            "global_nonproxy_eom_to_lagrangian_inversion": "OPEN_THEOREM_REQUIRED",
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": [
            "nonproxy_covariant_metric_spinor_eom_export",
            "global_helmholtz_integrability_theorem",
            "counterterm_flow_renormalization_closure_theorem",
            "BRST_nilpotency_and_Cutkosky_unitarity_theorem",
            "background_family_independence_theorem",
            "strict_selector_source_or_symmetry_breaking_premise_for_QW2191",
        ],
        "next_honest_step": "Podnieść P1701 traceability matrix z reduced-proxy do pełnego nonproxy bundle i domknąć theorem-level QG closures bez false-pass.",
        "lay_summary": "Mamy już spójny ślad od kernela strict do pełnego lagranżianu i równań ruchu oraz lokalny test odwrotnej zgodności. Nadal brakuje globalnych dowodów nonproxy, dlatego końcowy status pozostaje otwarty.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
