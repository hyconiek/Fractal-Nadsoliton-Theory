#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1694 = GEN / "p1694_s644_strict_kernel_to_full_lagrangian_bidirectional_map_witness.json"
IN1699 = GEN / "p1699_s649_strict_full_lagrangian_sector_eom_bundle_spinor_variational_fix_checkpoint.json"
IN1727 = GEN / "p1727_s677_strict_componentwise_metric_residual_coefficient_vector_stub_checkpoint.json"
OUT = GEN / "p1728_s678_strict_full_lagrangian_non_skeleton_and_bidirectional_closure_gap_checkpoint.json"


def _inject(template: str, repl: dict[str, str]) -> str:
    out = template
    for k, v in repl.items():
        out = out.replace(k, v)
    return out


def main() -> None:
    p1694 = json.loads(IN1694.read_text(encoding="utf-8"))
    p1699 = json.loads(IN1699.read_text(encoding="utf-8"))
    p1727 = json.loads(IN1727.read_text(encoding="utf-8"))

    coeff = p1694.get("forward_coefficient_map_symbolic", {})

    repl = {
        "M_Pl^2": coeff.get("Mpl2", "M_Pl^2"),
        "c_R2": coeff.get("cR2", "c_R2"),
        "c_Ric2": coeff.get("cRic2", "c_Ric2"),
        "c_Riem2": coeff.get("cRiem2", "c_Riem2"),
        "μ_H^2": coeff.get("muH2", "μ_H^2"),
        "λ_H": coeff.get("lambdaH", "λ_H"),
        "ξ_HR": coeff.get("xiHR", "ξ_HR"),
    }

    anchor = p1694.get("full_lagrangian_density_anchor", {})
    full_lagrangian_instantiated = {k: _inject(v, repl) for k, v in anchor.items()}

    payload = {
        "checkpoint": "P1728_S678",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": "K_strict -> coefficients -> full instantiated nonskeleton L_total -> EOM bundle -> bidirectional closure gap map",
        "kernel": p1694.get("kernel", "K_strict(d)=cos(omega*d+phi)/(1+beta*d^eta)"),
        "kernel_input": p1694.get("kernel_input", {}),
        "forward_coefficient_map_symbolic": coeff,
        "full_lagrangian_density_nonskeleton_instantiated": full_lagrangian_instantiated,
        "eom_bundle_anchor": p1699.get("reduced_bundle_with_spinor_fix", {}),
        "metric_residual_anchor": p1727.get("residual_coefficient_vector_stub", {}),
        "bidirectional_witness_map": {
            "forward_kernel_to_coefficients": "EXPORTED",
            "forward_coefficients_to_full_lagrangian": "EXPORTED_NONSKELETON_INSTANTIATED",
            "forward_lagrangian_to_eom": "PARTIAL_EXPORTED_NONPROXY_PENDING_FULL_METRIC_CERTIFICATE",
            "reverse_eom_to_variational_origin": "OPEN_THEOREM_REQUIRED",
            "reverse_variational_origin_to_kernel_identifiability": "OPEN_THEOREM_REQUIRED",
        },
        "strict_core_closure_gate": {
            "renormalization_counterterm_flow": "OPEN_WITNESS_REQUIRED",
            "unitarity_cutkosky_full_sector": "OPEN_WITNESS_REQUIRED",
            "background_independence_family": "OPEN_WITNESS_REQUIRED",
            "brst_nilpotency_cohomology": "OPEN_WITNESS_REQUIRED",
            "qw2191_selector_source": "OPEN_SELECTOR_PREMISE_OR_THEOREM_REQUIRED",
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Wyeksportować pierwszy niezerowy/zerowy jawny EL_g-E_{μν} residual vector (k_B1..k_C2) dla konkretnej rodziny metryk oraz dołączyć theorem-level witness map dla renormalizacji/unitarności/background-independence.",
        "lay_summary": "Mamy już pełny strict tor do jawnego lagranżianu (bez szkieletu) i do pakietu równań ruchu. Brakuje jeszcze dowodów domknięcia kwantowej grawitacji i pełnego kroku wstecznego z równań do źródła wariacyjnego.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
