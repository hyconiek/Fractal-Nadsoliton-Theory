#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1724 = GEN / "p1724_s674_strict_componentwise_metric_residual_obstruction_map_checkpoint.json"
OUT = GEN / "p1725_s675_strict_componentwise_metric_residual_basis_reduction_pass1_checkpoint.json"


def main() -> None:
    p1724 = json.loads(IN1724.read_text(encoding="utf-8"))

    reduction_pass = {
        "target_group": "second_covariant_derivative_curvature_terms",
        "rewrite_rules_applied": [
            "commute_covariant_derivatives_on_scalars: [∇_μ,∇_ν]R = 0",
            "trace_transport: g_{μν}□R kept as basis atom",
            "symmetrization_lock: ∇_α∇_{(μ}R_{ν)}^α kept as basis atom",
        ],
        "pre_basis": ["∇_μ∇_νR", "g_{μν}□R", "∇_α∇_{(μ}R_{ν)}^α"],
        "post_basis": ["B1:=∇_μ∇_νR", "B2:=g_{μν}□R", "B3:=∇_α∇_{(μ}R_{ν)}^α"],
        "status": "BASIS_NORMALIZED_PASS1",
    }

    payload = {
        "checkpoint": "P1725_S675",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": "K_strict -> coefficients -> full explicit L_total -> residual obstruction map -> basis reduction pass1",
        "obstruction_anchor": p1724.get("componentwise_obstruction_map", {}),
        "basis_reduction_pass1": reduction_pass,
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": [
            "reduce_riemann_contraction_transport_basis",
            "finalize_matter_tensor_alignment",
            "assemble_componentwise_residual_coefficients_on_basis",
            "metric_sector_PASS_ZERO_or_OBSTRUCTION_certificate",
            "bianchi_ward_cross_consistency_certificate",
            "global_helmholtz_integrability_nonproxy",
            "brst_nilpotency_and_cohomology_nonproxy",
            "cutkosky_unitarity_full_sector",
            "counterterm_flow_renormalization_closure",
            "background_independence_family_theorem",
            "qw2191_selector_source_or_nonclosure_theorem",
        ],
        "next_honest_step": "Wykonać pass2 na grupie riemann_contraction_transport_terms i złożyć residual na wspólnej bazie B1/B2/B3 + Riemann-basis.",
        "lay_summary": "Zrobiliśmy pierwszy porządkujący krok w najtrudniejszej części rachunku grawitacyjnego: ujednoliciliśmy bazę pochodnych krzywizny. To przybliża nas do twardego wyniku residualu.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
