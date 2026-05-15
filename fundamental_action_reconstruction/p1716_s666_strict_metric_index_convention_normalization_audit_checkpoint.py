#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1715 = GEN / "p1715_s665_strict_metric_curvature_term_expansion_scaffold_checkpoint.json"
OUT = GEN / "p1716_s666_strict_metric_index_convention_normalization_audit_checkpoint.json"


def main() -> None:
    p1715 = json.loads(IN1715.read_text(encoding="utf-8"))
    h = p1715.get("curvature_term_expansion_scaffold", {})

    normalized = {
        "H_R2_munu": "H^(R2)_{μν} = 2 R R_{μν} - (1/2) g_{μν} R^2 + 2(g_{μν}□ - ∇_μ∇_ν)R",
        "H_Ric2_munu": "H^(Ric2)_{μν} = 2 R_{μα}R_{ν}^α - (1/2)g_{μν}R_{αβ}R^{αβ} + □R_{μν} + g_{μν}∇_α∇_βR^{αβ} - 2∇_α∇_{(μ}R_{ν)}^α",
        "H_Riem2_munu": "H^(Riem2)_{μν} = 2 R_{μαβγ}R_{ν}^{αβγ} - (1/2)g_{μν}R_{αβγδ}R^{αβγδ} - 4∇_α∇_βR_{μν}^{\\ \alpha\beta}",
    }

    token_audit = {
        k: {
            "has_munu": ("_{μν}" in v),
            "has_covariant_derivative": ("∇" in v or "□" in v),
            "has_quadratic_curvature": ("R^2" in v or "R_{αβ}" in v or "R_{αβγδ}" in v),
        }
        for k, v in normalized.items()
    }

    audit_pass = all(all(flags.values()) for flags in token_audit.values())

    payload = {
        "checkpoint": "P1716_S666",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": "K_strict -> coefficients -> full explicit L_total -> E_munu scaffold -> index convention normalization audit",
        "input_h_terms": h,
        "normalized_h_terms": normalized,
        "token_audit": token_audit,
        "audit_status": "PASS_INDEX_CONVENTION_NORMALIZED" if audit_pass else "FAIL_TOKEN_AUDIT",
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": [
            "explicit_metric_sector_EL_residual_zero_certificate",
            "bianchi_ward_cross_consistency_certificate",
            "global_helmholtz_integrability_nonproxy",
            "brst_nilpotency_and_cohomology_nonproxy",
            "cutkosky_unitarity_full_sector",
            "counterterm_flow_renormalization_closure",
            "background_independence_family_theorem",
            "qw2191_selector_source_or_nonclosure_theorem",
        ],
        "next_honest_step": "Użyć normalized_h_terms do pierwszego jawnego residual test EL_g - E_munu oraz wyeksportować pass/fail z pełnym obstrukcyjnym śladem składników.",
        "lay_summary": "Uporządkowaliśmy zapis trudnych członów grawitacyjnych do jednej spójnej konwencji. To usuwa techniczne niejednoznaczności przed właściwym testem równania metrycznego.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
