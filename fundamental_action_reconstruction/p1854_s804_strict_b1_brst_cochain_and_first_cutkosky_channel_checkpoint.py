#!/usr/bin/env python3
"""P1854 S804 strict B1 BRST cochain and first Cutkosky channel checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    p1853 = load("p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json")
    p1852 = load("p1852_s802_strict_b1_brst_anomaly_and_cutkosky_seed_witness_checkpoint.json")

    eval_coeffs = ((p1853.get("b1_symbolic_evaluation") or {}).get("evaluated_coefficients") or {})

    # First explicit BRST anomaly-cochain assembly on B1 (still pre-theorem).
    brst_cochain_b1 = {
        "anomaly_polynomial": "A_B1 = k1*c*R^2 + k2*c*Ricci^2 + k3*c*Riemann^2 + k4*c*G_GB + k5*c*Tr(F^2)",
        "cochain_coeff_mapping": {
            "k1": (eval_coeffs.get("a_R2") or {}).get("symbolic", "A_B1_k1_from_a_R2"),
            "k2": (eval_coeffs.get("a_Ric2") or {}).get("symbolic", "A_B1_k2_from_a_Ric2"),
            "k3": (eval_coeffs.get("a_Riem2") or {}).get("symbolic", "A_B1_k3_from_a_Riem2"),
            "k4": (eval_coeffs.get("a_GB") or {}).get("symbolic", "A_B1_k4_from_a_GB"),
            "k5": "A_B1_k5_from_gauge_fermion_triangle_sector",
        },
        "target": "A_B1 == 0",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    cutkosky_first_channel = {
        "channel": "graviton -> gauge_gauge",
        "optical_identity_form": "Disc M(grav->grav)|_gauge_loop = Integral dPi_gg * M(grav->gg) * M*(grav->gg)",
        "normalization": "MSbar_B1_seed",
        "residue_positivity_requirement": "rho_gg(s) >= 0 on physical cut",
        "first_quantitative_proxy": {
            "input_from_b1_coefficients": {
                "a_R2_num": (eval_coeffs.get("a_R2") or {}).get("numeric_20d", "OPEN"),
                "a_Ric2_num": (eval_coeffs.get("a_Ric2") or {}).get("numeric_20d", "OPEN"),
            },
            "proxy_statement": "If projected gauge-loop discontinuity coefficient has positive sign under declared conventions, channel remains unitarity-compatible at seed level.",
            "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        },
    }

    out = {
        "packet_id": "P1854",
        "stage_id": "S804",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1853_present": "b1_symbolic_evaluation" in p1853,
            "p1852_present": "brst_anomaly_seed_contract" in p1852,
        },
        "strict_chain_extension": "K_strict -> B1 evaluated coefficients -> BRST cochain coefficients -> first Cutkosky channel proxy",
        "brst_cochain_b1": brst_cochain_b1,
        "cutkosky_first_channel": cutkosky_first_channel,
        "proven": "First explicit B1 cochain-coefficient mapping and first unitarity channel proxy are exported on strict-only lane.",
        "open": "A_B1 exact cancellation and full discontinuity integral evaluation remain open theorem tasks.",
        "false_pass_risk": "Cochain mapping and channel proxy are not equivalent to full BRST nilpotency or Cutkosky theorem discharge.",
        "next_honest_step": "Compute k5 from explicit gauge-fermion anomaly sector and evaluate full graviton->gauge_gauge discontinuity integral with positivity certificate.",
        "lay_explanation": "Mamy już pierwszy konkretny szkic liczenia anomalii i unitarności na jednym kanale, ale nadal trzeba policzyć pełne całki i pokazać ścisły dowód.",
    }

    path = GEN / "p1854_s804_strict_b1_brst_cochain_and_first_cutkosky_channel_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
