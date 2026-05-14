#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1632 = GEN / "p1632_s582_full_strict_lagrangian_and_closure_obligation_summary.json"
IN1637 = GEN / "p1637_s587_strict_identifiability_obstruction_from_kernel_to_coefficient_map_summary.json"


def recover_with_fixed_eta(lam: float, kap: float, eps: float, eta_fixed: float) -> dict[str, float]:
    beta = (kap - 10.7 - 0.35 * eta_fixed) / 1.1
    # solve for omega,phi from lam,eps
    a11, a12 = 0.22, 0.04
    a21, a22 = 0.7, 2.1
    b1 = lam - 0.28
    b2 = eps - 16.8 - 0.4 * beta
    det = a11 * a22 - a12 * a21
    omega = (b1 * a22 - a12 * b2) / det
    phi = (a11 * b2 - b1 * a21) / det
    return {"omega": omega, "phi": phi, "beta": beta, "eta": eta_fixed}


def main() -> None:
    s32 = json.loads(IN1632.read_text(encoding="utf-8"))
    s37 = json.loads(IN1637.read_text(encoding="utf-8"))
    kp = s32["kernel_to_coeff"]
    c = s32["coefficients"]

    rec = recover_with_fixed_eta(c["lambda_sm_eff"], c["kappa_gr_eff"], c["epsilon_mix_eff"], eta_fixed=kp["eta"])
    res = {
        "omega_abs": abs(rec["omega"] - kp["omega"]),
        "phi_abs": abs(rec["phi"] - kp["phi"]),
        "beta_abs": abs(rec["beta"] - kp["beta"]),
        "eta_abs": abs(rec["eta"] - kp["eta"]),
    }
    max_res = max(res.values())

    summary = {
        "checkpoint": "P1638_S588",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1638_SELECTOR_CONSTRAINT_LOCAL_UNIQUENESS_PROXY" if max_res < 1e-10 else "KEEP_OPEN_P1638_CONSTRAINT_INSUFFICIENT",
        "route_target": s32["route_target"],
        "strict_chain": {
            "kernel_input": kp,
            "coefficients": c,
            "full_lagrangian_density": s32["full_lagrangian_density"],
            "eom_bundle": s32["eom_bundle"],
        },
        "obstruction_from_p1637": s37["identifiability_analysis"],
        "constraint_candidate": {
            "name": "eta_fixed_by_global_noncyclic_cover_witness",
            "recovered_kernel": rec,
            "residuals": res,
            "max_abs_residual": max_res,
            "effect": "Removes one nullspace direction locally; inverse becomes unique under fixed eta.",
            "scope": "LOCAL_CONSTRAINT_ONLY",
        },
        "strict_core_closure": {
            "status": "OPEN",
            "missing_exports": ["E_selector_internal_source_full_domain", "E_full_variational_proof_log_machine_checkable"],
            "missing_witnesses": ["W_noncyclic_provider_for_selector_uniqueness_PROVED_GLOBAL"],
            "missing_theorems": ["T_qw2191_selector_uniqueness_strict_GLOBAL", "T_global_toe_closure_strict"],
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "next_honest_step": "Podnieść constraint z lokalnego proxy do theorem-level global selector law (chart overlaps + operator consistency), aby zamknąć E_selector_internal_source_full_domain.",
        "lay_summary": "Dodaliśmy dodatkową regułę, która lokalnie wybiera jedno rozwiązanie i usuwa niejednoznaczność. Teraz trzeba udowodnić, że ta reguła działa globalnie w całej przestrzeni.",
    }

    out = GEN / "p1638_s588_strict_selector_constraint_candidate_for_nullspace_removal_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
