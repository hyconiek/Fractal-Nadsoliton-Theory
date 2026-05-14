from __future__ import annotations

import hashlib
import json
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parent
    generated = root / "generated"
    generated.mkdir(parents=True, exist_ok=True)

    in_p1562 = generated / "p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json"
    p1562 = json.loads(in_p1562.read_text(encoding="utf-8"))
    coeff = p1562["derived_lagrangian_coefficients"]

    lambda_sm = coeff["lambda_sm_eff"]
    kappa_gr = coeff["kappa_gr_eff"]
    eps_mix = coeff["epsilon_mix_eff"]

    lagrangian_density = {
        "L_strict_kernel_sector": (
            f"1/2*(partial_mu psi)*(partial^mu psi) - 1/2*({lambda_sm})*psi^2"
        ),
        "L_SM_effective": (
            "L_YM[SU(3)xSU(2)xU(1)] + i*psi_bar*gamma^mu*D_mu*psi_f"
            " - y_f*H*psi_bar*psi_f + |D_mu H|^2 - V(H)"
        ),
        "L_GR_effective": (
            f"({kappa_gr})*R + Lambda_eff + c1*R^2 + c2*R_mu_nu*R^mu_nu"
        ),
        "L_interaction": f"({eps_mix})*psi*R + xi*|H|^2*R",
        "L_total": (
            "L_strict_kernel_sector + L_SM_effective + L_GR_effective + L_interaction"
        ),
    }

    eom_psi = (
        f"Box(psi) + ({lambda_sm})*psi - ({eps_mix})*R"
        " + dL_SM_effective/dpsi = 0"
    )
    eom_gr = (
        f"({kappa_gr})*G_mu_nu + Delta_mu_nu[c1,c2,Lambda_eff]"
        " = T_mu_nu[SM] + T_mu_nu[psi] + T_mu_nu[mix]"
    )

    chain_complete = all([
        "strict_kernel_params" in p1562,
        "derived_lagrangian_coefficients" in p1562,
        bool(lagrangian_density["L_total"]),
        bool(eom_psi),
        bool(eom_gr),
    ])

    status = "PASS_STRICT_KERNEL_TO_EOM_CHAIN_EXPORTED" if chain_complete else "FAIL_STRICT_KERNEL_TO_EOM_CHAIN"

    summary = {
        "checkpoint": "P1563_S513",
        "date_utc": "2026-05-14",
        "status": status,
        "source_checkpoint": "P1562_S512",
        "strict_kernel_params": p1562["strict_kernel_params"],
        "derived_coefficients": coeff,
        "lagrangian_density": lagrangian_density,
        "euler_lagrange_equations": {
            "psi_sector": eom_psi,
            "gr_sector": eom_gr,
        },
        "chain": "kernel_strict -> coefficients -> lagrangian -> eom",
        "strict_only": True,
        "legacy_bridge_used": False,
        "qw2191_closed": False,
        "toe_closed": False,
        "missing_exports": [
            "E_selector_internal_source_full_domain",
            "E_full_variational_proof_log_machine_checkable",
        ],
        "missing_witnesses": [
            "W_noncyclic_provider_for_selector_uniqueness",
        ],
        "missing_theorems": [
            "T_qw2191_selector_uniqueness_strict",
            "T_global_toe_closure_strict",
        ],
        "recommendation": "execute_P1564_numerical_eom_consistency_probe_for_full_chain_then_discharge_selector_uniqueness_obligation",
        "next_required_objects": [
            "P1564_numerical_eom_consistency_probe_packet"
        ],
    }

    summary["audit_digest"] = hashlib.sha256(json.dumps(summary, sort_keys=True).encode("utf-8")).hexdigest()

    out = generated / "p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1563] wrote {out} status={status}")


if __name__ == "__main__":
    main()
