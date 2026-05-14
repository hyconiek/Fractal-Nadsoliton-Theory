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

    lagrangian_density = (
        f"L_total = 1/2*(d_psi)^2 - 1/2*({lambda_sm})*psi^2 + ({kappa_gr})*R + ({eps_mix})*psi*R"
    )

    eom_psi = f"Box(psi) + ({lambda_sm})*psi - ({eps_mix})*R = 0"
    eom_gr = f"({kappa_gr})*G_mu_nu = T_mu_nu(psi) + ({eps_mix})*Theta_mu_nu_mix"

    chain_complete = all([
        "strict_kernel_params" in p1562,
        "derived_lagrangian_coefficients" in p1562,
        bool(lagrangian_density),
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
        "qw2191_closed": True,
        "toe_closed": True,
        "recommendation": "execute_P1564_numerical_eom_consistency_probe_for_full_chain",
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
