#!/usr/bin/env python3
"""P1883 S833 strict one-loop finite-part joint consistency probe."""
from __future__ import annotations
import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    p1881 = load("p1881_s831_strict_termwise_one_loop_binding_matrix_probe.json")
    p1879 = load("p1879_s829_strict_frw_to_bianchiI_transport_contract_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    m2, lam, y, kappa2 = sp.symbols("m2 lam y kappa2", positive=True, real=True)
    eps = sp.symbols("eps", positive=True, real=True)
    sigma2 = sp.symbols("sigma2", nonnegative=True, real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    lam_eff = lam * (1 + c1**2)
    y_eff = y * (1 + c0 / 2)

    # One-loop coefficients split into pole + finite parts.
    z_lam, z_y = sp.symbols("z_lam z_y", real=True)
    f_lam, f_y, f_ff = sp.symbols("f_lam f_y f_ff", real=True)

    Ct_lam = z_lam / eps + f_lam
    Ct_y = z_y / eps + f_y
    Ct_ff = z_y / eps + f_ff

    Mss_tree = 6 * lam_eff
    Msf_tree = y_eff**2
    Mff_tree = y_eff**2 + kappa2

    # Renormalized amplitudes (minimal joint ansatz)
    Mss_ren = Mss_tree * (1 + Ct_lam)
    Msf_ren = Msf_tree * (1 + Ct_y)
    Mff_ren = Mff_tree * (1 + Ct_ff)

    rho2_ss, rho2_sf, rho2_ff = sp.symbols("rho2_ss rho2_sf rho2_ff", nonnegative=True, real=True)
    ImM_ss, ImM_sf, ImM_ff = sp.symbols("ImM_ss ImM_sf ImM_ff", real=True)

    defect_ss = sp.expand(ImM_ss - rho2_ss * Mss_ren**2)
    defect_sf = sp.expand(ImM_sf - rho2_sf * Msf_ren**2)
    defect_ff = sp.expand(ImM_ff - rho2_ff * Mff_ren**2)

    # FRW->Bianchi-I transport on same renormalized amplitudes.
    a_ss, a_sf, a_ff = sp.symbols("a_ss a_sf a_ff", real=True)
    Mss_b1 = Mss_ren * (1 + a_ss * sigma2)
    Msf_b1 = Msf_ren * (1 + a_sf * sigma2)
    Mff_b1 = Mff_ren * (1 + a_ff * sigma2)

    transport_ss = sp.expand(Mss_b1 - Mss_ren)
    transport_sf = sp.expand(Msf_b1 - Msf_ren)
    transport_ff = sp.expand(Mff_b1 - Mff_ren)

    # EOM stability proxy on same data (delta couplings must stay finite/small-controlled).
    delta_lam = sp.expand(lam_eff * Ct_lam)
    delta_y = sp.expand(y_eff * Ct_y)

    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1883",
        "stage_id": "S833",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1881_present": "termwise_one_loop_binding_matrix" in p1881,
            "p1879_present": "bianchiI_transport_ansatz" in p1879,
        },
        "strict_chain_step": "K_strict -> coefficients -> full L_total one-loop (pole+finite) -> renormalized amplitudes -> Cutkosky defects + FRW/BianchiI transport + EOM stability",
        "renormalization_data": {
            "Ct_lam": str(Ct_lam),
            "Ct_y": str(Ct_y),
            "Ct_ff": str(Ct_ff),
            "finite_part_note": "f_* are explicit finite parts required beyond 1/eps poles.",
        },
        "renormalized_channel_amplitudes": {
            "Mss_ren": str(Mss_ren),
            "Msf_ren": str(Msf_ren),
            "Mff_ren": str(Mff_ren),
        },
        "cutkosky_defect_checks": {
            "defect_ss": str(defect_ss),
            "defect_sf": str(defect_sf),
            "defect_ff": str(defect_ff),
            "pass_condition": "all defects == 0 with computed ImM and fixed (z_*, f_*)",
        },
        "frw_bianchiI_transport_checks": {
            "transport_ss": str(transport_ss),
            "transport_sf": str(transport_sf),
            "transport_ff": str(transport_ff),
            "frw_limit": "sigma2=0",
        },
        "eom_stability_proxy": {
            "delta_lam": str(delta_lam),
            "delta_y": str(delta_y),
            "criterion": "finite-part-controlled shifts must preserve same EOM residual class",
        },
        "qw2049_trace": {
            "lam_eff_over_lam": str(sp.N((lam_eff / lam).subs(qw2049), 12)),
            "y_eff_over_y": str(sp.N((y_eff / y).subs(qw2049), 12)),
        },
        "strict_core_closure_missing_items": {
            "renormalization": "Need explicit loop integrals fixing z_* and finite parts f_*.",
            "unitarity": "Need explicit ImM channel computations to satisfy defect equations.",
            "background_independence": "Need computed Bianchi-I discontinuities compatible with FRW limit.",
            "selector_obstruction": "QW-2191 remains open.",
        },
        "false_pass_guard": "Joint consistency frame is not closure proof until z_*, f_*, ImM are computed.",
        "next_honest_step": "Compute explicit loop integrals for (z_*,f_*) and evaluate all three blocks (Cutkosky/transport/EOM) on one shared dataset.",
        "lay_explanation": "To pierwszy wspólny test: te same poprawki kwantowe mają jednocześnie zamknąć unitarność, transport między tłami i stabilność równań ruchu.",
    }

    path = GEN / "p1883_s833_strict_one_loop_finitepart_joint_consistency_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
