#!/usr/bin/env python3
"""P1884 S834 strict explicit loop-integral stub evaluation on shared QG checks."""
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
    p1883 = load("p1883_s833_strict_one_loop_finitepart_joint_consistency_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    lam, y, kappa2 = sp.symbols("lam y kappa2", positive=True, real=True)
    eps, mu2, s = sp.symbols("eps mu2 s", positive=True, real=True)
    sigma2 = sp.symbols("sigma2", nonnegative=True, real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    lam_eff = lam * (1 + c1**2)
    y_eff = y * (1 + c0 / 2)

    # Explicit 1-loop integral stubs (MS-like pole + finite log part)
    I_lam = 1 / (16 * sp.pi**2) * (1 / eps + sp.log(mu2 / s))
    I_y = 1 / (16 * sp.pi**2) * (1 / eps + sp.log(mu2 / s) + sp.Rational(1, 2))

    Ct_lam = lam_eff * I_lam
    Ct_y = y_eff**2 * I_y

    Mss_tree = 6 * lam_eff
    Msf_tree = y_eff**2
    Mff_tree = y_eff**2 + kappa2

    Mss_ren = sp.expand(Mss_tree * (1 + Ct_lam))
    Msf_ren = sp.expand(Msf_tree * (1 + Ct_y))
    Mff_ren = sp.expand(Mff_tree * (1 + Ct_y))

    rho2_ss, rho2_sf, rho2_ff = sp.symbols("rho2_ss rho2_sf rho2_ff", nonnegative=True, real=True)
    ImM_ss, ImM_sf, ImM_ff = sp.symbols("ImM_ss ImM_sf ImM_ff", real=True)

    defect_ss = sp.expand(ImM_ss - rho2_ss * Mss_ren**2)
    defect_sf = sp.expand(ImM_sf - rho2_sf * Msf_ren**2)
    defect_ff = sp.expand(ImM_ff - rho2_ff * Mff_ren**2)

    a_ss, a_sf, a_ff = sp.symbols("a_ss a_sf a_ff", real=True)
    transport_ss = sp.expand(Mss_ren * (1 + a_ss * sigma2) - Mss_ren)
    transport_sf = sp.expand(Msf_ren * (1 + a_sf * sigma2) - Msf_ren)
    transport_ff = sp.expand(Mff_ren * (1 + a_ff * sigma2) - Mff_ren)

    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1884",
        "stage_id": "S834",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {"p1883_present": "renormalized_channel_amplitudes" in p1883},
        "strict_chain_step": "K_strict -> coefficients -> explicit loop-integral stubs -> renormalized amplitudes -> Cutkosky/transport joint checks",
        "explicit_loop_integral_stubs": {
            "I_lam": str(I_lam),
            "I_y": str(I_y),
            "Ct_lam": str(Ct_lam),
            "Ct_y": str(Ct_y),
            "note": "Stubs include explicit finite log parts; full diagrammatic coefficients still pending.",
        },
        "renormalized_amplitudes": {
            "Mss_ren": str(Mss_ren),
            "Msf_ren": str(Msf_ren),
            "Mff_ren": str(Mff_ren),
        },
        "joint_checks": {
            "cutkosky_defects": {"ss": str(defect_ss), "sf": str(defect_sf), "ff": str(defect_ff)},
            "transport_deltas": {"ss": str(transport_ss), "sf": str(transport_sf), "ff": str(transport_ff)},
            "frw_limit": "sigma2=0",
        },
        "qw2049_trace": {
            "lam_eff_over_lam": str(sp.N((lam_eff / lam).subs(qw2049), 12)),
            "y_eff_over_y": str(sp.N((y_eff / y).subs(qw2049), 12)),
        },
        "strict_core_closure_missing_items": {
            "renormalization": "Need full diagrammatic combinatorics and finite parts beyond stubs.",
            "unitarity": "Need explicit ImM from the same loop data.",
            "background_independence": "Need Bianchi-I loop discontinuities consistent with FRW stub limit.",
            "selector_obstruction": "QW-2191 remains open.",
        },
        "false_pass_guard": "Explicit integral stubs are not full computed witness set.",
        "next_honest_step": "Replace stub integrals with full diagram-resolved values and solve defects/deltas jointly.",
        "lay_explanation": "Dodaliśmy już jawny kształt całek z częścią skończoną (logarytmy), ale to nadal wersja robocza przed pełnym rachunkiem diagramów.",
    }

    path = GEN / "p1884_s834_strict_explicit_loop_integral_stub_evaluation_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
