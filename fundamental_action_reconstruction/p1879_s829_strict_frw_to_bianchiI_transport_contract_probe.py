#!/usr/bin/env python3
"""P1879 S829 strict FRW->Bianchi-I transport contract for ImM/defect system."""
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
    p1878 = load("p1878_s828_strict_one_loop_imaginary_amplitude_binding_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    lam, y, kappa2 = sp.symbols("lam y kappa2", positive=True, real=True)
    sigma2 = sp.symbols("sigma2", nonnegative=True, real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    lam_eff = lam * (1 + c1**2)
    y_eff = y * (1 + c0 / 2)

    # FRW tree amplitudes.
    Mss_frw = 6 * lam_eff
    Msf_frw = y_eff**2
    Mff_frw = y_eff**2 + kappa2

    # Minimal Bianchi-I transport deformation (anisotropy correction contract).
    a_ss, a_sf, a_ff = sp.symbols("a_ss a_sf a_ff", real=True)
    Mss_b1 = Mss_frw * (1 + a_ss * sigma2)
    Msf_b1 = Msf_frw * (1 + a_sf * sigma2)
    Mff_b1 = Mff_frw * (1 + a_ff * sigma2)

    transport_def_ss = sp.expand(Mss_b1 - Mss_frw)
    transport_def_sf = sp.expand(Msf_b1 - Msf_frw)
    transport_def_ff = sp.expand(Mff_b1 - Mff_frw)

    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1879",
        "stage_id": "S829",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {"p1878_present": "channel_defect_equations" in p1878},
        "strict_chain_step": "K_strict -> coefficients -> full-L/EOM channel amplitudes -> FRW->BianchiI transport contract",
        "effective_couplings": {"lam_eff": str(lam_eff), "y_eff": str(y_eff)},
        "frw_channel_amplitudes": {
            "Mss_frw": str(Mss_frw),
            "Msf_frw": str(Msf_frw),
            "Mff_frw": str(Mff_frw),
        },
        "bianchiI_transport_ansatz": {
            "Mss_b1": str(Mss_b1),
            "Msf_b1": str(Msf_b1),
            "Mff_b1": str(Mff_b1),
            "anisotropy_parameter": "sigma2",
            "note": "a_* must be fixed by explicit anisotropic one-loop integrals in same renormalization scheme.",
        },
        "transport_defect_equations": {
            "delta_ss": str(transport_def_ss),
            "delta_sf": str(transport_def_sf),
            "delta_ff": str(transport_def_ff),
            "frw_limit_condition": "sigma2=0 => all deltas vanish",
        },
        "qw2049_trace": {
            "Mss_frw_over_lam": str(sp.N((Mss_frw / lam).subs(qw2049), 12)),
            "Msf_frw_over_y2": str(sp.N((Msf_frw / y**2).subs(qw2049), 12)),
        },
        "strict_core_closure_missing_items": {
            "background_independence": "Need explicit Bianchi-I loop/discontinuity integrals fixing a_* and matching FRW limit.",
            "unitarity": "Need optical-defect closure on Bianchi-I and compatibility with FRW channel closure.",
            "renormalization": "Need same-scheme counterterm map that absorbs anisotropic divergences.",
            "selector_obstruction": "QW-2191 remains open.",
        },
        "false_pass_guard": "Transport ansatz is a contract, not a proven background-independence theorem.",
        "next_honest_step": "Compute anisotropic one-loop ImM/channel integrals to fix a_* and test Bianchi-I optical defects against FRW limit.",
        "lay_explanation": "To pierwszy techniczny krok do niezależności od tła: sprawdzamy, jak wyniki z kosmologii FRW mają przenosić się na lekko anizotropowe tło Bianchi-I.",
    }

    path = GEN / "p1879_s829_strict_frw_to_bianchiI_transport_contract_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
