#!/usr/bin/env python3
"""P1887 S837 strict nontrivial transport branch probe with shared loop data."""
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
    p1886 = load("p1886_s836_strict_joint_equation_symbolic_solution_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    lam, y, kappa2 = sp.symbols("lam y kappa2", positive=True, real=True)
    eps, mu2, s = sp.symbols("eps mu2 s", positive=True, real=True)
    sigma2 = sp.symbols("sigma2", nonnegative=True, real=True)
    nu = sp.symbols("nu", real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    lam_eff = lam * (1 + c1**2)
    y_eff = y * (1 + c0 / 2)

    I = 1 / (16 * sp.pi**2) * (1 / eps + sp.log(mu2 / s))
    Mss = 6 * lam_eff * (1 + lam_eff * I)
    Msf = y_eff**2 * (1 + y_eff**2 * I)
    Mff = (y_eff**2 + kappa2) * (1 + y_eff**2 * I)

    # Nontrivial branch: shared anisotropy coefficient nu (instead of a*=0)
    Mss_b1 = Mss * (1 + nu * sigma2)
    Msf_b1 = Msf * (1 + nu * sigma2)
    Mff_b1 = Mff * (1 + nu * sigma2)

    delta_ss = sp.expand(Mss_b1 - Mss)
    delta_sf = sp.expand(Msf_b1 - Msf)
    delta_ff = sp.expand(Mff_b1 - Mff)

    consistency_ratio = {
        "delta_ss_over_Mss": str(sp.simplify(delta_ss / Mss)),
        "delta_sf_over_Msf": str(sp.simplify(delta_sf / Msf)),
        "delta_ff_over_Mff": str(sp.simplify(delta_ff / Mff)),
    }

    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1887",
        "stage_id": "S837",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {"p1886_present": "symbolic_solution_candidate" in p1886},
        "strict_chain_step": "K_strict -> coefficients -> renormalized amplitudes -> nontrivial FRW/BianchiI transport branch",
        "nontrivial_transport_branch": {
            "Mss_b1": str(Mss_b1),
            "Msf_b1": str(Msf_b1),
            "Mff_b1": str(Mff_b1),
            "branch_parameter": "nu",
            "branch_note": "Shared nonzero nu explores nontrivial transport beyond minimal a*=0 branch.",
        },
        "transport_deltas": {
            "delta_ss": str(delta_ss),
            "delta_sf": str(delta_sf),
            "delta_ff": str(delta_ff),
            "relative_form": consistency_ratio,
        },
        "qw2049_trace": {
            "lam_eff_over_lam": str(sp.N((lam_eff / lam).subs(qw2049), 12)),
            "y_eff_over_y": str(sp.N((y_eff / y).subs(qw2049), 12)),
        },
        "strict_core_closure_missing_items": {
            "background_independence": "Need theorem/witness that nonzero nu branch preserves joint closure conditions.",
            "unitarity": "Need Cutkosky defects recomputed on nontrivial transport branch.",
            "renormalization": "Need same-scheme finite-part consistency on transported amplitudes.",
            "selector_obstruction": "QW-2191 remains open.",
        },
        "false_pass_guard": "Nontrivial branch construction is exploratory and not a closure theorem.",
        "next_honest_step": "Solve joint Cutkosky+transport equations on nu!=0 branch with explicit loop/discontinuity inputs.",
        "lay_explanation": "To krok dalej niż trywialne przeniesienie: sprawdzamy, czy teoria może pozostać spójna także przy niezerowej deformacji tła.",
    }

    path = GEN / "p1887_s837_strict_nontrivial_transport_branch_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
