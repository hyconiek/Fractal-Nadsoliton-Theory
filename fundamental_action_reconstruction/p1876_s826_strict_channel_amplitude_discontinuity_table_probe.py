#!/usr/bin/env python3
"""P1876 S826 strict channel amplitude/discontinuity table and closure obligations probe."""
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
    p1875 = load("p1875_s825_strict_one_loop_counterterm_and_cutkosky_witness_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    lam, y, kappa2 = sp.symbols("lam y kappa2", positive=True, real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    lam_eff = sp.simplify(lam * (1 + c1**2))
    y_eff = sp.simplify(y * (1 + c0 / 2))

    # Strict-side tree proxies for selected channels.
    M_tree_ss_ss = 6 * lam_eff
    M_tree_sf_sf = y_eff**2
    M_tree_ff_ff = y_eff**2 + kappa2

    rho2_ss, rho2_sf, rho2_ff = sp.symbols("rho2_ss rho2_sf rho2_ff", nonnegative=True, real=True)
    ImM_ss, ImM_sf, ImM_ff = sp.symbols("ImM_ss ImM_sf ImM_ff", real=True)

    optical_defect_ss = ImM_ss - rho2_ss * M_tree_ss_ss**2
    optical_defect_sf = ImM_sf - rho2_sf * M_tree_sf_sf**2
    optical_defect_ff = ImM_ff - rho2_ff * M_tree_ff_ff**2

    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1876",
        "stage_id": "S826",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {"p1875_present": "one_loop_counterterm_dictionary_candidate" in p1875},
        "strict_chain_step": "K_strict -> effective couplings -> channel tree amplitudes -> optical discontinuity defects",
        "effective_couplings": {"lam_eff": str(lam_eff), "y_eff": str(y_eff)},
        "channel_amplitude_table": {
            "scalar_scalar_to_scalar_scalar": str(M_tree_ss_ss),
            "scalar_fermion_to_scalar_fermion": str(M_tree_sf_sf),
            "fermion_fermion_to_fermion_fermion": str(M_tree_ff_ff),
        },
        "cutkosky_discontinuity_table": {
            "defect_ss": str(optical_defect_ss),
            "defect_sf": str(optical_defect_sf),
            "defect_ff": str(optical_defect_ff),
            "strict_pass_condition": "defect_ss=defect_sf=defect_ff=0 after explicit one-loop discontinuity export",
        },
        "qw2049_trace": {
            "M_tree_ss_ss_qw2049_over_lam": str(sp.N((M_tree_ss_ss / lam).subs(qw2049), 12)),
            "M_tree_sf_sf_qw2049_over_y2": str(sp.N((M_tree_sf_sf / y**2).subs(qw2049), 12)),
            "M_tree_ff_ff_qw2049": str(sp.N(M_tree_ff_ff.subs({**qw2049, y: 1.0, kappa2: 1.0}), 12)),
        },
        "background_independence_obligation_matrix": {
            "frw_to_bianchiI_transport": "OPEN_MISSING_THEOREM",
            "frw_to_static_spherical_transport": "OPEN_MISSING_THEOREM",
            "atlas_cocycle_consistency": "OPEN_MISSING_WITNESS",
        },
        "strict_core_closure_missing_items": {
            "renormalization": "Need channel-matched loop integrals fixing finite parts and counterterm closure in same scheme.",
            "unitarity": "Need computed ImM_{ss,sf,ff} that set all defects to zero.",
            "background_independence": "Need transport theorem proving witness persistence across atlas charts.",
            "selector_obstruction": "QW-2191 remains open.",
        },
        "false_pass_guard": "Channel table and defect equations are obligations, not final closure witnesses.",
        "next_honest_step": "Export explicit one-loop ImM tables per channel and evaluate defect_ss/sf/ff with the same renormalization scheme.",
        "lay_explanation": "To krok od ogólnej formuły do konkretnych kanałów rozpraszania: zapisujemy, co dokładnie trzeba policzyć, żeby potwierdzić unitarność i spójność kwantową teorii.",
    }

    path = GEN / "p1876_s826_strict_channel_amplitude_discontinuity_table_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
