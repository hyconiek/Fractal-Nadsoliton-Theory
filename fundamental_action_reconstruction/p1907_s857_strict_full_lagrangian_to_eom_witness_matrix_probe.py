#!/usr/bin/env python3
"""P1907 S857 strict full-Lagrangian to EOM witness matrix probe."""
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
    p1906 = load("p1906_s856_strict_c1_gr_diagram_inventory_stub_probe.json")

    out = {
        "packet_id": "P1907",
        "stage_id": "S857",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1906_present": "diagram_inventory_stub" in p1906,
            "p1906_next_missing_section": p1906.get("strict_core_closure_missing_items", {}).get("next_missing_section"),
        },
        "strict_chain_forward": [
            "K_strict(omega,phi,beta,eta)",
            "coefficient_exports: {lambda_i, y_f, xi_H, kappa, gauge couplings}",
            "full L_SM+L_GR (non-skeleton, termwise)",
            "EOM set: {matter, gauge, Higgs, metric}",
            "closure witnesses: {renormalization, unitarity, background_independence}",
        ],
        "strict_chain_reverse": [
            "target observables + Ward/Slavnov-Taylor contracts",
            "EOM consistency + counterterm closure",
            "termwise Lagrangian coefficient constraints",
            "strict coefficient table constraints",
            "kernel tuple admissibility region",
        ],
        "full_lagrangian_term_registry_non_skeleton": {
            "sm_gauge_sector": ["-1/4 F^2_SU3", "-1/4 F^2_SU2", "-1/4 F^2_U1"],
            "sm_fermion_kinetic_sector": ["i psi_bar gamma^mu D_mu psi (all generations)"],
            "higgs_sector": ["|D_mu H|^2", "-mu_H^2 H^dagger H", "-lambda_H (H^dagger H)^2"],
            "yukawa_sector": ["-y_u qbar H_tilde u", "-y_d qbar H d", "-y_e lbar H e + h.c."],
            "gravity_sector": ["(M_Pl^2/2) R", "Lambda", "R^2", "R_{mu nu}R^{mu nu}", "R_{mu nu rho sigma}R^{mu nu rho sigma}"],
            "nonminimal_mixed_sector": ["-xi_H H^dagger H R", "dimension-6 EFT placeholders"],
            "note": "Registry is explicit and non-skeleton at sector level; coefficient tables and loop-normalization witnesses still missing.",
        },
        "eom_export_matrix": [
            {"field": "H", "source_terms": ["D^2 H", "mu_H^2 H", "lambda_H(H^dagger H)H", "xi_H R H", "Yukawa backreaction"], "status": "OPEN_SYMBOLIC"},
            {"field": "A_mu^a", "source_terms": ["D_nu F^{nu mu}", "matter currents", "ghost/BRST constraints"], "status": "OPEN_SYMBOLIC"},
            {"field": "psi_f", "source_terms": ["i gamma^mu D_mu psi_f", "Yukawa mass mixing"], "status": "OPEN_SYMBOLIC"},
            {"field": "g_{mu nu}", "source_terms": ["G_{mu nu}", "Lambda g_{mu nu}", "higher-curvature tensors", "T_{mu nu}^{SM}", "xi_H improvement"], "status": "OPEN_SYMBOLIC"},
        ],
        "missing_witness_exports": {
            "renormalization": "MISSING: explicit divergent coefficient table per diagram id + counterterm cancellation proof",
            "unitarity": "MISSING: Cutkosky/ImM channel-by-channel equality with common scheme",
            "background_independence": "MISSING: FRW/Bianchi-I transport witness beyond ansatz to covariant closure",
            "selector_qw2191": "MISSING: strict-core selector source or explicit symmetry-breaking premise export",
        },
        "false_pass_guard": "Full-sector registry and symbolic EOM scaffolding are not closure; numeric/symbolic witness exports remain mandatory.",
        "next_honest_step": "Export P1908 divergent_coefficients_table_v1 with explicit loop integral basis, renormalization scheme tags, and first PASS/FAIL counterterm rows tied to P1906 diagram ids.",
        "lay_explanation": "Mamy już pełną mapę składników teorii (jądro -> współczynniki -> pełny lagranżian -> równania ruchu), ale nadal brakuje twardych obliczeń pętlowych, które rozstrzygają czy teoria naprawdę zamyka się kwantowo i grawitacyjnie.",
    }

    path = GEN / "p1907_s857_strict_full_lagrangian_to_eom_witness_matrix_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
