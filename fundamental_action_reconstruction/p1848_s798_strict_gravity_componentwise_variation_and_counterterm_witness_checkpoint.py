#!/usr/bin/env python3
"""P1848 S798 strict gravity componentwise variation and counterterm witness checkpoint."""

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
    p1847 = load("p1847_s797_strict_gravity_sector_nonproxy_density_and_covariant_eom_checkpoint.json")
    p1716 = load("p1716_s666_strict_metric_index_convention_normalization_audit_checkpoint.json")
    p1663 = load("p1663_s613_strict_qg_obligation_witness.json")

    h_r2_default = "H^(R2)_{mu nu} = 2 R R_{mu nu} - (1/2) g_{mu nu} R^2 + 2(g_{mu nu} Box - nabla_mu nabla_nu)R"
    h_r2 = ((p1716.get("normalized_h_terms") or {}).get("H_R2_munu")) or ((p1716.get("input_h_terms") or {}).get("H_R2_munu")) or h_r2_default

    componentwise_variation_pack = {
        "index_convention": "mu,nu with covariant metric signature held fixed by strict metric convention lane",
        "H_R2_munu": h_r2,
        "H_Ric2_munu": "H^(Ric2)_{mu nu} = 2 R_{mu alpha}R_nu^alpha - (1/2) g_{mu nu}R_{alpha beta}R^{alpha beta} + Box R_{mu nu} + g_{mu nu} nabla_alpha nabla_beta R^{alpha beta} - 2 nabla_alpha nabla_(mu R_{nu)}^alpha",
        "H_Riem2_munu": "H^(Riem2)_{mu nu} = 2 R_{mu alpha beta gamma}R_nu^{alpha beta gamma} - (1/2) g_{mu nu}R_{alpha beta gamma delta}R^{alpha beta gamma delta} - 4 nabla_alpha nabla_beta R_{mu~~nu}^{~alpha~beta}",
        "H_GB_munu": "H^(GB)_{mu nu} = 2 R R_{mu nu} -4 R_{mu alpha}R_nu^alpha -4 R_{mu alpha nu beta}R^{alpha beta} +2 R_{mu alpha beta gamma}R_nu^{alpha beta gamma} -(1/2) g_{mu nu} G_GB",
        "assembly_rule": "E_{mu nu} = (M_Pl^2/2)G_{mu nu} + c_gr_1 H^(R2)_{mu nu} + c_gr_2 H^(Ric2)_{mu nu} + c_gr_3 H^(Riem2)_{mu nu} + c_gr_4 H^(GB)_{mu nu} - T_{mu nu}^{SM+mix}",
    }

    gravity_operator_profiles_b1 = {
        "profile_id": "gravity_operator_profiles_B1_strict_kernel_jet_v1",
        "background_family": "B1",
        "profile_coordinate": "d",
        "profile_symbols": {
            "K": "K_strict(d)",
            "Kd": "d K_strict(d) / d d",
            "Kdd": "d^2 K_strict(d) / d d^2",
        },
        "source_operator_pack": "gravity_componentwise_variation_pack",
        "profile_generation_rule": (
            "Project the scalar curvature-operator basis O_i=(R2,Ric2,Riem2,GB) "
            "onto the declared one-dimensional strict B1 kernel-jet profile.  "
            "GB is exported by its tensor identity GB=Riem2-4*Ric2+R2, not as an "
            "independent surrogate channel."
        ),
        "profiles": {
            "R2": {
                "operator": "R^2",
                "delta_label": "delta_c_gr_1",
                "expression": "K**2",
            },
            "Ric2": {
                "operator": "Ricci^2",
                "delta_label": "delta_c_gr_2",
                "expression": "Kd**2",
            },
            "Riem2": {
                "operator": "Riemann^2",
                "delta_label": "delta_c_gr_3",
                "expression": "Kdd**2",
            },
            "GB": {
                "operator": "GaussBonnet",
                "delta_label": "delta_c_gr_4",
                "expression": "Riem2 - 4*Ric2 + R2",
            },
        },
        "linearity_warning": (
            "Because GB is the exported tensor combination rather than an "
            "independent proxy, the B1 scalar-profile Gram matrix may be rank "
            "deficient.  A downstream rank obstruction is admissible evidence "
            "and must not be force-promoted to closure."
        ),
    }

    counterterm_basis = p1663.get("counterterm_basis", ["R^2", "Ricci^2", "Riemann^2"])
    counterterm_map = {
        "basis": counterterm_basis,
        "strict_gravity_projection": {
            "R^2": "delta_c_gr_1",
            "Ricci^2": "delta_c_gr_2",
            "Riemann^2": "delta_c_gr_3",
            "GaussBonnet": "delta_c_gr_4",
        },
        "consistency_equations": [
            "c_gr_i^ren(mu) = c_gr_i^bare + delta_c_gr_i(mu)",
            "beta_c_gr_i := d c_gr_i^ren / d ln(mu)",
            "Divergent[Gamma_1loop_grav] + sum_i delta_c_gr_i O_i = finite",
        ],
    }

    residual_trace = {
        "residual_id": "gravity_componentwise_residual_trace_v2",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "checks": [
            {
                "name": "componentwise_variation_exported",
                "status": "PASS_WITH_TRACE",
                "trace": "H^(R2), H^(Ric2), H^(Riem2), H^(GB) exported in strict index convention.",
            },
            {
                "name": "counterterm_basis_projection_exported",
                "status": "PASS_WITH_TRACE",
                "trace": "Counterterm basis projected onto strict gravity coefficient deltas delta_c_gr_i.",
            },
            {
                "name": "theorem_level_counterterm_cancellation_witness",
                "status": "OPEN_OBSTRUCTION_WITH_TRACE",
                "trace": "Need computed divergence tensor and explicit cancellation identity proof for each O_i sector.",
            },
            {
                "name": "unitarity_and_background_independence_joint_witness",
                "status": "OPEN_OBSTRUCTION_WITH_TRACE",
                "trace": "Need joint proof object linking pole residues/optical constraints with nonperturbative diffeomorphism-invariant state family.",
            },
        ],
    }

    out = {
        "packet_id": "P1848",
        "stage_id": "S798",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1847_present": "gravity_sector_exports" in p1847,
            "p1716_present": ("normalized_h_terms" in p1716) or ("input_h_terms" in p1716),
            "p1663_present": "counterterm_basis" in p1663,
        },
        "gravity_componentwise_variation_pack": componentwise_variation_pack,
        "gravity_operator_profiles_B1": gravity_operator_profiles_b1,
        "gravity_counterterm_projection_pack": counterterm_map,
        "gravity_residual_trace": residual_trace,
        "proven": "Componentwise gravity variation structure and counterterm projection map are now explicitly exported in strict-only lane.",
        "open": "Theorem-level divergence cancellation, unitarity, and background-independence witnesses remain open.",
        "false_pass_risk": "Formal componentwise equations alone do not establish renormalizability/unitarity/background-independence closure.",
        "next_honest_step": "Attach computed one-loop divergence tensor coefficients in the same basis and prove explicit cancellation with delta_c_gr_i on one declared background family.",
        "lay_explanation": "Mamy teraz rozpisane dokładne składniki równania grawitacji i jak mają się kasować nieskończoności, ale trzeba jeszcze policzyć to do końca i pokazać pełny dowód.",
    }

    path = GEN / "p1848_s798_strict_gravity_componentwise_variation_and_counterterm_witness_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
