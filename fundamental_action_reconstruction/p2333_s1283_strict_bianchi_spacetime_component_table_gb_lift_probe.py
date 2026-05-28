#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import scipy.linalg as la
import sympy as sp

from p2332_s1282_strict_bianchi_spatial_tensor_component_table_gb_lift_probe import (
    CHANNEL_ORDER,
    GB_NULL_VECTOR,
    build_spatial_components,
    coefficient_vector,
    grep_hits as inherited_grep_hits,
    load_json,
    numeric_matrix,
    sha256_file,
    sha256_json,
    sha256_text,
    symbol_context,
)

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2333_s1283_strict_bianchi_spacetime_component_table_gb_lift_probe.json"
MD = GEN / "p2333_s1283_strict_bianchi_spacetime_component_table_gb_lift_probe.md"

SOURCE_FILES = {
    "P1853_B1_COEFFICIENTS": GEN / "p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json",
    "P1981_R2_LAPSE": GEN / "p1981_s931_strict_adm_bianchi_r2_lapse_variation_obligation.json",
    "P1982_RICCI2_LAPSE": GEN / "p1982_s932_strict_adm_bianchi_ricci2_lapse_variation_obligation.json",
    "P1983_RIEMANN2_LAPSE": GEN / "p1983_s933_strict_adm_bianchi_riemann2_lapse_variation_obligation.json",
    "P1984_GB_LAPSE": GEN / "p1984_s934_strict_adm_bianchi_gauss_bonnet_lapse_cancellation_witness.json",
    "P2332_SPATIAL_TABLE": GEN / "p2332_s1282_strict_bianchi_spatial_tensor_component_table_gb_lift_probe.json",
}
COMPONENTS = ("E_lapse_N", "E_spatial_1", "E_spatial_2", "E_spatial_3")


def parse_lapse_components(artifacts: dict[str, dict[str, Any]], local: dict[str, sp.Symbol]) -> dict[str, sp.Expr]:
    q_shear = local["sigma1"] ** 2 + local["sigma1"] * local["sigma2"] + local["sigma2"] ** 2
    qd_shear = 2 * local["sigma1"] * local["dsigma1"] + local["sigma1"] * local["dsigma2"] + local["sigma2"] * local["dsigma1"] + 2 * local["sigma2"] * local["dsigma2"]
    q_subs = {local["Q"]: q_shear, sp.Symbol("Qd"): qd_shear}
    e_r2 = sp.sympify(
        artifacts["P1981_R2_LAPSE"].get("r2_lapse_euler_operator", {}).get("EL_N_difference", "0"),
        locals=local,
    ).subs(q_subs)
    e_ric2 = sp.sympify(
        artifacts["P1982_RICCI2_LAPSE"].get("ricci2_lapse_euler_operator", {}).get("EL_N_difference", "0"),
        locals=local,
    ).subs(q_subs)
    e_riem2 = sp.sympify(
        artifacts["P1983_RIEMANN2_LAPSE"].get("riemann2_lapse_euler_operator", {}).get("EL_N_difference", "0"),
        locals=local,
    ).subs(q_subs)
    scale = local["N"] ** 8 / local["V"]
    e_gb = sp.factor(e_riem2 - 4 * e_ric2 + e_r2)
    return {"R2": sp.factor(e_r2 * scale), "Ric2": sp.factor(e_ric2 * scale), "Riem2": sp.factor(e_riem2 * scale), "GB": sp.factor(e_gb * scale)}


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    local, poly_vars = symbol_context()
    q_shear = local["sigma1"] ** 2 + local["sigma1"] * local["sigma2"] + local["sigma2"] ** 2

    p1853 = artifacts["P1853_B1_COEFFICIENTS"]
    p1981 = artifacts["P1981_R2_LAPSE"]
    p1982 = artifacts["P1982_RICCI2_LAPSE"]
    p1983 = artifacts["P1983_RIEMANN2_LAPSE"]
    coeffs = p1853.get("b1_symbolic_evaluation", {}).get("evaluated_coefficients", {})
    strict_coefficients = {
        "R2": sp.sympify(coeffs.get("a_R2", {}).get("symbolic", "0"), locals=local),
        "Ric2": sp.sympify(coeffs.get("a_Ric2", {}).get("symbolic", "0"), locals=local),
        "Riem2": sp.sympify(coeffs.get("a_Riem2", {}).get("symbolic", "0"), locals=local),
        "GB": sp.sympify(coeffs.get("a_GB", {}).get("symbolic", "0"), locals=local),
    }

    densities = {
        "R2": sp.sympify(p1981.get("r2_lapse_euler_operator", {}).get("density_difference_NV_R2", "0"), locals=local).subs({local["Q"]: q_shear}),
        "Ric2": sp.sympify(p1982.get("ricci2_lapse_euler_operator", {}).get("density_difference_NV_Ricci2", "0"), locals=local),
        "Riem2": sp.sympify(p1983.get("riemann2_lapse_euler_operator", {}).get("density_difference_NV_Riemann2", "0"), locals=local),
    }
    densities["GB"] = sp.factor(densities["Riem2"] - 4 * densities["Ric2"] + densities["R2"])

    lapse_components = parse_lapse_components(artifacts, local)
    spatial_components = {channel: build_spatial_components(densities[channel], local) for channel in CHANNEL_ORDER}
    component_table = {channel: [lapse_components[channel], *spatial_components[channel]] for channel in CHANNEL_ORDER}
    gb_component_residuals = [
        sp.simplify(component_table["GB"][idx] - component_table["Riem2"][idx] + 4 * component_table["Ric2"][idx] - component_table["R2"][idx])
        for idx in range(len(COMPONENTS))
    ]
    gb_component_identity_zero = all(residual == 0 for residual in gb_component_residuals)
    tracefree_spatial_zero_by_channel = {
        channel: bool(sp.simplify(sum(component_table[channel][1:])) == 0) for channel in CHANNEL_ORDER
    }

    monomial_set: set[tuple[int, ...]] = set()
    for components in component_table.values():
        for component in components:
            monomial_set.update(monomial for monomial, _coeff in sp.Poly(sp.expand(component), *poly_vars, domain="EX").terms())
    monomial_basis = sorted(monomial_set)
    coefficient_vectors = {
        channel: coefficient_vector(component_table[channel], monomial_basis, poly_vars) for channel in CHANNEL_ORDER
    }
    channel_matrix = numeric_matrix(coefficient_vectors)
    gram = channel_matrix @ channel_matrix.T
    singular_values = np.linalg.svd(gram, compute_uv=False)
    max_sv = float(np.max(singular_values)) if singular_values.size else 0.0
    rank_tolerance = max(gram.shape) * max_sv * np.finfo(float).eps * 100.0 if max_sv else 0.0
    numeric_rank = int(np.sum(singular_values > rank_tolerance)) if max_sv else 0
    numeric_nullity = int(gram.shape[1] - numeric_rank)
    null_residual = gram @ GB_NULL_VECTOR
    null_residual_l2 = float(la.norm(null_residual, ord=2))

    target_vector_exact = [
        sp.expand(sum(strict_coefficients[channel] * coefficient_vectors[channel][idx] for channel in CHANNEL_ORDER))
        for idx in range(len(coefficient_vectors["R2"]))
    ]
    target_vector_numeric = np.array([float(sp.N(value, 50)) for value in target_vector_exact], dtype=float)
    coeff_numeric = np.array([float(sp.N(strict_coefficients[channel], 50)) for channel in CHANNEL_ORDER], dtype=float)
    reconstruction_residual_l2 = float(la.norm(target_vector_numeric - coeff_numeric @ channel_matrix, ord=2))
    lsq_solution, *_ = np.linalg.lstsq(channel_matrix.T, target_vector_numeric, rcond=None)
    lsq_residual_l2 = float(la.norm(channel_matrix.T @ lsq_solution - target_vector_numeric, ord=2))

    component_rows = []
    for channel in CHANNEL_ORDER:
        for component_name, expression in zip(COMPONENTS, component_table[channel]):
            component_rows.append({
                "channel": channel,
                "component": component_name,
                "expression_sha256": sha256_text(str(expression)),
                "expression_excerpt": str(expression)[:260],
                "term_count": len(sp.Poly(sp.expand(expression), *poly_vars, domain="EX").terms()),
            })

    table_summary = {
        "background_family_witness": "ADM_BianchiI_lapse_plus_tracefree_spatial_EOM",
        "component_basis": list(COMPONENTS),
        "channels": list(CHANNEL_ORDER),
        "component_entry_count": len(component_rows),
        "monomial_basis_size_per_component": len(monomial_basis),
        "same_basis_vector_length": len(target_vector_exact),
        "tracefree_spatial_zero_by_channel": tracefree_spatial_zero_by_channel,
        "all_spatial_tracefree_sums_zero": all(tracefree_spatial_zero_by_channel.values()),
        "gb_component_identity_zero": gb_component_identity_zero,
    }
    same_basis_divergence_target = {
        "target_definition": "sum_i a_i * E_channel_i in concatenated lapse+tracefree-spatial Bianchi-I coefficient basis",
        "strict_coefficients": {channel: str(strict_coefficients[channel]) for channel in CHANNEL_ORDER},
        "target_vector_sha256": sha256_text(json.dumps([str(value) for value in target_vector_exact], sort_keys=True)),
        "direct_reconstruction_residual_l2": reconstruction_residual_l2,
        "least_squares_residual_l2": lsq_residual_l2,
        "least_squares_solution_R2_Ric2_Riem2_GB": lsq_solution.tolist(),
        "coefficient_null_ambiguity": "GB null vector (-1,4,-1,1) remains exact in lapse+tracefree-spatial Bianchi-I basis",
    }
    gb_lift_rank_test = {
        "channel_gram_matrix": gram.tolist(),
        "singular_values": singular_values.tolist(),
        "rank_tolerance": rank_tolerance,
        "numeric_rank": numeric_rank,
        "numeric_nullity": numeric_nullity,
        "null_vector_R2_Ric2_Riem2_GB": GB_NULL_VECTOR.tolist(),
        "null_vector_residual_l2": null_residual_l2,
        "gb_dependence_lifted": numeric_rank == len(CHANNEL_ORDER),
        "verdict": "GB_DEPENDENCE_PERSISTS_IN_SPACETIME_BIANCHI_WITNESS" if numeric_rank < len(CHANNEL_ORDER) else "GB_DEPENDENCE_LIFTED_IN_THIS_WITNESS",
    }
    selection_omission_note = {
        "lay_metaphor": "Choosing one valley is simultaneously not-choosing the other equivalent valleys.",
        "formal_analogue": "A signed selector would define both a selected branch and a complementary excluded branch set within the degenerate orbit.",
        "current_status": "Conceptually useful for future QW-2191 work, but not used in this renormalization proof and not a selector discharge.",
    }

    theorem_export = {
        "theorem_name": "P2333 lapse plus tracefree spatial Bianchi-I component table and full local GB-lift test",
        "claim": "Adding the lapse/time component to the tracefree spatial ADM/Bianchi-I component table yields a same-basis 16-entry local spacetime witness and divergence target, but the Gauss-Bonnet channel remains linearly dependent: the channel Gram rank is 2 with nullity 2 and includes null vector (-1,4,-1,1). Therefore this local spacetime Bianchi-I witness still does not promote quotient-scope renormalization to full/global closure.",
        "proof_witnesses": [
            "E_lapse_N components are loaded from P1981/P1982/P1983 and GB is assembled in the same basis.",
            "Tracefree spatial components are derived by the P2297 Euler rule and combined with lapse components.",
            "The same-basis divergence target reconstructs from strict coefficients with small numerical residual.",
            "The channel Gram matrix rank is 2, nullity is 2, and the GB null vector residual is small.",
        ],
        "scope_limits": [
            "This is ADM/Bianchi-I local spacetime component evidence, not a proof over all backgrounds.",
            "A global atlas/transport theorem is still required before full renormalization closure.",
            "The selection-omission valley intuition is recorded only as an interpretation note; no QW-2191 selector is exported.",
        ],
        "strict_guardrails": {
            "strict_kernel_only": True,
            "no_legacy_kernel_role_transfer": True,
            "no_selector_premise_added": True,
            "no_qw2191_discharge_claimed": True,
            "no_g1_g3_update_claimed": True,
            "no_toe_closure_claimed": True,
        },
    }

    hits = inherited_grep_hits()
    probe = {
        "probe_id": "P2333_S1283_STRICT_BIANCHI_SPACETIME_COMPONENT_TABLE_GB_LIFT",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {"search_patterns": "inherited_from_P2332", "hit_count": len(hits), "top_hits": hits[:20]},
        "component_table_summary": table_summary,
        "component_table_rows": component_rows,
        "same_basis_divergence_target": same_basis_divergence_target,
        "gb_lift_rank_test": gb_lift_rank_test,
        "selection_omission_coupling_note": selection_omission_note,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }
    gatekeeper_checks = {
        "p1853_coefficients_loaded": p1853.get("packet_id") == "P1853",
        "p1981_loaded": p1981.get("packet_id") == "P1981",
        "p1982_loaded": p1982.get("packet_id") == "P1982",
        "p1983_loaded": p1983.get("packet_id") == "P1983",
        "component_entry_count_16": len(component_rows) == 16,
        "all_spatial_tracefree_sums_zero": table_summary["all_spatial_tracefree_sums_zero"],
        "gb_component_identity_zero": gb_component_identity_zero,
        "same_basis_target_reconstructs": reconstruction_residual_l2 < 1e-10 and lsq_residual_l2 < 1e-10,
        "gb_lift_rank_is_2": numeric_rank == 2,
        "gb_lift_nullity_is_2": numeric_nullity == 2,
        "gb_null_residual_small": null_residual_l2 < 1e-6,
        "full_global_renormalization_not_claimed": True,
        "selection_omission_note_not_promoted_to_selector": True,
        "no_legacy_kernel_role_transfer": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_g1_g3_update_claimed": True,
        "no_toe_closure_claimed": True,
    }
    payload = {
        "schema_version": "p2333_s1283_v1",
        "packet_id": "P2333",
        "stage_id": "S1283",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_BIANCHI_SPACETIME_COMPONENT_TABLE_BUILT_GB_LIFT_STILL_NEGATIVE",
        "result_kind": "STRICT_BIANCHI_SPACETIME_COMPONENT_TABLE_AND_GB_LIFT_RANK_TEST_NO_GLOBAL_RENORMALIZATION_CLAIM",
        "strict_bianchi_spacetime_component_table_gb_lift_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": "Add a genuine second atlas/background transport theorem beyond local ADM/Bianchi-I, or prove that GB topological dependence remains quotient-only across the admitted strict background atlas.",
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2333 Bianchi-I spacetime component table and GB-lift rank test\n\n"
        "Status: lapse+tracefree-spatial table built; local GB lift remains negative.\n\n"
        f"- Component entries: `{len(component_rows)}`.\n"
        f"- Same-basis vector length: `{len(target_vector_exact)}`.\n"
        f"- Channel Gram rank/nullity: `{numeric_rank}/{numeric_nullity}`.\n"
        f"- GB null residual L2: `{null_residual_l2:.3e}`.\n"
        f"- Target reconstruction residual L2: `{reconstruction_residual_l2:.3e}`.\n"
        "- Valley intuition: choosing one branch also excludes the other degenerate branches, but this is only an interpretation note here.\n"
        "- No full/global renormalization, no QW-2191 discharge, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
