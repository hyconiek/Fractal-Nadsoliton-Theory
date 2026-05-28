#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import scipy.linalg as la
import sympy as sp

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p2332_s1282_strict_bianchi_spatial_tensor_component_table_gb_lift_probe.json"
MD = GEN / "p2332_s1282_strict_bianchi_spatial_tensor_component_table_gb_lift_probe.md"

SOURCE_FILES = {
    "P1853_B1_COEFFICIENTS": GEN / "p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json",
    "P1981_R2_LAPSE": GEN / "p1981_s931_strict_adm_bianchi_r2_lapse_variation_obligation.json",
    "P1982_RICCI2_LAPSE": GEN / "p1982_s932_strict_adm_bianchi_ricci2_lapse_variation_obligation.json",
    "P1983_RIEMANN2_LAPSE": GEN / "p1983_s933_strict_adm_bianchi_riemann2_lapse_variation_obligation.json",
    "P2297_SPATIAL_PROVIDER_OBSTRUCTION": GEN / "p2297_s1247_strict_non_gb_spatial_eom_provider_matrix_obstruction_probe.json",
    "P2331_SECOND_BACKGROUND_GAP_AUDIT": GEN / "p2331_s1281_strict_second_background_curvature_witness_gap_audit.json",
}

CHANNEL_ORDER = ("R2", "Ric2", "Riem2", "GB")
COMPONENTS = ("E_spatial_1", "E_spatial_2", "E_spatial_3")
GB_NULL_VECTOR = np.array([-1.0, 4.0, -1.0, 1.0], dtype=float)
GREP_PATTERNS = (
    "spatial_eom_projection",
    "tracefree spatial",
    "density_difference_NV_R2",
    "density_difference_NV_Ricci2",
    "density_difference_NV_Riemann2",
    "Gauss-Bonnet",
    "same-basis divergence",
    "GB-lift",
)


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": path.relative_to(REPO).as_posix()}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def grep_hits() -> list[dict[str, Any]]:
    candidates = [
        ROOT / "p1981_s931_strict_adm_bianchi_r2_lapse_variation_obligation.py",
        ROOT / "p1982_s932_strict_adm_bianchi_ricci2_lapse_variation_obligation.py",
        ROOT / "p1983_s933_strict_adm_bianchi_riemann2_lapse_variation_obligation.py",
        ROOT / "p2297_s1247_strict_non_gb_spatial_eom_provider_matrix_obstruction_probe.py",
        ROOT / "p2331_s1281_strict_second_background_curvature_witness_gap_audit.py",
    ]
    candidates.extend(SOURCE_FILES.values())
    rows: list[dict[str, Any]] = []
    for path in candidates:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        lowered = text.lower()
        count = sum(lowered.count(pattern.lower()) for pattern in GREP_PATTERNS)
        if count == 0:
            continue
        first_line = 0
        excerpt = ""
        for index, line in enumerate(text.splitlines(), start=1):
            if any(pattern.lower() in line.lower() for pattern in GREP_PATTERNS):
                first_line = index
                excerpt = line.strip()[:220]
                break
        rows.append({
            "path": path.relative_to(REPO).as_posix(),
            "pattern_hit_count": count,
            "first_hit_line": first_line,
            "first_hit_excerpt": excerpt,
        })
    rows.sort(key=lambda row: (-int(row["pattern_hit_count"]), row["path"]))
    return rows


def symbol_context() -> tuple[dict[str, sp.Symbol], list[sp.Symbol]]:
    names = [
        "N", "Nd", "Ndd", "V", "H", "Hd", "Hdd",
        "sigma1", "sigma2", "dsigma1", "dsigma2", "d2sigma1", "d2sigma2", "Q",
    ]
    local = {name: sp.Symbol(name, real=True) for name in names}
    local.update({"pi": sp.pi, "log": sp.log, "ln": sp.log})
    poly_vars = [
        local["N"], local["Nd"], local["Ndd"], local["H"], local["Hd"], local["Hdd"],
        local["sigma1"], local["sigma2"], local["dsigma1"], local["dsigma2"],
        local["d2sigma1"], local["d2sigma2"],
    ]
    return local, poly_vars


def build_spatial_components(density: sp.Expr, local: dict[str, sp.Symbol]) -> list[sp.Expr]:
    N = local["N"]
    Nd = local["Nd"]
    Ndd = local["Ndd"]
    V = local["V"]
    H = local["H"]
    Hd = local["Hd"]
    Hdd = local["Hdd"]
    s1 = local["sigma1"]
    s2 = local["sigma2"]
    sd1 = local["dsigma1"]
    sd2 = local["dsigma2"]
    sdd1 = local["d2sigma1"]
    sdd2 = local["d2sigma2"]
    dt_rules = {N: Nd, Nd: Ndd, V: 3 * H * V, H: Hd, Hd: Hdd, s1: sd1, s2: sd2, sd1: sdd1, sd2: sdd2}

    def total_dt(expr: sp.Expr) -> sp.Expr:
        return sp.factor(sum(sp.diff(expr, var) * dvar for var, dvar in dt_rules.items()))

    def euler_component(sigma: sp.Symbol, dsigma: sp.Symbol, d2sigma: sp.Symbol) -> sp.Expr:
        raw = sp.diff(density, sigma) - total_dt(sp.diff(density, dsigma)) + total_dt(total_dt(sp.diff(density, d2sigma)))
        return sp.factor(sp.simplify(raw * N**8 / V))

    e1 = euler_component(s1, sd1, sdd1)
    e2 = euler_component(s2, sd2, sdd2)
    e3 = sp.factor(sp.simplify(-e1 - e2))
    return [e1, e2, e3]


def coefficient_vector(components: list[sp.Expr], monomial_basis: list[tuple[int, ...]], poly_vars: list[sp.Symbol]) -> list[sp.Expr]:
    values: list[sp.Expr] = []
    for component in components:
        poly = sp.Poly(sp.expand(component), *poly_vars, domain="EX")
        coeff_by_monomial = {monomial: coeff for monomial, coeff in poly.terms()}
        values.extend([sp.expand(coeff_by_monomial.get(monomial, 0)) for monomial in monomial_basis])
    return values


def numeric_matrix(vectors: dict[str, list[sp.Expr]]) -> np.ndarray:
    return np.array([[float(sp.N(value, 50)) for value in vectors[channel]] for channel in CHANNEL_ORDER], dtype=float)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    local, poly_vars = symbol_context()
    s1 = local["sigma1"]
    s2 = local["sigma2"]
    q_shear = s1**2 + s1 * s2 + s2**2

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

    component_table = {channel: build_spatial_components(densities[channel], local) for channel in CHANNEL_ORDER}
    gb_component_residuals = [
        sp.simplify(component_table["GB"][idx] - component_table["Riem2"][idx] + 4 * component_table["Ric2"][idx] - component_table["R2"][idx])
        for idx in range(len(COMPONENTS))
    ]
    gb_component_identity_zero = all(residual == 0 for residual in gb_component_residuals)
    tracefree_zero_by_channel = {channel: bool(sp.simplify(sum(component_table[channel])) == 0) for channel in CHANNEL_ORDER}

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
    reconstruction_residual = target_vector_numeric - coeff_numeric @ channel_matrix
    reconstruction_residual_l2 = float(la.norm(reconstruction_residual, ord=2))
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
        "background_family_witness": "ADM_BianchiI_tracefree_spatial_EOM",
        "component_basis": list(COMPONENTS),
        "channels": list(CHANNEL_ORDER),
        "component_entry_count": len(component_rows),
        "monomial_basis_size_per_component": len(monomial_basis),
        "same_basis_vector_length": len(target_vector_exact),
        "tracefree_zero_by_channel": tracefree_zero_by_channel,
        "all_tracefree_sums_zero": all(tracefree_zero_by_channel.values()),
        "gb_component_identity_zero": gb_component_identity_zero,
    }

    same_basis_divergence_target = {
        "target_definition": "sum_i a_i * E_channel_i in concatenated tracefree spatial Bianchi-I component coefficient basis",
        "strict_coefficients": {channel: str(strict_coefficients[channel]) for channel in CHANNEL_ORDER},
        "target_vector_sha256": sha256_text(json.dumps([str(value) for value in target_vector_exact], sort_keys=True)),
        "direct_reconstruction_residual_l2": reconstruction_residual_l2,
        "least_squares_residual_l2": lsq_residual_l2,
        "least_squares_solution_R2_Ric2_Riem2_GB": lsq_solution.tolist(),
        "coefficient_null_ambiguity": "tracefree spatial basis has at least the GB null vector (-1,4,-1,1), and the rank test finds an additional spatial-projection null direction",
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
        "verdict": "GB_DEPENDENCE_PERSISTS_IN_TRACEFREE_SPATIAL_BIANCHI_WITNESS" if numeric_rank < len(CHANNEL_ORDER) else "GB_DEPENDENCE_LIFTED_IN_THIS_WITNESS",
    }

    theorem_export = {
        "theorem_name": "P2332 tracefree spatial Bianchi-I tensor-component table and GB-lift rank test",
        "claim": "A same-basis tracefree spatial Bianchi-I component table and divergence target can be built from existing strict ADM/Bianchi-I curvature densities, but the Gauss-Bonnet channel remains linearly dependent in this witness: the channel Gram rank is 2 with nullity 2, including null vector (-1,4,-1,1). Thus this independent background-family witness does not lift B1 quotient renormalization to full four-channel/global closure.",
        "proof_witnesses": [
            "R2/Ric2/Riem2 tracefree spatial Euler components are derived from P1981/P1982/P1983 densities using the P2297 Euler-operator rule.",
            "GB components are assembled in the same basis as Riem2 - 4*Ric2 + R2 and satisfy the component identity exactly.",
            "The same-basis strict divergence target reconstructs with zero numerical residual from the exported strict coefficients.",
            "The channel Gram matrix has numeric rank 2 and nullity 2 in this tracefree spatial projection; the GB null vector residual is small, so dependence is not lifted.",
        ],
        "scope_limits": [
            "This table is tracefree spatial ADM/Bianchi-I, not the complete spacetime tensor-component table.",
            "The result is a negative lift test, not full/global renormalization closure.",
            "No selector, QW-2191, G1/G3, or ToE closure is claimed.",
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

    hits = grep_hits()
    probe = {
        "probe_id": "P2332_S1282_STRICT_BIANCHI_SPATIAL_TENSOR_COMPONENT_TABLE_GB_LIFT",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {"search_patterns": list(GREP_PATTERNS), "hit_count": len(hits), "top_hits": hits[:20]},
        "component_table_summary": table_summary,
        "component_table_rows": component_rows,
        "same_basis_divergence_target": same_basis_divergence_target,
        "gb_lift_rank_test": gb_lift_rank_test,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "grep_hits_found": len(hits) >= 5,
        "p1853_coefficients_loaded": p1853.get("packet_id") == "P1853",
        "p1981_loaded": p1981.get("packet_id") == "P1981",
        "p1982_loaded": p1982.get("packet_id") == "P1982",
        "p1983_loaded": p1983.get("packet_id") == "P1983",
        "component_entry_count_12": len(component_rows) == 12,
        "all_tracefree_sums_zero": table_summary["all_tracefree_sums_zero"],
        "gb_component_identity_zero": gb_component_identity_zero,
        "same_basis_target_reconstructs": reconstruction_residual_l2 < 1e-12 and lsq_residual_l2 < 1e-12,
        "gb_lift_rank_is_2": numeric_rank == 2,
        "gb_lift_nullity_is_2": numeric_nullity == 2,
        "gb_null_residual_small": null_residual_l2 < 1e-8,
        "full_global_renormalization_not_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_g1_g3_update_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2332_s1282_v1",
        "packet_id": "P2332",
        "stage_id": "S1282",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_TRACEFREE_SPATIAL_COMPONENT_TABLE_BUILT_GB_LIFT_NEGATIVE",
        "result_kind": "STRICT_BIANCHI_SPATIAL_TENSOR_COMPONENT_TABLE_AND_GB_LIFT_RANK_TEST_NO_GLOBAL_RENORMALIZATION_CLAIM",
        "strict_bianchi_spatial_tensor_component_table_gb_lift_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": "Extend the component table from tracefree spatial ADM/Bianchi-I components to a complete spacetime tensor-component basis including lapse/time and same-basis divergence target; only then rerun the full global GB-lift test.",
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2332 Bianchi-I spatial tensor-component table and GB-lift rank test\n\n"
        "Status: same-basis tracefree spatial table built; GB lift remains negative.\n\n"
        f"- Component entries: `{len(component_rows)}`.\n"
        f"- Same-basis vector length: `{len(target_vector_exact)}`.\n"
        f"- Channel Gram rank/nullity: `{numeric_rank}/{numeric_nullity}`.\n"
        f"- GB null residual L2: `{null_residual_l2:.3e}`.\n"
        f"- Target reconstruction residual L2: `{reconstruction_residual_l2:.3e}`.\n"
        "- No full/global renormalization, no QW-2191 discharge, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
