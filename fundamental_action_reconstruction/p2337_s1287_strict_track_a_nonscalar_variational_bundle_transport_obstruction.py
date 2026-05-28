#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

from p2332_s1282_strict_bianchi_spatial_tensor_component_table_gb_lift_probe import (
    build_spatial_components,
    load_json,
    sha256_file,
    symbol_context,
)
from p2333_s1283_strict_bianchi_spacetime_component_table_gb_lift_probe import parse_lapse_components

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2337_s1287_strict_track_a_nonscalar_variational_bundle_transport_obstruction.json"
MD = GEN / "p2337_s1287_strict_track_a_nonscalar_variational_bundle_transport_obstruction.md"

SOURCE_FILES = {
    "P1981_R2_LAPSE": GEN / "p1981_s931_strict_adm_bianchi_r2_lapse_variation_obligation.json",
    "P1982_RICCI2_LAPSE": GEN / "p1982_s932_strict_adm_bianchi_ricci2_lapse_variation_obligation.json",
    "P1983_RIEMANN2_LAPSE": GEN / "p1983_s933_strict_adm_bianchi_riemann2_lapse_variation_obligation.json",
    "P2333_BIANCHI_SPACETIME_TABLE": GEN / "p2333_s1283_strict_bianchi_spacetime_component_table_gb_lift_probe.json",
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2336_TRACK_A_SCALAR_TRANSPORT": GEN / "p2336_s1286_strict_track_a_second_atlas_finite_part_transport_theorem.json",
}
COMPONENTS = ("E_lapse_N", "E_spatial_1", "E_spatial_2", "E_spatial_3")
TRACK_A_CHANNELS = ("R2", "Ric2")
ANISOTROPY_ZERO_NAMES = ("sigma1", "sigma2", "dsigma1", "dsigma2", "d2sigma1", "d2sigma2")


def sha256_json_local(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def grep_hits() -> list[str]:
    candidates = [
        ROOT / "p2333_s1283_strict_bianchi_spacetime_component_table_gb_lift_probe.py",
        ROOT / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.py",
        ROOT / "p2336_s1286_strict_track_a_second_atlas_finite_part_transport_theorem.py",
        ROOT / "p2337_s1287_strict_track_a_nonscalar_variational_bundle_transport_obstruction.py",
        GEN / "p2333_s1283_strict_bianchi_spacetime_component_table_gb_lift_probe.json",
        GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
        GEN / "p2336_s1286_strict_track_a_second_atlas_finite_part_transport_theorem.json",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "P2336|P2335|P2333|Track A|non-scalar|nonscalar|variational-bundle|tensor-component|FRW restriction|Bianchi|scalar lift|GB ledger|QW-2191",
        *existing,
    ]
    proc = subprocess.run(cmd, cwd=ROOT.parent, text=True, capture_output=True, check=False)
    return proc.stdout.splitlines()[:120]


def track_a_coefficients(p2335: dict[str, Any], local: dict[str, sp.Symbol]) -> tuple[sp.Expr, sp.Expr]:
    probe = p2335.get("strict_two_track_renormalization_ledger_probe", {})
    ledger = probe.get("two_track_ledger", {})
    track_a = ledger.get("track_A_non_gb_eom_transportable_quotient", {})
    coeffs = track_a.get("coefficients", {})
    return sp.sympify(coeffs.get("b_EOM_R2", "0"), locals=local), sp.sympify(
        coeffs.get("b_EOM_Ric2", "0"), locals=local
    )


def build_track_a_components(
    artifacts: dict[str, dict[str, Any]], local: dict[str, sp.Symbol]
) -> tuple[list[sp.Expr], dict[str, list[sp.Expr]], tuple[sp.Expr, sp.Expr]]:
    q_shear = local["sigma1"] ** 2 + local["sigma1"] * local["sigma2"] + local["sigma2"] ** 2
    p1981 = artifacts["P1981_R2_LAPSE"]
    p1982 = artifacts["P1982_RICCI2_LAPSE"]
    p2335 = artifacts["P2335_TWO_TRACK_LEDGER"]

    densities = {
        "R2": sp.sympify(
            p1981.get("r2_lapse_euler_operator", {}).get("density_difference_NV_R2", "0"), locals=local
        ).subs({local["Q"]: q_shear}),
        "Ric2": sp.sympify(
            p1982.get("ricci2_lapse_euler_operator", {}).get("density_difference_NV_Ricci2", "0"), locals=local
        ),
    }
    lapse_components = parse_lapse_components(artifacts, local)
    channel_components = {
        channel: [lapse_components[channel], *build_spatial_components(densities[channel], local)]
        for channel in TRACK_A_CHANNELS
    }
    b_r2, b_ric2 = track_a_coefficients(p2335, local)
    track_a_components = [
        sp.factor(sp.simplify(b_r2 * channel_components["R2"][idx] + b_ric2 * channel_components["Ric2"][idx]))
        for idx in range(len(COMPONENTS))
    ]
    return track_a_components, channel_components, (b_r2, b_ric2)


def component_term_counts(components: list[sp.Expr], poly_vars: list[sp.Symbol]) -> list[int]:
    return [len(sp.Poly(sp.expand(component), *poly_vars, domain="EX").terms()) for component in components]


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    local, poly_vars = symbol_context()
    track_a_components, _channel_components, (b_r2, b_ric2) = build_track_a_components(artifacts, local)

    frw_zero_subs = {local[name]: sp.Integer(0) for name in ANISOTROPY_ZERO_NAMES}
    frw_restriction = [sp.factor(sp.simplify(component.subs(frw_zero_subs))) for component in track_a_components]
    scalar_lift_prediction_from_frw = [sp.Integer(0) for _component in track_a_components]
    scalar_lift_residual = [
        sp.factor(sp.simplify(track_a_components[idx] - scalar_lift_prediction_from_frw[idx]))
        for idx in range(len(track_a_components))
    ]
    residual_term_counts = component_term_counts(scalar_lift_residual, poly_vars)
    residual_nonzero_by_component = [count > 0 for count in residual_term_counts]
    residual_symbolically_nonzero = any(residual_nonzero_by_component)
    frw_restriction_zero = all(value == 0 for value in frw_restriction)
    tracefree_spatial_sum = sp.factor(sp.simplify(sum(track_a_components[1:])))

    sample_subs = {
        local["N"]: sp.Rational(1, 1),
        local["Nd"]: sp.Rational(3, 100),
        local["Ndd"]: sp.Rational(-1, 50),
        local["H"]: sp.Rational(7, 10),
        local["Hd"]: sp.Rational(11, 100),
        local["Hdd"]: sp.Rational(-3, 100),
        local["sigma1"]: sp.Rational(1, 5),
        local["sigma2"]: sp.Rational(-1, 10),
        local["dsigma1"]: sp.Rational(1, 20),
        local["dsigma2"]: sp.Rational(-1, 50),
        local["d2sigma1"]: sp.Rational(1, 100),
        local["d2sigma2"]: sp.Rational(-3, 200),
    }
    sample_values = [float(sp.N(component.subs(sample_subs), 40)) for component in scalar_lift_residual]
    sample_l2 = float(np.linalg.norm(np.array(sample_values, dtype=float), ord=2))

    p2333 = artifacts["P2333_BIANCHI_SPACETIME_TABLE"]
    p2335 = artifacts["P2335_TWO_TRACK_LEDGER"]
    p2336 = artifacts["P2336_TRACK_A_SCALAR_TRANSPORT"]
    p2336_checks = p2336.get("gatekeeper_checks", {})
    dependencies = {
        "p2333_component_table_loaded": p2333.get("packet_id") == "P2333"
        and (p2333.get("gatekeeper_checks") or {}).get("component_entry_count_16") is True,
        "p2335_two_track_ledger_loaded": p2335.get("result_kind")
        == "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
        "p2336_track_a_scalar_transport_loaded": p2336.get("result_kind")
        == "PASS_TRACK_A_LOCAL_SECOND_ATLAS_FINITE_PART_TRANSPORT__NO_TRACK_B_OR_GLOBAL_CLOSURE",
        "p2336_full_tensor_bundle_not_claimed": p2336_checks.get("no_full_tensor_bundle_transport_claimed") is True,
        "p2336_track_b_not_transported": p2336_checks.get("track_b_not_transported") is True,
    }

    theorem_export = {
        "theorem_name": "P2337 Track-A non-scalar variational-bundle scalar-lift obstruction",
        "claim": (
            "The P2336 scalar finite-part transport of the two Track-A coefficients cannot be silently promoted to a "
            "non-scalar variational-bundle/tensor-component transport theorem.  In the P2333 ADM/Bianchi-I component "
            "basis the Track-A anisotropic component vector vanishes on the FRW restriction, so any pure scalar lift from "
            "that FRW anchor predicts zero components; the actual Bianchi-I Track-A component residual is symbolically nonzero."
        ),
        "proof_witnesses": [
            "Track-A components are rebuilt from P1981/P1982 lapse data and the P2333/P2297 spatial Euler rule.",
            "The FRW anisotropy-zero restriction of the Track-A Bianchi difference vector is exactly zero.",
            "The same Track-A vector has nonzero symbolic residual terms away from the FRW restriction.",
            "A rational anisotropic sample gives a nonzero residual L2 norm.",
        ],
        "not_licensed": [
            "full non-scalar variational-bundle transport theorem",
            "global background-independence theorem",
            "transport or normalization of the Track B GB topological ledger",
            "independent a_GB measurement",
            "boundary/topological-number normalization theorem",
            "full/global renormalization closure",
            "QW-2191 selector discharge",
            "G1/G3 update",
            "ToE closure",
        ],
        "next_honest_step": (
            "Do not reuse the scalar finite-part multiplier as a tensor-bundle theorem.  Either construct a genuinely "
            "non-scalar connection/provider that maps the Track-A component residuals, or switch to the separate Track-B "
            "boundary/topological-number normalization witness."
        ),
    }

    component_rows = []
    for idx, component_name in enumerate(COMPONENTS):
        residual = scalar_lift_residual[idx]
        component_rows.append(
            {
                "component": component_name,
                "frw_restriction": str(frw_restriction[idx]),
                "scalar_lift_prediction": str(scalar_lift_prediction_from_frw[idx]),
                "residual_term_count": residual_term_counts[idx],
                "residual_nonzero": residual_nonzero_by_component[idx],
                "residual_sha256": hashlib.sha256(str(sp.factor(residual)).encode("utf-8")).hexdigest(),
                "residual_excerpt": str(residual)[:260],
            }
        )

    probe = {
        "probe_id": "P2337_S1287_STRICT_TRACK_A_NONSCALAR_VARIATIONAL_BUNDLE_TRANSPORT_OBSTRUCTION",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2336 scalar Track-A transport versus P2333 tensor-component/non-scalar upgrade boundary",
            "top_hits": grep_hits(),
        },
        "track_A_coefficients": {"b_EOM_R2": str(b_r2), "b_EOM_Ric2": str(b_ric2)},
        "component_basis": list(COMPONENTS),
        "scalar_lift_obstruction": {
            "frw_restriction_rule": {name: "0" for name in ANISOTROPY_ZERO_NAMES},
            "frw_restriction_zero": frw_restriction_zero,
            "scalar_lift_model_tested": "pure P2336 scalar multiplier applied to the FRW Track-A component anchor",
            "scalar_lift_prediction_reason": "FRW Track-A anisotropic difference anchor is zero, so a scalar multiplier predicts zero component residuals.",
            "residual_symbolically_nonzero": residual_symbolically_nonzero,
            "residual_term_counts": residual_term_counts,
            "tracefree_spatial_sum": str(tracefree_spatial_sum),
            "tracefree_spatial_sum_zero": tracefree_spatial_sum == 0,
            "sample_values": sample_values,
            "sample_l2_norm": sample_l2,
            "component_rows": component_rows,
        },
        "current_export_dependencies": dependencies,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json_local(theorem_export),
    }

    gatekeeper_checks = {
        "repo_grep_audit_recorded": len(probe["repo_grep_audit"]["top_hits"]) > 0,
        "p2333_component_table_loaded": dependencies["p2333_component_table_loaded"],
        "p2335_two_track_ledger_loaded": dependencies["p2335_two_track_ledger_loaded"],
        "p2336_track_a_scalar_transport_loaded": dependencies["p2336_track_a_scalar_transport_loaded"],
        "p2336_full_tensor_bundle_limit_preserved": dependencies["p2336_full_tensor_bundle_not_claimed"],
        "p2336_track_b_boundary_preserved": dependencies["p2336_track_b_not_transported"],
        "frw_restriction_zero": frw_restriction_zero,
        "scalar_lift_residual_symbolically_nonzero": residual_symbolically_nonzero,
        "sample_l2_norm_nonzero": sample_l2 > 1e-6,
        "tracefree_spatial_sum_zero": tracefree_spatial_sum == 0,
        "no_scalar_to_tensor_bundle_promotion_claimed": True,
        "no_track_b_transport_claimed": True,
        "no_independent_a_gb_claimed": True,
        "no_boundary_topological_number_witness_claimed": True,
        "no_background_globalization_claimed": True,
        "no_full_global_renormalization_claimed": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_g1_g3_update_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2337_s1287_v1",
        "packet_id": "P2337",
        "stage_id": "S1287",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "result_kind": "STRICT_TRACK_A_NONSCALAR_VARIATIONAL_BUNDLE_SCALAR_LIFT_OBSTRUCTION_NO_GLOBAL_CLOSURE",
        "strict_track_a_nonscalar_variational_bundle_transport_obstruction_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2337 strict Track-A non-scalar variational-bundle transport obstruction\n\n"
        "Status: scalar finite-part Track-A transport cannot be promoted to tensor-component transport by silent scalar lift.\n\n"
        f"- Track A coefficients: `b_R2 = {b_r2}`, `b_Ric2 = {b_ric2}`.\n"
        f"- FRW anisotropy-zero restriction is zero: `{frw_restriction_zero}`.\n"
        f"- Scalar-lift residual term counts by component: `{residual_term_counts}`.\n"
        f"- Rational anisotropic sample residual L2: `{sample_l2:.6e}`.\n"
        f"- Tracefree spatial sum remains zero: `{tracefree_spatial_sum == 0}`.\n"
        "- No full non-scalar bundle transport, no Track B GB ledger transport, no independent `a_GB`, no full/global renormalization, no QW-2191 discharge, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
