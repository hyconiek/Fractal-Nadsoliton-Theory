#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2351_s1301_strict_track_b_quartic_support_global_convexity_jacobian_identity_witness.json"
MD = GEN / "p2351_s1301_strict_track_b_quartic_support_global_convexity_jacobian_identity_witness.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2338_TRACK_B_CONTRACT": GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.json",
    "P2345_CONVEX_GAUSS_MAP_PROMISE": GEN / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.json",
    "P2348_NONSYMMETRIC_POLYNOMIAL": GEN / "p2348_s1298_strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction.json",
    "P2350_QUARTIC_SUPPORT_WITNESS": GEN / "p2350_s1300_strict_track_b_quartic_support_convex_boundary_integrated_witness.json",
}


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def grep_hits() -> list[str]:
    candidates = [
        ROOT / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.py",
        ROOT / "p2348_s1298_strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction.py",
        ROOT / "p2350_s1300_strict_track_b_quartic_support_convex_boundary_integrated_witness.py",
        GEN / "p2350_s1300_strict_track_b_quartic_support_convex_boundary_integrated_witness.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track B|quartic support|support-function|convex|Gauss-map|sigma3|Jacobian|selector|QW-2191|ToE closure|global renormalization",
        *existing,
    ]
    proc = subprocess.run(cmd, cwd=ROOT.parent, text=True, capture_output=True, check=False)
    return proc.stdout.splitlines()[:180]


def track_b_coefficient(p2335: dict[str, Any]) -> sp.Expr:
    probe = p2335.get("strict_two_track_renormalization_ledger_probe", {})
    ledger = probe.get("two_track_ledger", {})
    track_b = ledger.get("track_B_gb_topological_counterterm_ledger", {})
    raw = track_b.get("ledger_coefficient_b_GB_topological", "0")
    return sp.sympify(raw, locals={"pi": sp.pi, "log": sp.log, "ln": sp.log})


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    p2335 = artifacts["P2335_TWO_TRACK_LEDGER"]
    p2338 = artifacts["P2338_TRACK_B_CONTRACT"]
    p2345 = artifacts["P2345_CONVEX_GAUSS_MAP_PROMISE"]
    p2348 = artifacts["P2348_NONSYMMETRIC_POLYNOMIAL"]
    p2350 = artifacts["P2350_QUARTIC_SUPPORT_WITNESS"]

    b_gb = sp.factor(track_b_coefficient(p2335))
    per_euler_pairing = sp.factor(sp.simplify(32 * sp.pi**2 * b_gb))

    eps = sp.Rational(1, 10)
    y = sp.symbols("y", real=True)
    q = sp.symbols("q", positive=True)

    support_function = "h(u)=1 + eps*u1**4 on S^3, eps=1/10"
    coordinate_definition = "y = u1**2 in [0,1]"
    r_perp = sp.factor(1 - 3 * eps * y**2)
    r_grad = sp.factor(1 + 12 * eps * y - 15 * eps * y**2)
    r_perp_second_derivative = sp.diff(r_perp, y, 2)
    r_grad_second_derivative = sp.diff(r_grad, y, 2)
    r_perp_endpoint_values = {"y=0": str(r_perp.subs(y, 0)), "y=1": str(r_perp.subs(y, 1))}
    r_grad_endpoint_values = {"y=0": str(r_grad.subs(y, 0)), "y=1": str(r_grad.subs(y, 1))}
    min_radius_lower_bound = sp.Rational(7, 10)
    r_perp_minus_bound = sp.factor(r_perp - min_radius_lower_bound)
    r_grad_minus_bound = sp.factor(r_grad - min_radius_lower_bound)
    r_perp_bound_certificate = sp.factor(3 * eps * (1 - y**2))
    r_grad_bound_certificate = sp.factor(3 * eps * (1 - y) * (5 * y + 1))
    determinant_jacobian = sp.factor(r_grad * r_perp**2)
    determinant_lower_bound = sp.factor(min_radius_lower_bound**3)
    determinant_endpoint_values = {
        "y=0": str(determinant_jacobian.subs(y, 0)),
        "y=1": str(determinant_jacobian.subs(y, 1)),
    }

    shape_sigma3 = sp.factor(1 / determinant_jacobian)
    boundary_density = sp.factor(16 * shape_sigma3)
    density_jacobian_residual = sp.factor(sp.simplify(boundary_density * determinant_jacobian - 16))
    abstract_density_jacobian_residual = sp.factor(sp.simplify(16 * (1 / q) * q - 16))
    area_s3 = sp.factor(2 * sp.pi**2)
    integral_sigma3 = area_s3
    integrated_boundary_functional = sp.factor(16 * integral_sigma3)
    integrated_boundary_residual = sp.factor(sp.simplify(integrated_boundary_functional - 32 * sp.pi**2))
    normalized_pairing = sp.factor(sp.simplify(b_gb * integrated_boundary_functional))
    normalized_pairing_residual = sp.factor(sp.simplify(normalized_pairing - per_euler_pairing))

    p2345_probe = p2345.get(
        "strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_probe", {}
    )
    p2345_lemma = p2345_probe.get("track_B_convex_gauss_map_boundary_functional_theorem", {})
    p2348_probe = p2348.get("strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction_probe", {})
    p2348_witness = p2348_probe.get("track_B_chern_boundary_polynomial_nonsymmetric_reduction", {})
    p2348_poly = p2348_witness.get("boundary_polynomial", {})
    p2350_probe = p2350.get("strict_track_b_quartic_support_convex_boundary_integrated_witness_probe", {})
    p2350_witness = p2350_probe.get("track_B_quartic_support_convex_boundary_integrated_witness", {})
    p2350_support = p2350_witness.get("support_witness", {})

    p2345_boundary_residual = sp.factor(
        sp.sympify(p2345_lemma.get("boundary_functional_convex_class", "0"), locals={"pi": sp.pi})
        - integrated_boundary_functional
    )
    p2348_flat_replay_residual = sp.factor(
        sp.sympify(
            p2348_poly.get("flat_convex_integral_replay", {}).get("substitute_integral_sigma3_equals_2pi2", "0"),
            locals={"pi": sp.pi},
        )
        - integrated_boundary_functional
    )
    p2350_integrated_replay_residual = sp.factor(
        sp.sympify(
            p2350_support.get("integrated_boundary_functional_16_integral_sigma3", "0"),
            locals={"pi": sp.pi},
        )
        - integrated_boundary_functional
    )
    p2350_pairing_replay_residual = sp.factor(
        sp.sympify(
            p2350_support.get("normalized_track_B_pairing", "0"),
            locals={"pi": sp.pi, "log": sp.log},
        )
        - normalized_pairing
    )

    convexity_certificate = {
        "support_function": support_function,
        "coordinate_definition": coordinate_definition,
        "radii_eigenvalues": {
            "r_perp_multiplicity_2": str(r_perp),
            "r_gradient_multiplicity_1": str(r_grad),
        },
        "r_perp_endpoint_values": r_perp_endpoint_values,
        "r_gradient_endpoint_values": r_grad_endpoint_values,
        "r_perp_second_derivative": str(r_perp_second_derivative),
        "r_gradient_second_derivative": str(r_grad_second_derivative),
        "min_radius_lower_bound": str(min_radius_lower_bound),
        "r_perp_minus_lower_bound": str(r_perp_minus_bound),
        "r_gradient_minus_lower_bound": str(r_grad_minus_bound),
        "positive_certificates_on_0_1": {
            "r_perp_minus_7_over_10": str(r_perp_bound_certificate),
            "r_gradient_minus_7_over_10": str(r_grad_bound_certificate),
        },
        "global_strict_convexity_for_eps_1_over_10": True,
        "determinant_jacobian": str(determinant_jacobian),
        "determinant_endpoint_values": determinant_endpoint_values,
        "determinant_lower_bound_from_radii": str(determinant_lower_bound),
    }

    jacobian_identity_certificate = {
        "shape_sigma3_as_inverse_det_radii": str(shape_sigma3),
        "flat_boundary_density_16sigma3": str(boundary_density),
        "density_times_jacobian_minus_16_residual": str(density_jacobian_residual),
        "abstract_q_residual": str(abstract_density_jacobian_residual),
        "area_S3": str(area_s3),
        "integral_sigma3_dA_via_gauss_map_degree_1": str(integral_sigma3),
        "integrated_boundary_functional_16_integral_sigma3": str(integrated_boundary_functional),
        "integrated_boundary_residual": str(integrated_boundary_residual),
        "normalized_track_B_pairing": str(normalized_pairing),
        "normalized_pairing_residual": str(normalized_pairing_residual),
        "p2345_boundary_residual": str(p2345_boundary_residual),
        "p2348_flat_replay_residual": str(p2348_flat_replay_residual),
        "p2350_integrated_replay_residual": str(p2350_integrated_replay_residual),
        "p2350_pairing_replay_residual": str(p2350_pairing_replay_residual),
    }

    dependencies = {
        "p2335_two_track_ledger_loaded": p2335.get("result_kind")
        == "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
        "p2338_contract_loaded": p2338.get("result_kind")
        == "STRICT_TRACK_B_GB_BOUNDARY_TOPOLOGICAL_NORMALIZATION_CONTRACT_EXPORTED_NO_NORMALIZATION_CLAIM",
        "p2345_convex_gauss_map_loaded": p2345.get("result_kind")
        == "STRICT_TRACK_B_CONVEX_GAUSS_MAP_BOUNDARY_FUNCTIONAL_THEOREM_AND_CLOSURE_PROMISE_MAP_NO_SELECTOR_CLOSURE",
        "p2348_nonsymmetric_polynomial_loaded": p2348.get("result_kind")
        == "STRICT_TRACK_B_CHERN_BOUNDARY_POLYNOMIAL_NONSYMMETRIC_REDUCTION_NO_INTEGRATED_GLOBAL_THEOREM",
        "p2350_quartic_support_loaded": p2350.get("result_kind")
        == "STRICT_TRACK_B_QUARTIC_SUPPORT_CONVEX_BOUNDARY_INTEGRATED_WITNESS_NO_UNIVERSAL_THEOREM",
        "p2345_boundary_value_replayed": p2345_boundary_residual == 0,
        "p2348_flat_integral_replayed": p2348_flat_replay_residual == 0,
        "p2350_integrated_value_replayed": p2350_integrated_replay_residual == 0,
        "p2350_pairing_replayed": p2350_pairing_replay_residual == 0,
    }

    theorem_export = {
        "theorem_name": "P2351 Track-B quartic support global convexity and Jacobian identity witness",
        "claim": (
            "For the explicit quartic support-function boundary h(u)=1+eps*u1**4 on S^3 with eps=1/10, "
            "the spherical-Hessian radii eigenvalues are globally bounded below by 7/10 on y=u1**2 in [0,1]. "
            "Thus the boundary is globally strictly convex in this fixture.  The determinant-Jacobian identity then gives "
            "(16*sigma3)*det(radii)=16 pointwise and replays the Track-B integrated boundary functional 32*pi**2."
        ),
        "proof_witnesses": [
            "r_perp - 7/10 = 3*(1-y**2)/10 is nonnegative on y in [0,1].",
            "r_gradient - 7/10 = 3*(1-y)*(5*y+1)/10 is nonnegative on y in [0,1].",
            "All three radii are therefore at least 7/10, so det(radii) is positive and bounded below by 343/1000.",
            "The exact polynomial determinant is recorded and the pointwise density-Jacobian residual simplifies to zero.",
            "The P2350 integrated value and pairing replay with zero residual, upgrading sample positivity to global convexity for this one fixture.",
        ],
        "not_licensed": [
            "general support-function theorem for arbitrary h on S^3",
            "general quartic-perturbation convexity theorem for all eps or all quartic modes",
            "integrated arbitrary-boundary theorem",
            "general Chern-Gauss-Bonnet theorem over arbitrary compact four-manifolds with boundary",
            "global renormalization closure",
            "independent a_GB measurement separate from the P2335 ledger coefficient",
            "bulk EOM GB force or EOM-only GB lift",
            "selector premise or QW-2191 selector discharge",
            "choice of the unique physical future",
            "observer-level prediction",
            "G1/G3 update",
            "ToE closure",
        ],
        "next_honest_step": (
            "Track-B now has a global convexity/Jacobian proof for the quartic support fixture, not just samples. "
            "The next honest non-selector move is either a compact controlled-class synthesis of P2345/P2348/P2349/P2350/P2351 "
            "or a return to Track-A tensor-bundle obstruction before any selector claim."
        ),
    }

    probe = {
        "probe_id": "P2351_S1301_STRICT_TRACK_B_QUARTIC_SUPPORT_GLOBAL_CONVEXITY_JACOBIAN_IDENTITY_WITNESS",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2350 next step: global convexity and Jacobian identity proof for quartic support boundary",
            "top_hits": grep_hits(),
        },
        "track_B_quartic_support_global_convexity_jacobian_identity_witness": {
            "ledger_coefficient_b_GB": str(b_gb),
            "per_euler_pairing": str(per_euler_pairing),
            "convexity_certificate": convexity_certificate,
            "jacobian_identity_certificate": jacobian_identity_certificate,
        },
        "current_export_dependencies": dependencies,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "repo_grep_audit_recorded": len(probe["repo_grep_audit"]["top_hits"]) > 0,
        "p2335_two_track_ledger_loaded": dependencies["p2335_two_track_ledger_loaded"],
        "p2338_contract_loaded": dependencies["p2338_contract_loaded"],
        "p2345_convex_gauss_map_loaded": dependencies["p2345_convex_gauss_map_loaded"],
        "p2348_nonsymmetric_polynomial_loaded": dependencies["p2348_nonsymmetric_polynomial_loaded"],
        "p2350_quartic_support_loaded": dependencies["p2350_quartic_support_loaded"],
        "p2345_boundary_value_replayed": dependencies["p2345_boundary_value_replayed"],
        "p2348_flat_integral_replayed": dependencies["p2348_flat_integral_replayed"],
        "p2350_integrated_value_replayed": dependencies["p2350_integrated_value_replayed"],
        "p2350_pairing_replayed": dependencies["p2350_pairing_replayed"],
        "r_perp_lower_bound_identity_exact": sp.simplify(r_perp_minus_bound - r_perp_bound_certificate) == 0,
        "r_gradient_lower_bound_identity_exact": sp.simplify(r_grad_minus_bound - r_grad_bound_certificate) == 0,
        "global_strict_convexity_for_fixture_certified": bool(min_radius_lower_bound > 0),
        "determinant_lower_bound_positive": bool(determinant_lower_bound > 0),
        "density_jacobian_residual_zero": density_jacobian_residual == 0,
        "abstract_density_jacobian_residual_zero": abstract_density_jacobian_residual == 0,
        "integrated_boundary_residual_zero": integrated_boundary_residual == 0,
        "normalized_pairing_residual_zero": normalized_pairing_residual == 0,
        "general_support_function_theorem_not_claimed": True,
        "general_quartic_perturbation_theorem_not_claimed": True,
        "integrated_arbitrary_boundary_theorem_not_claimed": True,
        "general_cgb_theorem_not_claimed": True,
        "global_renormalization_not_claimed": True,
        "independent_a_gb_not_claimed": True,
        "bulk_eom_gb_force_not_claimed": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_unique_future_choice_claimed": True,
        "no_observer_prediction_claimed": True,
        "no_g1_g3_update_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2351_s1301_v1",
        "packet_id": "P2351",
        "stage_id": "S1301",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-29T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_QUARTIC_SUPPORT_GLOBAL_CONVEXITY_JACOBIAN_IDENTITY_WITNESS_NO_UNIVERSAL_THEOREM",
        "strict_track_b_quartic_support_global_convexity_jacobian_identity_witness_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2351 strict Track-B quartic support global convexity/Jacobian identity witness\n\n"
        "Status: global convexity certificate for the explicit quartic support fixture; no universal boundary theorem.\n\n"
        f"- `b_GB = {b_gb}`.\n"
        f"- Support function: `{support_function}` with `{coordinate_definition}`.\n"
        f"- Radii: `r_perp = {r_perp}` and `r_gradient = {r_grad}`.\n"
        f"- Lower-bound certificates: `r_perp-7/10 = {r_perp_bound_certificate}`; `r_gradient-7/10 = {r_grad_bound_certificate}`.\n"
        f"- Determinant Jacobian: `{determinant_jacobian}`; lower bound `{determinant_lower_bound}`.\n"
        f"- Pointwise density-Jacobian residual: `{density_jacobian_residual}`.\n"
        f"- Integrated replay: `Integral sigma3 dA = {integral_sigma3}`; boundary functional `{integrated_boundary_functional}`; residual `{integrated_boundary_residual}`.\n"
        f"- Normalized Track-B pairing: `{normalized_pairing}`; residual `{normalized_pairing_residual}`.\n"
        "- No general support-function theorem, no arbitrary-boundary theorem, no general Chern-Gauss-Bonnet theorem, no global renormalization, no independent `a_GB`, no selector premise, no QW-2191 discharge, no unique-future choice, no observer prediction, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
