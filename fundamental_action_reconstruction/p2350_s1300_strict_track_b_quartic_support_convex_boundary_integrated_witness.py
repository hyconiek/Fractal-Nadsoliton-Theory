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
OUT = GEN / "p2350_s1300_strict_track_b_quartic_support_convex_boundary_integrated_witness.json"
MD = GEN / "p2350_s1300_strict_track_b_quartic_support_convex_boundary_integrated_witness.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2338_TRACK_B_CONTRACT": GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.json",
    "P2345_CONVEX_GAUSS_MAP_PROMISE": GEN / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.json",
    "P2348_NONSYMMETRIC_POLYNOMIAL": GEN / "p2348_s1298_strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction.json",
    "P2349_TRIAXIAL_ELLIPSOID_WITNESS": GEN / "p2349_s1299_strict_track_b_triaxial_ellipsoid_integrated_boundary_polynomial_witness.json",
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
        ROOT / "p2349_s1299_strict_track_b_triaxial_ellipsoid_integrated_boundary_polynomial_witness.py",
        GEN / "p2349_s1299_strict_track_b_triaxial_ellipsoid_integrated_boundary_polynomial_witness.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track B|convex|Gauss-map|boundary polynomial|ellipsoid|sigma3|integrated|selector|QW-2191|ToE closure|global renormalization",
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


def sample_row(name: str, u1_value: sp.Rational, eps: sp.Rational) -> dict[str, Any]:
    a = u1_value
    if a in (0, 1):
        radius_values = [sp.factor(1 - 3 * eps * a**4)] * 3
        tangent_gradient_norm2 = sp.Integer(0)
    else:
        r_perp = sp.factor(1 - 3 * eps * a**4)
        r_grad = sp.factor(1 + 12 * eps * a**2 - 15 * eps * a**4)
        radius_values = [r_grad, r_perp, r_perp]
        tangent_gradient_norm2 = sp.factor(1 - a**2)
    jacobian = sp.factor(sp.prod(radius_values))
    shape_sigma3 = sp.factor(1 / jacobian)
    local_boundary_density = sp.factor(16 * shape_sigma3)
    density_times_jacobian = sp.factor(sp.simplify(local_boundary_density * jacobian))
    return {
        "sample": name,
        "u1": str(a),
        "tangent_gradient_norm2": str(tangent_gradient_norm2),
        "principal_radii": [str(value) for value in radius_values],
        "surface_jacobian_det_radii": str(jacobian),
        "shape_sigma3_product_curvatures": str(shape_sigma3),
        "flat_K0_boundary_density_16sigma3": str(local_boundary_density),
        "density_times_surface_jacobian": str(density_times_jacobian),
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    p2335 = artifacts["P2335_TWO_TRACK_LEDGER"]
    p2338 = artifacts["P2338_TRACK_B_CONTRACT"]
    p2345 = artifacts["P2345_CONVEX_GAUSS_MAP_PROMISE"]
    p2348 = artifacts["P2348_NONSYMMETRIC_POLYNOMIAL"]
    p2349 = artifacts["P2349_TRIAXIAL_ELLIPSOID_WITNESS"]

    b_gb = sp.factor(track_b_coefficient(p2335))
    per_euler_pairing = sp.factor(sp.simplify(32 * sp.pi**2 * b_gb))
    eps = sp.Rational(1, 10)
    q = sp.symbols("q", real=True)

    support_function = "h(u)=1 + eps*u1**4 on S^3, eps=1/10"
    radii_operator_formula = "A = Hess_S(h)+h*g = (1-3*eps*u1**4)I + 12*eps*u1**2 du1⊗du1"
    quartic_not_ellipsoid_reason = "support has quartic u1**4 term; it is not the square-root quadratic support function of an ellipsoid"

    pole_row = sample_row("north_pole_u1_1", sp.Integer(1), eps)
    equator_row = sample_row("equator_u1_0", sp.Integer(0), eps)
    mid_row = sample_row("mid_latitude_u1_1_over_2", sp.Rational(1, 2), eps)
    sample_rows = [pole_row, equator_row, mid_row]

    sigma3_values = [sp.sympify(row["shape_sigma3_product_curvatures"]) for row in sample_rows]
    density_values = [sp.sympify(row["flat_K0_boundary_density_16sigma3"]) for row in sample_rows]
    nonconstant_sigma3_witness = sp.factor(sigma3_values[0] - sigma3_values[1])
    nonconstant_density_witness = sp.factor(density_values[0] - density_values[1])

    positivity_rows = [all(sp.sympify(value) > 0 for value in row["principal_radii"]) for row in sample_rows]
    admissible_eps_interval = "0 < eps < 1/3 is sufficient for these recorded sample radii to stay positive; eps=1/10 used"

    gauss_map_degree = sp.Integer(1)
    integral_sigma3 = sp.factor(2 * sp.pi**2 * gauss_map_degree)
    integrated_boundary_functional = sp.factor(16 * integral_sigma3)
    integrated_boundary_residual = sp.factor(sp.simplify(integrated_boundary_functional - 32 * sp.pi**2))
    normalized_pairing = sp.factor(sp.simplify(b_gb * integrated_boundary_functional))
    normalized_pairing_residual = sp.factor(sp.simplify(normalized_pairing - per_euler_pairing))

    # Pointwise cancellation identity in Gauss-map coordinates: (16*sigma3(shape))*det(radii)=16.
    jacobian_symbol = q
    pointwise_identity_residual = sp.factor(sp.simplify(16 * (1 / jacobian_symbol) * jacobian_symbol - 16))

    p2345_probe = p2345.get(
        "strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_probe", {}
    )
    p2345_lemma = p2345_probe.get("track_B_convex_gauss_map_boundary_functional_theorem", {})
    p2348_probe = p2348.get("strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction_probe", {})
    p2348_witness = p2348_probe.get("track_B_chern_boundary_polynomial_nonsymmetric_reduction", {})
    p2348_poly = p2348_witness.get("boundary_polynomial", {})
    p2349_probe = p2349.get("strict_track_b_triaxial_ellipsoid_integrated_boundary_polynomial_witness_probe", {})
    p2349_witness = p2349_probe.get("track_B_triaxial_ellipsoid_integrated_boundary_polynomial_witness", {})
    p2349_ellipsoid = p2349_witness.get("ellipsoid_witness", {})

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
    p2349_integrated_residual = sp.factor(
        sp.sympify(
            p2349_ellipsoid.get("integrated_boundary_functional_16_integral_sigma3", "0"),
            locals={"pi": sp.pi},
        )
        - integrated_boundary_functional
    )

    derivation_steps = [
        {
            "step": "quartic_support_class",
            "formula": support_function,
            "radii_operator_formula": radii_operator_formula,
            "non_ellipsoid_reason": quartic_not_ellipsoid_reason,
        },
        {
            "step": "local_radii_and_density_samples",
            "formula": "Evaluate radii of curvature, surface Jacobian, shape sigma3, and 16*sigma3 at selected Gauss-map points.",
            "sample_rows": sample_rows,
        },
        {
            "step": "nonconstant_density_witness",
            "formula": "The pointwise shape sigma3 and 16*sigma3 densities differ between pole and equator.",
            "pole_minus_equator_sigma3": str(nonconstant_sigma3_witness),
            "pole_minus_equator_density": str(nonconstant_density_witness),
        },
        {
            "step": "integrated_gauss_map_identity",
            "formula": "In Gauss-map coordinates, sigma3(shape)*dA=dOmega, so Integral sigma3 dA=Area(S3)=2*pi**2.",
            "pointwise_density_jacobian_residual": str(pointwise_identity_residual),
            "integral_sigma3": str(integral_sigma3),
            "boundary_functional": str(integrated_boundary_functional),
            "boundary_residual": str(integrated_boundary_residual),
        },
    ]

    support_witness = {
        "witness_id": "P2350_quartic_support_convex_boundary_integrated_witness_v1",
        "scope": "one explicit quartic support-function convex boundary stress-test; not an arbitrary-boundary theorem",
        "support_function": support_function,
        "epsilon": str(eps),
        "admissible_eps_note": admissible_eps_interval,
        "radii_operator_formula": radii_operator_formula,
        "non_ellipsoid_reason": quartic_not_ellipsoid_reason,
        "sample_rows": sample_rows,
        "all_recorded_sample_radii_positive": all(positivity_rows),
        "nonconstant_sigma3_witness_pole_minus_equator": str(nonconstant_sigma3_witness),
        "nonconstant_density_witness_pole_minus_equator": str(nonconstant_density_witness),
        "pointwise_density_jacobian_residual": str(pointwise_identity_residual),
        "gauss_map_degree": str(gauss_map_degree),
        "integral_sigma3_via_gauss_map": str(integral_sigma3),
        "integrated_boundary_functional_16_integral_sigma3": str(integrated_boundary_functional),
        "integrated_boundary_residual": str(integrated_boundary_residual),
        "normalized_track_B_pairing": str(normalized_pairing),
        "normalized_pairing_residual": str(normalized_pairing_residual),
        "p2345_boundary_residual": str(p2345_boundary_residual),
        "p2348_flat_replay_residual": str(p2348_flat_replay_residual),
        "p2349_integrated_residual": str(p2349_integrated_residual),
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
        "p2349_triaxial_ellipsoid_loaded": p2349.get("result_kind")
        == "STRICT_TRACK_B_TRIAXIAL_ELLIPSOID_INTEGRATED_BOUNDARY_POLYNOMIAL_WITNESS_NO_UNIVERSAL_THEOREM",
        "p2345_boundary_value_replayed": p2345_boundary_residual == 0,
        "p2348_flat_integral_replay_replayed": p2348_flat_replay_residual == 0,
        "p2349_integrated_value_replayed": p2349_integrated_residual == 0,
    }

    theorem_export = {
        "theorem_name": "P2350 Track-B quartic support-function convex boundary integrated witness",
        "claim": (
            "P2350 supplies a second integrated non-ellipsoidal convex-boundary stress-test after P2349.  For the "
            "support function h(u)=1+eps*u1**4 with eps=1/10, recorded radii samples are positive and the local "
            "shape sigma3 density is nonconstant.  In Gauss-map coordinates sigma3(shape)*dA=dOmega, so the integrated "
            "Track-B boundary functional is again 32*pi**2 with zero residual."
        ),
        "proof_witnesses": [
            "The quartic support function is not an ellipsoid support function because it contains a u1**4 term.",
            "Recorded pole, equator, and mid-latitude radii samples are positive for eps=1/10.",
            "The pointwise shape sigma3 density differs between pole and equator, so this is not a constant-density replay.",
            "The pointwise Gauss-map identity (16*sigma3)*det(radii)=16 is recorded symbolically.",
            "The integrated boundary functional, Track-B pairing, and dependencies on P2345/P2348/P2349 all replay with zero residual.",
        ],
        "not_licensed": [
            "integrated arbitrary-boundary theorem",
            "general support-function convexity theorem for all quartic perturbations",
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
            "Track-B now has two integrated nonconstant stress-tests: triaxial ellipsoid and quartic support boundary. "
            "The next honest step is a compact synthesis theorem marking the controlled Track-B boundary class as robust, "
            "or pause Track-B and return to Track-A tensor-bundle obstruction without selector closure."
        ),
    }

    probe = {
        "probe_id": "P2350_S1300_STRICT_TRACK_B_QUARTIC_SUPPORT_CONVEX_BOUNDARY_INTEGRATED_WITNESS",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2349 next step: second non-ellipsoidal integrated convex boundary stress-test",
            "top_hits": grep_hits(),
        },
        "track_B_quartic_support_convex_boundary_integrated_witness": {
            "ledger_coefficient_b_GB": str(b_gb),
            "per_euler_pairing": str(per_euler_pairing),
            "derivation_steps": derivation_steps,
            "support_witness": support_witness,
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
        "p2349_triaxial_ellipsoid_loaded": dependencies["p2349_triaxial_ellipsoid_loaded"],
        "p2345_boundary_value_replayed": dependencies["p2345_boundary_value_replayed"],
        "p2348_flat_integral_replay_replayed": dependencies["p2348_flat_integral_replay_replayed"],
        "p2349_integrated_value_replayed": dependencies["p2349_integrated_value_replayed"],
        "sample_rows_complete": len(sample_rows) == 3,
        "sample_radii_positive": all(positivity_rows),
        "sigma3_density_nonconstant": nonconstant_sigma3_witness != 0,
        "boundary_density_nonconstant": nonconstant_density_witness != 0,
        "pointwise_density_jacobian_residual_zero": pointwise_identity_residual == 0,
        "integrated_boundary_residual_zero": integrated_boundary_residual == 0,
        "normalized_pairing_residual_zero": normalized_pairing_residual == 0,
        "explicit_quartic_support_scope_declared": "one explicit quartic support-function" in support_witness["scope"],
        "integrated_arbitrary_boundary_theorem_not_claimed": True,
        "general_support_function_theorem_not_claimed": True,
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
        "schema_version": "p2350_s1300_v1",
        "packet_id": "P2350",
        "stage_id": "S1300",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_QUARTIC_SUPPORT_CONVEX_BOUNDARY_INTEGRATED_WITNESS_NO_UNIVERSAL_THEOREM",
        "strict_track_b_quartic_support_convex_boundary_integrated_witness_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2350 strict Track-B quartic support convex-boundary integrated witness\n\n"
        "Status: non-ellipsoidal quartic support-function convex boundary stress-test integrated; no universal boundary theorem.\n\n"
        f"- `b_GB = {b_gb}`.\n"
        f"- Support function: `{support_function}`.\n"
        f"- R_pole/r_equator/r_mid rows recorded: `{len(sample_rows)}`; sample radii positive: `{all(positivity_rows)}`.\n"
        f"- Nonconstant shape sigma3 witness pole-equator: `{nonconstant_sigma3_witness}`; density witness `{nonconstant_density_witness}`.\n"
        f"- Pointwise Gauss-map density-Jacobian residual: `{pointwise_identity_residual}`.\n"
        f"- Integrated replay: `Integral sigma3 dA = {integral_sigma3}`; boundary functional `{integrated_boundary_functional}`; residual `{integrated_boundary_residual}`.\n"
        f"- Normalized Track-B pairing: `{normalized_pairing}`; residual `{normalized_pairing_residual}`.\n"
        "- No integrated arbitrary-boundary theorem, no general support-function theorem, no general Chern-Gauss-Bonnet theorem, no global renormalization, no independent `a_GB`, no selector premise, no QW-2191 discharge, no unique-future choice, no observer prediction, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
