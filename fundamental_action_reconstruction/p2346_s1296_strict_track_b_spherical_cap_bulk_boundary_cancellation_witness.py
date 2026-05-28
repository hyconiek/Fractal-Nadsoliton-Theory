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
OUT = GEN / "p2346_s1296_strict_track_b_spherical_cap_bulk_boundary_cancellation_witness.json"
MD = GEN / "p2346_s1296_strict_track_b_spherical_cap_bulk_boundary_cancellation_witness.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2338_TRACK_B_CONTRACT": GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.json",
    "P2341_D4_BOUNDARY_FIXTURE": GEN / "p2341_s1291_strict_track_b_d4_boundary_correction_fixture_witness.json",
    "P2345_CONVEX_GAUSS_MAP_PROMISE": GEN / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.json",
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
        ROOT / "p2341_s1291_strict_track_b_d4_boundary_correction_fixture_witness.py",
        ROOT / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.py",
        GEN / "p2341_s1291_strict_track_b_d4_boundary_correction_fixture_witness.md",
        GEN / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track B|boundary functional|convex|Gauss-map|bulk|boundary|global renormalization|selector|QW-2191|ToE closure",
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
    p2341 = artifacts["P2341_D4_BOUNDARY_FIXTURE"]
    p2345 = artifacts["P2345_CONVEX_GAUSS_MAP_PROMISE"]

    c = sp.symbols("c", real=True)
    b_gb = sp.factor(track_b_coefficient(p2335))
    per_euler_pairing = sp.factor(sp.simplify(32 * sp.pi**2 * b_gb))

    chi_cap = sp.Integer(1)
    unit_s4_gb_density = sp.Integer(24)
    cap_volume = sp.factor(2 * sp.pi**2 * (sp.Rational(2, 3) - c + c**3 / 3))
    bulk_gb_integral = sp.factor(sp.simplify(unit_s4_gb_density * cap_volume))
    boundary_transgression = sp.factor(16 * sp.pi**2 * c * (3 - c**2))
    bulk_plus_boundary = sp.factor(sp.simplify(bulk_gb_integral + boundary_transgression))
    target_topological_number = sp.factor(32 * sp.pi**2 * chi_cap)
    cancellation_residual = sp.factor(sp.simplify(bulk_plus_boundary - target_topological_number))
    c_derivative_residual = sp.factor(sp.diff(bulk_plus_boundary, c))

    boundary_area_factor = sp.factor(2 * sp.pi**2 * (1 - c**2) ** sp.Rational(3, 2))
    principal_curvature_parameter = "k = c/sqrt(1-c**2)"
    curvature_adjusted_boundary_density = "8*(3*k + 2*k**3) dA on geodesic S^3 boundary"

    flat_limit_bulk = sp.factor(bulk_gb_integral.subs(c, 1))
    flat_limit_boundary = sp.factor(boundary_transgression.subs(c, 1))
    hemisphere_bulk = sp.factor(bulk_gb_integral.subs(c, 0))
    hemisphere_boundary = sp.factor(boundary_transgression.subs(c, 0))
    complement_sample_bulk = sp.factor(bulk_gb_integral.subs(c, -sp.Rational(1, 2)))
    complement_sample_boundary = sp.factor(boundary_transgression.subs(c, -sp.Rational(1, 2)))
    complement_sample_total = sp.factor(sp.simplify(complement_sample_bulk + complement_sample_boundary))

    normalized_cap_pairing = sp.factor(sp.simplify(b_gb * bulk_plus_boundary))
    normalized_pairing_residual = sp.factor(sp.simplify(normalized_cap_pairing - per_euler_pairing * chi_cap))

    p2341_probe = p2341.get("strict_track_b_d4_boundary_correction_fixture_witness_probe", {})
    p2341_witness = p2341_probe.get("track_B_D4_boundary_fixture_witness", {})
    p2345_probe = p2345.get(
        "strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_probe", {}
    )
    p2345_lemma = p2345_probe.get("track_B_convex_gauss_map_boundary_functional_theorem", {})

    derivation_steps = [
        {
            "step": "curved_bulk_fixture",
            "formula": "Unit S^4 has GB density 24 and geodesic cap volume 2*pi**2*(2/3 - c + c**3/3), c=cos(rho).",
            "gb_density": str(unit_s4_gb_density),
            "cap_volume": str(cap_volume),
            "bulk_gb_integral": str(bulk_gb_integral),
        },
        {
            "step": "boundary_transgression_fixture",
            "formula": "For this one-parameter spherical-cap class, boundary transgression is 16*pi**2*c*(3-c**2).",
            "boundary_area_factor": str(boundary_area_factor),
            "principal_curvature_parameter": principal_curvature_parameter,
            "curvature_adjusted_boundary_density": curvature_adjusted_boundary_density,
            "boundary_transgression": str(boundary_transgression),
        },
        {
            "step": "radius_cancellation",
            "formula": "bulk(c) + boundary(c) = 32*pi**2 and d/dc of the total is zero.",
            "bulk_plus_boundary": str(bulk_plus_boundary),
            "target_topological_number": str(target_topological_number),
            "cancellation_residual": str(cancellation_residual),
            "c_derivative_residual": str(c_derivative_residual),
        },
        {
            "step": "flat_limit_and_hemisphere_checks",
            "formula": "c=1 rederives the flat D4 boundary value; c=0 moves all charge into the bulk hemisphere.",
            "flat_limit_bulk": str(flat_limit_bulk),
            "flat_limit_boundary": str(flat_limit_boundary),
            "hemisphere_bulk": str(hemisphere_bulk),
            "hemisphere_boundary": str(hemisphere_boundary),
        },
    ]

    spherical_cap_class = {
        "class_id": "unit_S4_geodesic_four_ball_cap_c_in_minus1_1",
        "scope": "one-parameter curved bulk-plus-boundary spherical-cap class only; not arbitrary curved manifolds",
        "radius_parameter": "rho in (0, pi), c = cos(rho)",
        "chi_cap": str(chi_cap),
        "unit_s4_gb_density": str(unit_s4_gb_density),
        "cap_volume": str(cap_volume),
        "bulk_gb_integral": str(bulk_gb_integral),
        "boundary_transgression": str(boundary_transgression),
        "bulk_plus_boundary": str(bulk_plus_boundary),
        "target_topological_number_32pi2_chi": str(target_topological_number),
        "cancellation_residual": str(cancellation_residual),
        "c_derivative_residual": str(c_derivative_residual),
        "flat_limit_bulk_c1": str(flat_limit_bulk),
        "flat_limit_boundary_c1": str(flat_limit_boundary),
        "hemisphere_bulk_c0": str(hemisphere_bulk),
        "hemisphere_boundary_c0": str(hemisphere_boundary),
        "sample_c_minus_half_bulk": str(complement_sample_bulk),
        "sample_c_minus_half_boundary": str(complement_sample_boundary),
        "sample_c_minus_half_total": str(complement_sample_total),
        "normalized_track_B_pairing": str(normalized_cap_pairing),
        "normalized_pairing_residual": str(normalized_pairing_residual),
    }

    dependencies = {
        "p2335_two_track_ledger_loaded": p2335.get("result_kind")
        == "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
        "p2338_contract_loaded": p2338.get("result_kind")
        == "STRICT_TRACK_B_GB_BOUNDARY_TOPOLOGICAL_NORMALIZATION_CONTRACT_EXPORTED_NO_NORMALIZATION_CLAIM",
        "p2341_d4_boundary_fixture_loaded": p2341.get("result_kind")
        == "STRICT_TRACK_B_D4_NONZERO_BOUNDARY_FIXTURE_WITNESS_NO_UNIVERSAL_BOUNDARY_THEOREM",
        "p2345_convex_gauss_map_loaded": p2345.get("result_kind")
        == "STRICT_TRACK_B_CONVEX_GAUSS_MAP_BOUNDARY_FUNCTIONAL_THEOREM_AND_CLOSURE_PROMISE_MAP_NO_SELECTOR_CLOSURE",
        "p2341_flat_limit_boundary_rederived": p2341_witness.get("boundary_correction_fixture_value")
        == str(flat_limit_boundary),
        "p2345_pairing_rederived": p2345_lemma.get("normalized_track_B_pairing") == str(normalized_cap_pairing),
    }

    theorem_export = {
        "theorem_name": "P2346 Track-B spherical-cap curved bulk-boundary cancellation witness",
        "claim": (
            "P2346 is the first Track-B curved bulk-plus-boundary computation after the flat convex Gauss-map lane. "
            "For a unit S^4 geodesic four-ball cap with c=cos(rho), the bulk GB integral is "
            "24*2*pi**2*(2/3-c+c**3/3), the declared spherical-cap boundary transgression is "
            "16*pi**2*c*(3-c**2), and their sum is 32*pi**2 for every c in the one-parameter cap class."
        ),
        "proof_witnesses": [
            "The bulk term is computed from the constant-curvature unit S^4 GB density 24 and the cap volume integral.",
            "The boundary transgression is recorded as the curvature-adjusted spherical-cap boundary density integrated over geodesic S^3.",
            "The polynomial c-dependence cancels exactly: 32*pi**2 - 48*pi**2*c + 16*pi**2*c**3 + 48*pi**2*c - 16*pi**2*c**3 = 32*pi**2.",
            "The flat c=1 endpoint rederives the P2341 D4 boundary correction while the c=0 hemisphere has zero boundary term and bulk 32*pi**2.",
            "The Track-B ledger pairing residual against per-Euler normalization is zero.",
        ],
        "not_licensed": [
            "general curved bulk-plus-boundary theorem for arbitrary four-manifolds",
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
            "This is the strongest current Track-B closure promise because it adds an exact curved bulk-boundary "
            "cancellation to the earlier flat boundary computations.  Next either derive the spherical-cap boundary "
            "transgression from the full Chern boundary form, or return to Track-A tensor-bundle obstruction; do not move "
            "to selector closure until these structural proof obligations are exhausted."
        ),
    }

    probe = {
        "probe_id": "P2346_S1296_STRICT_TRACK_B_SPHERICAL_CAP_BULK_BOUNDARY_CANCELLATION_WITNESS",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2345 next step: controlled curved spherical-cap bulk-plus-boundary cancellation",
            "top_hits": grep_hits(),
        },
        "track_B_spherical_cap_bulk_boundary_cancellation_witness": {
            "ledger_coefficient_b_GB": str(b_gb),
            "per_euler_pairing": str(per_euler_pairing),
            "derivation_steps": derivation_steps,
            "spherical_cap_class": spherical_cap_class,
        },
        "current_export_dependencies": dependencies,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "repo_grep_audit_recorded": len(probe["repo_grep_audit"]["top_hits"]) > 0,
        "p2335_two_track_ledger_loaded": dependencies["p2335_two_track_ledger_loaded"],
        "p2338_contract_loaded": dependencies["p2338_contract_loaded"],
        "p2341_d4_boundary_fixture_loaded": dependencies["p2341_d4_boundary_fixture_loaded"],
        "p2345_convex_gauss_map_loaded": dependencies["p2345_convex_gauss_map_loaded"],
        "p2341_flat_limit_boundary_rederived": dependencies["p2341_flat_limit_boundary_rederived"],
        "p2345_pairing_rederived": dependencies["p2345_pairing_rederived"],
        "bulk_density_24_declared": unit_s4_gb_density == 24,
        "bulk_boundary_total_is_32pi2": cancellation_residual == 0,
        "radius_derivative_zero": c_derivative_residual == 0,
        "flat_limit_rederives_d4_boundary": flat_limit_bulk == 0 and flat_limit_boundary == 32 * sp.pi**2,
        "hemisphere_moves_charge_to_bulk": hemisphere_bulk == 32 * sp.pi**2 and hemisphere_boundary == 0,
        "sample_complement_cap_total_correct": complement_sample_total == 32 * sp.pi**2,
        "normalized_pairing_residual_zero": normalized_pairing_residual == 0,
        "one_parameter_curved_scope_declared": "one-parameter curved" in spherical_cap_class["scope"],
        "general_curved_theorem_not_claimed": True,
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
        "schema_version": "p2346_s1296_v1",
        "packet_id": "P2346",
        "stage_id": "S1296",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_SPHERICAL_CAP_CURVED_BULK_BOUNDARY_CANCELLATION_WITNESS_NO_GLOBAL_THEOREM",
        "strict_track_b_spherical_cap_bulk_boundary_cancellation_witness_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2346 strict Track-B spherical-cap curved bulk-boundary cancellation witness\n\n"
        "Status: one curved S4 geodesic-cap bulk-plus-boundary cancellation computed; no global curved theorem.\n\n"
        f"- `b_GB = {b_gb}`.\n"
        f"- Cap parameter: `c=cos(rho)`; bulk GB integral `{bulk_gb_integral}`.\n"
        f"- Boundary transgression: `{boundary_transgression}`.\n"
        f"- Bulk plus boundary: `{bulk_plus_boundary}`; residual against `32*pi**2*chi`: `{cancellation_residual}`; derivative residual `{c_derivative_residual}`.\n"
        f"- Flat limit `c=1`: bulk `{flat_limit_bulk}`, boundary `{flat_limit_boundary}`; hemisphere `c=0`: bulk `{hemisphere_bulk}`, boundary `{hemisphere_boundary}`.\n"
        f"- Normalized Track-B pairing: `{normalized_cap_pairing}`; residual `{normalized_pairing_residual}`.\n"
        "- No general curved bulk-boundary theorem, no general Chern-Gauss-Bonnet theorem, no global renormalization, no independent `a_GB`, no selector premise, no QW-2191 discharge, no unique-future choice, no observer prediction, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
