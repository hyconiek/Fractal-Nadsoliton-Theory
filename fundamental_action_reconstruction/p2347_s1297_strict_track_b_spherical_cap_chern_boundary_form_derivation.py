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
OUT = GEN / "p2347_s1297_strict_track_b_spherical_cap_chern_boundary_form_derivation.json"
MD = GEN / "p2347_s1297_strict_track_b_spherical_cap_chern_boundary_form_derivation.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2338_TRACK_B_CONTRACT": GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.json",
    "P2346_SPHERICAL_CAP_CANCELLATION": GEN / "p2346_s1296_strict_track_b_spherical_cap_bulk_boundary_cancellation_witness.json",
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
        ROOT / "p2346_s1296_strict_track_b_spherical_cap_bulk_boundary_cancellation_witness.py",
        GEN / "p2346_s1296_strict_track_b_spherical_cap_bulk_boundary_cancellation_witness.md",
        ROOT / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.py",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track B|spherical-cap|boundary transgression|Chern|bulk-boundary|selector|QW-2191|ToE closure|global renormalization",
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
    p2346 = artifacts["P2346_SPHERICAL_CAP_CANCELLATION"]

    c = sp.symbols("c", real=True)
    b_gb = sp.factor(track_b_coefficient(p2335))
    per_euler_pairing = sp.factor(sp.simplify(32 * sp.pi**2 * b_gb))

    s = sp.sqrt(1 - c**2)
    area_geodesic_s3 = 2 * sp.pi**2 * s**3
    geodesic_principal_curvature = c / s
    ambient_sectional_curvature = sp.Integer(1)

    chern_density_scalar = sp.factor(8 * (3 * ambient_sectional_curvature * geodesic_principal_curvature + 2 * geodesic_principal_curvature**3))
    chern_density_times_area_unsimplified = sp.factor(chern_density_scalar * area_geodesic_s3)
    chern_density_times_area = sp.factor(sp.simplify(chern_density_times_area_unsimplified))
    target_boundary_transgression = sp.factor(16 * sp.pi**2 * c * (3 - c**2))
    chern_boundary_residual = sp.factor(sp.simplify(chern_density_times_area - target_boundary_transgression))

    bulk_from_p2346_formula = sp.factor(16 * sp.pi**2 * (c - 1) ** 2 * (c + 2))
    total_with_derived_boundary = sp.factor(sp.simplify(bulk_from_p2346_formula + chern_density_times_area))
    total_residual = sp.factor(sp.simplify(total_with_derived_boundary - 32 * sp.pi**2))
    total_derivative_residual = sp.factor(sp.diff(total_with_derived_boundary, c))

    flat_limit_boundary = sp.factor(chern_density_times_area.subs(c, 1))
    hemisphere_boundary = sp.factor(chern_density_times_area.subs(c, 0))
    sample_c_half_boundary = sp.factor(chern_density_times_area.subs(c, sp.Rational(1, 2)))
    sample_c_half_bulk = sp.factor(bulk_from_p2346_formula.subs(c, sp.Rational(1, 2)))
    sample_c_half_total = sp.factor(sp.simplify(sample_c_half_bulk + sample_c_half_boundary))

    normalized_pairing = sp.factor(sp.simplify(b_gb * total_with_derived_boundary))
    normalized_pairing_residual = sp.factor(sp.simplify(normalized_pairing - per_euler_pairing))

    p2346_probe = p2346.get("strict_track_b_spherical_cap_bulk_boundary_cancellation_witness_probe", {})
    p2346_witness = p2346_probe.get("track_B_spherical_cap_bulk_boundary_cancellation_witness", {})
    p2346_cap = p2346_witness.get("spherical_cap_class", {})
    p2346_boundary_expr = sp.sympify(p2346_cap.get("boundary_transgression", "0"), locals={"pi": sp.pi, "c": c})
    p2346_total_expr = sp.sympify(p2346_cap.get("bulk_plus_boundary", "0"), locals={"pi": sp.pi, "c": c})
    p2346_boundary_residual = sp.factor(sp.simplify(chern_density_times_area - p2346_boundary_expr))
    p2346_total_residual = sp.factor(sp.simplify(total_with_derived_boundary - p2346_total_expr))

    derivation_steps = [
        {
            "step": "geodesic_boundary_geometry",
            "formula": "For a unit S^4 geodesic cap with c=cos(rho), boundary S^3 has area 2*pi**2*(1-c**2)**(3/2) and principal curvature k=c/sqrt(1-c**2).",
            "boundary_area": str(area_geodesic_s3),
            "principal_curvature": "c/sqrt(1-c**2)",
            "ambient_sectional_curvature": str(ambient_sectional_curvature),
        },
        {
            "step": "chern_boundary_scalar",
            "formula": "B_3 density scalar in this symmetric cap class is 8*(3*K_ambient*k + 2*k**3).",
            "density_scalar": str(chern_density_scalar),
        },
        {
            "step": "integrated_chern_boundary_form",
            "formula": "Integral_boundary B_3 = 8*(3*k+2*k**3)*Area(S^3_rho).",
            "density_times_area_unsimplified": str(chern_density_times_area_unsimplified),
            "integrated_boundary_transgression": str(chern_density_times_area),
            "residual_against_p2346_boundary": str(p2346_boundary_residual),
        },
        {
            "step": "bulk_boundary_cancellation_replay",
            "formula": "P2346 bulk(c) plus derived Chern boundary(c) is independent of c and equals 32*pi**2.",
            "bulk_from_p2346_formula": str(bulk_from_p2346_formula),
            "total_with_derived_boundary": str(total_with_derived_boundary),
            "total_residual": str(total_residual),
            "total_derivative_residual": str(total_derivative_residual),
        },
    ]

    chern_form_derivation = {
        "derivation_id": "P2347_spherical_cap_chern_boundary_form_derivation_v1",
        "scope": "symmetric unit-S4 geodesic-cap Chern boundary-form derivation only; not a general Chern-Gauss-Bonnet theorem",
        "parameter": "c = cos(rho), rho in (0, pi)",
        "boundary_area": str(area_geodesic_s3),
        "geodesic_principal_curvature": "c/sqrt(1-c**2)",
        "ambient_sectional_curvature": str(ambient_sectional_curvature),
        "chern_boundary_density_scalar": str(chern_density_scalar),
        "integrated_boundary_transgression": str(chern_density_times_area),
        "target_boundary_transgression": str(target_boundary_transgression),
        "chern_boundary_residual": str(chern_boundary_residual),
        "bulk_from_p2346_formula": str(bulk_from_p2346_formula),
        "total_with_derived_boundary": str(total_with_derived_boundary),
        "total_residual": str(total_residual),
        "total_derivative_residual": str(total_derivative_residual),
        "flat_limit_boundary_c1": str(flat_limit_boundary),
        "hemisphere_boundary_c0": str(hemisphere_boundary),
        "sample_c_half_bulk": str(sample_c_half_bulk),
        "sample_c_half_boundary": str(sample_c_half_boundary),
        "sample_c_half_total": str(sample_c_half_total),
        "normalized_track_B_pairing": str(normalized_pairing),
        "normalized_pairing_residual": str(normalized_pairing_residual),
        "p2346_boundary_residual": str(p2346_boundary_residual),
        "p2346_total_residual": str(p2346_total_residual),
    }

    dependencies = {
        "p2335_two_track_ledger_loaded": p2335.get("result_kind")
        == "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
        "p2338_contract_loaded": p2338.get("result_kind")
        == "STRICT_TRACK_B_GB_BOUNDARY_TOPOLOGICAL_NORMALIZATION_CONTRACT_EXPORTED_NO_NORMALIZATION_CLAIM",
        "p2346_spherical_cap_loaded": p2346.get("result_kind")
        == "STRICT_TRACK_B_SPHERICAL_CAP_CURVED_BULK_BOUNDARY_CANCELLATION_WITNESS_NO_GLOBAL_THEOREM",
        "p2346_boundary_transgression_derived": p2346_boundary_residual == 0,
        "p2346_total_replayed": p2346_total_residual == 0,
    }

    theorem_export = {
        "theorem_name": "P2347 Track-B spherical-cap Chern boundary-form derivation",
        "claim": (
            "P2347 derives the P2346 spherical-cap boundary transgression from the symmetric Chern boundary form. "
            "For a unit S^4 geodesic cap, the boundary has area 2*pi**2*(1-c**2)**(3/2), principal curvature "
            "k=c/sqrt(1-c**2), and the cap-class Chern scalar 8*(3*k+2*k**3).  Multiplying and simplifying gives "
            "16*pi**2*c*(3-c**2), exactly the P2346 boundary term, and replaying the P2346 bulk term gives total 32*pi**2."
        ),
        "proof_witnesses": [
            "The geodesic S^3 boundary area and equal-principal-curvature data are made explicit.",
            "The symmetric Chern boundary scalar 8*(3*k+2*k**3) is multiplied by the area element.",
            "The square-root factors cancel algebraically to produce 16*pi**2*c*(3-c**2).",
            "The derived boundary term exactly matches P2346's stored boundary transgression.",
            "The derived boundary term plus P2346 bulk term has zero c-derivative and residual zero against 32*pi**2.",
        ],
        "not_licensed": [
            "general Chern boundary form for arbitrary hypersurfaces",
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
            "The immediate P2346 debt is now discharged for the spherical-cap symmetry class.  The next honest "
            "proof/computational step is either derive the full non-symmetric Chern boundary polynomial and compare it "
            "to this cap reduction, or stop Track-B at controlled classes and return to the Track-A tensor-bundle obstruction."
        ),
    }

    probe = {
        "probe_id": "P2347_S1297_STRICT_TRACK_B_SPHERICAL_CAP_CHERN_BOUNDARY_FORM_DERIVATION",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2346 next step: derive spherical-cap boundary transgression from symmetric Chern boundary form",
            "top_hits": grep_hits(),
        },
        "track_B_spherical_cap_chern_boundary_form_derivation": {
            "ledger_coefficient_b_GB": str(b_gb),
            "per_euler_pairing": str(per_euler_pairing),
            "derivation_steps": derivation_steps,
            "chern_form_derivation": chern_form_derivation,
        },
        "current_export_dependencies": dependencies,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "repo_grep_audit_recorded": len(probe["repo_grep_audit"]["top_hits"]) > 0,
        "p2335_two_track_ledger_loaded": dependencies["p2335_two_track_ledger_loaded"],
        "p2338_contract_loaded": dependencies["p2338_contract_loaded"],
        "p2346_spherical_cap_loaded": dependencies["p2346_spherical_cap_loaded"],
        "p2346_boundary_transgression_derived": dependencies["p2346_boundary_transgression_derived"],
        "p2346_total_replayed": dependencies["p2346_total_replayed"],
        "chern_boundary_residual_zero": chern_boundary_residual == 0,
        "bulk_boundary_total_is_32pi2": total_residual == 0,
        "total_derivative_zero": total_derivative_residual == 0,
        "flat_limit_boundary_is_32pi2": flat_limit_boundary == 32 * sp.pi**2,
        "hemisphere_boundary_zero": hemisphere_boundary == 0,
        "sample_c_half_total_correct": sample_c_half_total == 32 * sp.pi**2,
        "normalized_pairing_residual_zero": normalized_pairing_residual == 0,
        "symmetric_cap_scope_declared": "symmetric unit-S4 geodesic-cap" in chern_form_derivation["scope"],
        "general_chern_boundary_form_not_claimed": True,
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
        "schema_version": "p2347_s1297_v1",
        "packet_id": "P2347",
        "stage_id": "S1297",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_SPHERICAL_CAP_CHERN_BOUNDARY_FORM_DERIVATION_NO_GENERAL_CGB_THEOREM",
        "strict_track_b_spherical_cap_chern_boundary_form_derivation_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2347 strict Track-B spherical-cap Chern boundary-form derivation\n\n"
        "Status: P2346 spherical-cap boundary transgression derived from the symmetric Chern boundary form; no general CGB theorem.\n\n"
        f"- `b_GB = {b_gb}`.\n"
        f"- Boundary geometry: area `{area_geodesic_s3}`, principal curvature `c/sqrt(1-c**2)`.\n"
        f"- Chern scalar: `8*(3*k+2*k**3)` gives integrated boundary term `{chern_density_times_area}`.\n"
        f"- Residual against P2346 boundary transgression: `{p2346_boundary_residual}`.\n"
        f"- P2346 bulk plus derived boundary: `{total_with_derived_boundary}`; residual `{total_residual}`; derivative residual `{total_derivative_residual}`.\n"
        f"- Flat limit boundary: `{flat_limit_boundary}`; hemisphere boundary: `{hemisphere_boundary}`; normalized pairing `{normalized_pairing}`.\n"
        "- No general Chern boundary form, no general Chern-Gauss-Bonnet theorem, no global renormalization, no independent `a_GB`, no selector premise, no QW-2191 discharge, no unique-future choice, no observer prediction, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
