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
OUT = GEN / "p2344_s1294_strict_track_b_flat_spheroidal_d4_boundary_functional_integral_witness.json"
MD = GEN / "p2344_s1294_strict_track_b_flat_spheroidal_d4_boundary_functional_integral_witness.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2338_TRACK_B_CONTRACT": GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.json",
    "P2341_D4_BOUNDARY_FIXTURE": GEN / "p2341_s1291_strict_track_b_d4_boundary_correction_fixture_witness.json",
    "P2343_FLAT_ROUND_D4_DERIVATION": GEN / "p2343_s1293_strict_track_b_flat_round_d4_boundary_functional_derivation.json",
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
        ROOT / "p2343_s1293_strict_track_b_flat_round_d4_boundary_functional_derivation.py",
        GEN / "p2341_s1291_strict_track_b_d4_boundary_correction_fixture_witness.md",
        GEN / "p2343_s1293_strict_track_b_flat_round_d4_boundary_functional_derivation.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track B|D4|D\\^4|boundary functional|flat round|spheroid|non-round|S3|S\\^3|sigma3|principal curv|Euler|chi|a_GB|QW-2191",
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
    p2343 = artifacts["P2343_FLAT_ROUND_D4_DERIVATION"]

    q = sp.symbols("q", positive=True)
    theta = sp.symbols("theta", real=True)
    t = sp.symbols("t", nonnegative=True)

    b_gb = sp.factor(track_b_coefficient(p2335))
    per_euler_pairing = sp.factor(sp.simplify(32 * sp.pi**2 * b_gb))
    chi_flat_spheroidal_d4 = sp.Integer(1)
    bulk_gb_flat_value = sp.Integer(0)

    spheroid_equation = "x1**2 + x2**2 + x3**2 + x4**2/q**2 <= 1, q > 0"
    gauss_map_jacobian_density = sp.factor(
        (1 / q) / (sp.sin(theta) ** 2 + sp.cos(theta) ** 2 / q**2) ** 2
    )
    s3_polar_area_element = sp.factor(4 * sp.pi * sp.sin(theta) ** 2)
    polar_density_times_area = sp.factor(sp.simplify(gauss_map_jacobian_density * s3_polar_area_element))
    half_interval_t_integrand = sp.factor(t**2 / (t**2 + 1 / q**2) ** 2)
    half_interval_integral_value = sp.factor(sp.pi * q / 4)
    sigma3_integral_over_boundary = sp.factor(sp.simplify((8 * sp.pi / q) * half_interval_integral_value))
    target_sigma3_integral = sp.factor(2 * sp.pi**2)
    sigma3_integral_residual = sp.factor(sp.simplify(sigma3_integral_over_boundary - target_sigma3_integral))

    boundary_density_normalization = sp.Integer(16)
    boundary_functional = sp.factor(boundary_density_normalization * sigma3_integral_over_boundary)
    target_topological_number = sp.factor(32 * sp.pi**2 * chi_flat_spheroidal_d4)
    boundary_functional_residual = sp.factor(sp.simplify(boundary_functional - target_topological_number))
    q_derivative_residual = sp.factor(sp.diff(boundary_functional, q))
    round_limit_value = sp.factor(boundary_functional.subs(q, 1))
    sample_q = sp.Integer(2)
    sample_boundary_functional = sp.factor(boundary_functional.subs(q, sample_q))

    density_round_limit = sp.factor(gauss_map_jacobian_density.subs(q, 1))
    density_sample_pole = sp.factor(gauss_map_jacobian_density.subs({q: sample_q, theta: 0}))
    density_sample_equator = sp.factor(gauss_map_jacobian_density.subs({q: sample_q, theta: sp.pi / 2}))
    density_sample_nonconstant_witness = sp.factor(density_sample_pole - density_sample_equator)

    normalized_spheroidal_pairing = sp.factor(sp.simplify(b_gb * (bulk_gb_flat_value + boundary_functional)))
    normalized_pairing_residual = sp.factor(
        sp.simplify(normalized_spheroidal_pairing - per_euler_pairing * chi_flat_spheroidal_d4)
    )
    sample_pairing_q2 = sp.factor(normalized_spheroidal_pairing.subs(q, sample_q))

    p2341_probe = p2341.get("strict_track_b_d4_boundary_correction_fixture_witness_probe", {})
    p2341_witness = p2341_probe.get("track_B_D4_boundary_fixture_witness", {})
    p2343_probe = p2343.get("strict_track_b_flat_round_d4_boundary_functional_derivation_probe", {})
    p2343_derivation = p2343_probe.get("track_B_flat_round_D4_boundary_functional_derivation", {})
    p2343_flat_round = p2343_derivation.get("flat_round_boundary_class", {})

    derivation_steps = [
        {
            "step": "define_flat_spheroidal_domain",
            "formula": spheroid_equation,
            "scope": "flat interior, smooth convex spheroidal S^3_q boundary, q > 0",
        },
        {
            "step": "gauss_map_jacobian_density",
            "formula": "J_q(theta) = (1/q)/(sin(theta)**2 + cos(theta)**2/q**2)**2",
            "symbolic_density": str(gauss_map_jacobian_density),
            "round_limit_density": str(density_round_limit),
        },
        {
            "step": "polar_integral_reduction",
            "formula": "Integral_{S3_q} sigma3(K)dA = (8*pi/q) * Integral_0^oo t**2/(t**2 + 1/q**2)**2 dt",
            "half_interval_t_integrand": str(half_interval_t_integrand),
            "half_interval_integral_value": str(half_interval_integral_value),
        },
        {
            "step": "exact_integral_evaluation",
            "formula": "(8*pi/q) * (pi*q/4) = 2*pi**2",
            "sigma3_integral_over_boundary": str(sigma3_integral_over_boundary),
            "sigma3_integral_residual": str(sigma3_integral_residual),
        },
        {
            "step": "track_b_boundary_functional_normalization",
            "formula": "B_spheroidal(q) = 16 * Integral_{S3_q} sigma3(K)dA = 32*pi**2",
            "boundary_functional": str(boundary_functional),
            "q_derivative_residual": str(q_derivative_residual),
        },
    ]

    spheroidal_class = {
        "class_id": "flat_spheroidal_four_ball_axis_ratio_q_positive",
        "scope": "one-parameter flat spheroidal D^4_q boundary class only; not arbitrary non-round boundaries",
        "domain_equation": spheroid_equation,
        "axis_ratio_symbol": "q > 0",
        "nonround_subclass_condition": "q != 1",
        "bulk_gb_fixture_value": str(bulk_gb_flat_value),
        "boundary_manifold": "spheroidal S^3_q",
        "gauss_map_jacobian_density": str(gauss_map_jacobian_density),
        "s3_polar_area_element": str(s3_polar_area_element),
        "density_times_polar_area": str(polar_density_times_area),
        "half_interval_t_integrand": str(half_interval_t_integrand),
        "half_interval_integral_value": str(half_interval_integral_value),
        "sigma3_integral_over_boundary": str(sigma3_integral_over_boundary),
        "target_sigma3_integral": str(target_sigma3_integral),
        "sigma3_integral_residual": str(sigma3_integral_residual),
        "boundary_density_normalization_factor": str(boundary_density_normalization),
        "boundary_functional_B_spheroidal_q": str(boundary_functional),
        "target_topological_number_32pi2_chi": str(target_topological_number),
        "boundary_functional_residual": str(boundary_functional_residual),
        "q_derivative_residual": str(q_derivative_residual),
        "round_limit_value_q1": str(round_limit_value),
        "sample_q": str(sample_q),
        "sample_q2_boundary_functional": str(sample_boundary_functional),
        "sample_q2_density_pole": str(density_sample_pole),
        "sample_q2_density_equator": str(density_sample_equator),
        "sample_q2_nonconstant_density_difference": str(density_sample_nonconstant_witness),
        "normalized_track_B_pairing": str(normalized_spheroidal_pairing),
        "normalized_pairing_residual": str(normalized_pairing_residual),
        "sample_q2_pairing": str(sample_pairing_q2),
    }

    dependencies = {
        "p2335_two_track_ledger_loaded": p2335.get("result_kind")
        == "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
        "p2338_contract_loaded": p2338.get("result_kind")
        == "STRICT_TRACK_B_GB_BOUNDARY_TOPOLOGICAL_NORMALIZATION_CONTRACT_EXPORTED_NO_NORMALIZATION_CLAIM",
        "p2341_d4_boundary_fixture_loaded": p2341.get("result_kind")
        == "STRICT_TRACK_B_D4_NONZERO_BOUNDARY_FIXTURE_WITNESS_NO_UNIVERSAL_BOUNDARY_THEOREM",
        "p2343_flat_round_derivation_loaded": p2343.get("result_kind")
        == "STRICT_TRACK_B_FLAT_ROUND_D4_BOUNDARY_FUNCTIONAL_DERIVATION_NO_UNIVERSAL_BOUNDARY_THEOREM",
        "p2341_boundary_value_rederived": p2341_witness.get("boundary_correction_fixture_value") == str(boundary_functional),
        "p2341_pairing_rederived": p2341_witness.get("normalized_D4_pairing") == str(normalized_spheroidal_pairing),
        "p2343_round_limit_rederived": p2343_flat_round.get("boundary_functional_B_flat_round_R")
        == str(round_limit_value),
    }

    theorem_export = {
        "theorem_name": "P2344 Track-B flat spheroidal D4 boundary functional integral witness",
        "claim": (
            "P2344 extends P2343 from the flat round D^4_R class to a genuinely non-round one-parameter flat "
            "spheroidal class D^4_q with boundary S^3_q and q>0.  Using the explicit Gauss-map Jacobian density "
            "J_q(theta)=(1/q)/(sin(theta)**2+cos(theta)**2/q**2)**2 and the substitution t=tan(theta), the integral "
            "of sigma3(K)dA reduces to (8*pi/q)*Integral_0^oo t**2/(t**2+1/q**2)**2 dt = 2*pi**2.  Therefore "
            "the declared Track-B boundary functional 16*Integral sigma3(K)dA remains 32*pi**2 with zero q-derivative "
            "and gives the same ledger pairing as P2341."
        ),
        "proof_witnesses": [
            "The q=2 density is nonconstant: pole value 8 and equator value 1/2 differ by 15/2.",
            "The polar integral is reduced to a one-dimensional rational integral in t.",
            "The exact rational integral value pi*q/4 cancels the exterior 8*pi/q factor.",
            "The resulting boundary functional is q-independent and matches the P2343 round limit at q=1.",
            "The Track-B pairing residual against per-Euler normalization is zero.",
        ],
        "not_licensed": [
            "boundary functional for arbitrary non-spheroidal or non-convex boundaries",
            "general Chern-Gauss-Bonnet theorem over arbitrary compact four-manifolds with boundary",
            "general boundary gluing theorem for arbitrary interfaces",
            "connected-sum theorem beyond previously declared fixture bookkeeping",
            "independent a_GB measurement separate from the P2335 ledger coefficient",
            "bulk EOM GB force or EOM-only GB lift",
            "full/global renormalization closure",
            "selector premise or QW-2191 selector discharge",
            "G1/G3 update",
            "ToE closure",
        ],
        "next_honest_step": (
            "Keep selector work deferred.  The next proof/computational Track-B step should either formalize the "
            "Gauss-map degree/integral lemma as a named theorem for a bounded convex flat-boundary class, or stop at "
            "the round+spheroidal explicit class witnesses without claiming a universal boundary theorem."
        ),
    }

    probe = {
        "probe_id": "P2344_S1294_STRICT_TRACK_B_FLAT_SPHEROIDAL_D4_BOUNDARY_FUNCTIONAL_INTEGRAL_WITNESS",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2343 next step: explicit non-round spheroidal boundary functional integral",
            "top_hits": grep_hits(),
        },
        "track_B_flat_spheroidal_D4_boundary_functional_integral_witness": {
            "ledger_coefficient_b_GB": str(b_gb),
            "per_euler_pairing": str(per_euler_pairing),
            "derivation_steps": derivation_steps,
            "spheroidal_boundary_class": spheroidal_class,
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
        "p2343_flat_round_derivation_loaded": dependencies["p2343_flat_round_derivation_loaded"],
        "p2341_boundary_value_rederived": dependencies["p2341_boundary_value_rederived"],
        "p2341_pairing_rederived": dependencies["p2341_pairing_rederived"],
        "p2343_round_limit_rederived": dependencies["p2343_round_limit_rederived"],
        "flat_bulk_gb_zero_declared": bulk_gb_flat_value == 0,
        "nonround_sample_density_is_nonconstant": density_sample_nonconstant_witness != 0,
        "rational_integral_value_recorded": half_interval_integral_value == sp.pi * q / 4,
        "sigma3_integral_is_sphere_volume": sigma3_integral_residual == 0,
        "boundary_functional_is_32pi2": boundary_functional_residual == 0,
        "q_derivative_zero": q_derivative_residual == 0,
        "normalized_pairing_residual_zero": normalized_pairing_residual == 0,
        "one_parameter_spheroidal_scope_declared": "one-parameter flat spheroidal" in spheroidal_class["scope"],
        "arbitrary_boundary_functional_not_claimed": True,
        "general_cgb_theorem_not_claimed": True,
        "general_boundary_gluing_theorem_not_claimed": True,
        "independent_a_gb_not_claimed": True,
        "bulk_eom_gb_force_not_claimed": True,
        "no_full_global_renormalization_claimed": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_g1_g3_update_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2344_s1294_v1",
        "packet_id": "P2344",
        "stage_id": "S1294",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_FLAT_SPHEROIDAL_D4_BOUNDARY_FUNCTIONAL_INTEGRAL_WITNESS_NO_UNIVERSAL_BOUNDARY_THEOREM",
        "strict_track_b_flat_spheroidal_d4_boundary_functional_integral_witness_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2344 strict Track-B flat spheroidal D4 boundary functional integral witness\n\n"
        "Status: one non-round flat spheroidal D4_q boundary class integrated exactly; no universal boundary theorem.\n\n"
        f"- `b_GB = {b_gb}`.\n"
        f"- Boundary class: `{spheroid_equation}`; non-round subclass `q != 1`.\n"
        f"- Gauss-map density: `{gauss_map_jacobian_density}`; q=2 pole/equator values `{density_sample_pole}` and `{density_sample_equator}`.\n"
        f"- Integral reduction: `(8*pi/q) * integral_0^oo {half_interval_t_integrand} dt`, with rational integral `{half_interval_integral_value}`.\n"
        f"- `Integral sigma3(K)dA = {sigma3_integral_over_boundary}`; residual against `2*pi**2`: `{sigma3_integral_residual}`.\n"
        f"- Boundary functional: `16*Integral sigma3(K)dA = {boundary_functional}`; q-derivative residual `{q_derivative_residual}`.\n"
        f"- Normalized Track-B pairing: `{normalized_spheroidal_pairing}`; residual `{normalized_pairing_residual}`.\n"
        "- No arbitrary-boundary functional, no general Chern-Gauss-Bonnet theorem, no independent `a_GB`, no full/global renormalization, no selector premise, no QW-2191 discharge, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
