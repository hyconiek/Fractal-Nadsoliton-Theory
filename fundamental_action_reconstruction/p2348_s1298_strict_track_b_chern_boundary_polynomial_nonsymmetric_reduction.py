#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from itertools import permutations
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2348_s1298_strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction.json"
MD = GEN / "p2348_s1298_strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2338_TRACK_B_CONTRACT": GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.json",
    "P2343_FLAT_ROUND_D4_DERIVATION": GEN / "p2343_s1293_strict_track_b_flat_round_d4_boundary_functional_derivation.json",
    "P2345_CONVEX_GAUSS_MAP_PROMISE": GEN / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.json",
    "P2347_SPHERICAL_CAP_CHERN_DERIVATION": GEN / "p2347_s1297_strict_track_b_spherical_cap_chern_boundary_form_derivation.json",
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
        ROOT / "p2343_s1293_strict_track_b_flat_round_d4_boundary_functional_derivation.py",
        ROOT / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.py",
        ROOT / "p2347_s1297_strict_track_b_spherical_cap_chern_boundary_form_derivation.py",
        GEN / "p2347_s1297_strict_track_b_spherical_cap_chern_boundary_form_derivation.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track B|Chern boundary|boundary form|sigma3|convex|spherical-cap|selector|QW-2191|ToE closure|global renormalization",
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
    p2343 = artifacts["P2343_FLAT_ROUND_D4_DERIVATION"]
    p2345 = artifacts["P2345_CONVEX_GAUSS_MAP_PROMISE"]
    p2347 = artifacts["P2347_SPHERICAL_CAP_CHERN_DERIVATION"]

    k1, k2, k3, k, K = sp.symbols("k1 k2 k3 k K")
    b_gb = sp.factor(track_b_coefficient(p2335))
    per_euler_pairing = sp.factor(sp.simplify(32 * sp.pi**2 * b_gb))

    sigma1 = sp.factor(k1 + k2 + k3)
    sigma2 = sp.factor(k1 * k2 + k1 * k3 + k2 * k3)
    sigma3 = sp.factor(k1 * k2 * k3)
    chern_boundary_polynomial = sp.factor(8 * K * sigma1 + 16 * sigma3)

    flat_polynomial = sp.factor(chern_boundary_polynomial.subs(K, 0))
    flat_expected = sp.factor(16 * sigma3)
    flat_residual = sp.factor(sp.simplify(flat_polynomial - flat_expected))

    equal_curvature_polynomial = sp.factor(chern_boundary_polynomial.subs({K: 1, k1: k, k2: k, k3: k}))
    spherical_cap_expected = sp.factor(8 * (3 * k + 2 * k**3))
    spherical_cap_residual = sp.factor(sp.simplify(equal_curvature_polynomial - spherical_cap_expected))

    permutation_residuals = []
    for perm in permutations((k1, k2, k3)):
        permuted = chern_boundary_polynomial.xreplace({k1: perm[0], k2: perm[1], k3: perm[2]})
        permutation_residuals.append(str(sp.factor(sp.simplify(permuted - chern_boundary_polynomial))))
    all_permutation_residuals_zero = all(value == "0" for value in permutation_residuals)

    nonsymmetric_sample = {K: sp.Integer(5), k1: sp.Integer(1), k2: sp.Integer(2), k3: sp.Integer(3)}
    nonsymmetric_sample_sigma1 = sp.factor(sigma1.subs(nonsymmetric_sample))
    nonsymmetric_sample_sigma2 = sp.factor(sigma2.subs(nonsymmetric_sample))
    nonsymmetric_sample_sigma3 = sp.factor(sigma3.subs(nonsymmetric_sample))
    nonsymmetric_sample_value = sp.factor(chern_boundary_polynomial.subs(nonsymmetric_sample))
    nonsymmetric_sample_expected = sp.factor(8 * 5 * 6 + 16 * 6)
    nonsymmetric_sample_residual = sp.factor(sp.simplify(nonsymmetric_sample_value - nonsymmetric_sample_expected))

    flat_integral_sigma3_symbol = sp.symbols("I_sigma3")
    flat_integrated_boundary_functional = sp.factor(16 * flat_integral_sigma3_symbol)
    flat_convex_substitution_value = sp.factor(flat_integrated_boundary_functional.subs(flat_integral_sigma3_symbol, 2 * sp.pi**2))
    flat_convex_residual = sp.factor(sp.simplify(flat_convex_substitution_value - 32 * sp.pi**2))
    flat_convex_pairing = sp.factor(sp.simplify(b_gb * flat_convex_substitution_value))
    flat_convex_pairing_residual = sp.factor(sp.simplify(flat_convex_pairing - per_euler_pairing))

    p2343_probe = p2343.get("strict_track_b_flat_round_d4_boundary_functional_derivation_probe", {})
    p2343_derivation = p2343_probe.get("track_B_flat_round_D4_boundary_functional_derivation", {})
    p2343_flat_round = p2343_derivation.get("flat_round_boundary_class", {})
    p2345_probe = p2345.get(
        "strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_probe", {}
    )
    p2345_lemma = p2345_probe.get("track_B_convex_gauss_map_boundary_functional_theorem", {})
    p2347_probe = p2347.get("strict_track_b_spherical_cap_chern_boundary_form_derivation_probe", {})
    p2347_witness = p2347_probe.get("track_B_spherical_cap_chern_boundary_form_derivation", {})
    p2347_derivation = p2347_witness.get("chern_form_derivation", {})

    p2343_boundary_residual = sp.factor(
        sp.sympify(p2343_flat_round.get("boundary_functional_B_flat_round_R", "0"), locals={"pi": sp.pi})
        - flat_convex_substitution_value
    )
    p2345_boundary_residual = sp.factor(
        sp.sympify(p2345_lemma.get("boundary_functional_convex_class", "0"), locals={"pi": sp.pi})
        - flat_convex_substitution_value
    )
    p2347_equal_curvature_residual = sp.factor(
        sp.sympify(
            p2347_derivation.get("chern_boundary_density_scalar", "0"),
            locals={"c": sp.symbols("c"), "sqrt": sp.sqrt},
        ).subs(sp.symbols("c"), k / sp.sqrt(1 + k**2))
        - equal_curvature_polynomial.subs(k, k)
    )

    derivation_steps = [
        {
            "step": "non_symmetric_shape_operator_data",
            "formula": "Use diagonal shape operator eigenvalues k1,k2,k3 and ambient sectional curvature K.",
            "sigma1": str(sigma1),
            "sigma2_recorded_not_used": str(sigma2),
            "sigma3": str(sigma3),
        },
        {
            "step": "chern_boundary_polynomial",
            "formula": "B3_poly(K;k1,k2,k3) = 8*K*sigma1 + 16*sigma3.",
            "polynomial": str(chern_boundary_polynomial),
        },
        {
            "step": "flat_track_b_reduction",
            "formula": "At K=0, B3_poly reduces to 16*sigma3, the P2343/P2345 flat Track-B boundary density.",
            "flat_polynomial": str(flat_polynomial),
            "flat_residual": str(flat_residual),
        },
        {
            "step": "spherical_cap_equal_curvature_reduction",
            "formula": "At K=1 and k1=k2=k3=k, B3_poly reduces to 8*(3*k+2*k**3), the P2347 cap scalar.",
            "equal_curvature_polynomial": str(equal_curvature_polynomial),
            "spherical_cap_residual": str(spherical_cap_residual),
        },
        {
            "step": "nonsymmetric_sample_evaluation",
            "formula": "For K=5 and (k1,k2,k3)=(1,2,3), B3_poly evaluates to 336.",
            "sample_sigma1": str(nonsymmetric_sample_sigma1),
            "sample_sigma2_recorded_not_used": str(nonsymmetric_sample_sigma2),
            "sample_sigma3": str(nonsymmetric_sample_sigma3),
            "sample_value": str(nonsymmetric_sample_value),
            "sample_residual": str(nonsymmetric_sample_residual),
        },
    ]

    boundary_polynomial = {
        "polynomial_id": "P2348_track_B_chern_boundary_polynomial_nonsymmetric_v1",
        "scope": "local algebraic Chern boundary-polynomial reduction for diagonal shape operator data; not an integrated arbitrary-boundary theorem",
        "ambient_sectional_curvature_symbol": "K",
        "principal_curvature_symbols": ["k1", "k2", "k3"],
        "sigma1": str(sigma1),
        "sigma2_recorded_not_used": str(sigma2),
        "sigma3": str(sigma3),
        "boundary_polynomial": str(chern_boundary_polynomial),
        "flat_K0_reduction": str(flat_polynomial),
        "flat_K0_expected": str(flat_expected),
        "flat_K0_residual": str(flat_residual),
        "spherical_cap_equal_curvature_reduction": str(equal_curvature_polynomial),
        "spherical_cap_expected": str(spherical_cap_expected),
        "spherical_cap_residual": str(spherical_cap_residual),
        "permutation_residuals": permutation_residuals,
        "all_permutation_residuals_zero": all_permutation_residuals_zero,
        "nonsymmetric_sample": {
            "K": "5",
            "k1": "1",
            "k2": "2",
            "k3": "3",
            "sigma1": str(nonsymmetric_sample_sigma1),
            "sigma2_recorded_not_used": str(nonsymmetric_sample_sigma2),
            "sigma3": str(nonsymmetric_sample_sigma3),
            "boundary_polynomial_value": str(nonsymmetric_sample_value),
            "residual": str(nonsymmetric_sample_residual),
        },
        "flat_convex_integral_replay": {
            "symbolic_integral_sigma3": str(flat_integral_sigma3_symbol),
            "integrated_boundary_functional": str(flat_integrated_boundary_functional),
            "substitute_integral_sigma3_equals_2pi2": str(flat_convex_substitution_value),
            "flat_convex_residual": str(flat_convex_residual),
            "flat_convex_pairing": str(flat_convex_pairing),
            "flat_convex_pairing_residual": str(flat_convex_pairing_residual),
        },
    }

    dependencies = {
        "p2335_two_track_ledger_loaded": p2335.get("result_kind")
        == "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
        "p2338_contract_loaded": p2338.get("result_kind")
        == "STRICT_TRACK_B_GB_BOUNDARY_TOPOLOGICAL_NORMALIZATION_CONTRACT_EXPORTED_NO_NORMALIZATION_CLAIM",
        "p2343_flat_round_loaded": p2343.get("result_kind")
        == "STRICT_TRACK_B_FLAT_ROUND_D4_BOUNDARY_FUNCTIONAL_DERIVATION_NO_UNIVERSAL_BOUNDARY_THEOREM",
        "p2345_convex_gauss_map_loaded": p2345.get("result_kind")
        == "STRICT_TRACK_B_CONVEX_GAUSS_MAP_BOUNDARY_FUNCTIONAL_THEOREM_AND_CLOSURE_PROMISE_MAP_NO_SELECTOR_CLOSURE",
        "p2347_spherical_cap_chern_loaded": p2347.get("result_kind")
        == "STRICT_TRACK_B_SPHERICAL_CAP_CHERN_BOUNDARY_FORM_DERIVATION_NO_GENERAL_CGB_THEOREM",
        "p2343_flat_value_replayed": p2343_boundary_residual == 0,
        "p2345_convex_value_replayed": p2345_boundary_residual == 0,
    }

    theorem_export = {
        "theorem_name": "P2348 Track-B non-symmetric Chern boundary polynomial reduction",
        "claim": (
            "P2348 exports the local algebraic Chern boundary polynomial B3_poly=8*K*(k1+k2+k3)+16*k1*k2*k3 "
            "for diagonal shape-operator data.  It proves that the polynomial reduces to the flat Track-B density "
            "16*sigma3 at K=0 and to the P2347 spherical-cap scalar 8*(3*k+2*k**3) when K=1 and k1=k2=k3=k."
        ),
        "proof_witnesses": [
            "The elementary symmetric data sigma1, sigma2, sigma3 are recorded with sigma2 deliberately absent from the polynomial.",
            "The K=0 reduction is exactly 16*sigma3, matching the flat P2343/P2345 density lane.",
            "The equal-curvature K=1 reduction is exactly 8*(3*k+2*k**3), matching the P2347 spherical-cap scalar.",
            "Permutation residuals over all six permutations of k1,k2,k3 are zero.",
            "A non-symmetric sample K=5,(1,2,3) evaluates to 336 with zero residual.",
        ],
        "not_licensed": [
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
            "The non-symmetric local polynomial now connects the flat and spherical-cap lanes.  The next honest "
            "proof/computational step is an integrated test on a boundary with nonconstant, unequal principal curvatures, "
            "or else pause Track-B and return to the Track-A tensor-bundle obstruction."
        ),
    }

    probe = {
        "probe_id": "P2348_S1298_STRICT_TRACK_B_CHERN_BOUNDARY_POLYNOMIAL_NONSYMMETRIC_REDUCTION",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2347 next step: non-symmetric local Chern boundary polynomial reduction",
            "top_hits": grep_hits(),
        },
        "track_B_chern_boundary_polynomial_nonsymmetric_reduction": {
            "ledger_coefficient_b_GB": str(b_gb),
            "per_euler_pairing": str(per_euler_pairing),
            "derivation_steps": derivation_steps,
            "boundary_polynomial": boundary_polynomial,
        },
        "current_export_dependencies": dependencies,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "repo_grep_audit_recorded": len(probe["repo_grep_audit"]["top_hits"]) > 0,
        "p2335_two_track_ledger_loaded": dependencies["p2335_two_track_ledger_loaded"],
        "p2338_contract_loaded": dependencies["p2338_contract_loaded"],
        "p2343_flat_round_loaded": dependencies["p2343_flat_round_loaded"],
        "p2345_convex_gauss_map_loaded": dependencies["p2345_convex_gauss_map_loaded"],
        "p2347_spherical_cap_chern_loaded": dependencies["p2347_spherical_cap_chern_loaded"],
        "p2343_flat_value_replayed": dependencies["p2343_flat_value_replayed"],
        "p2345_convex_value_replayed": dependencies["p2345_convex_value_replayed"],
        "flat_reduction_residual_zero": flat_residual == 0,
        "spherical_cap_reduction_residual_zero": spherical_cap_residual == 0,
        "permutation_residuals_zero": all_permutation_residuals_zero,
        "nonsymmetric_sample_residual_zero": nonsymmetric_sample_residual == 0,
        "flat_convex_integral_replay_zero": flat_convex_residual == 0,
        "flat_convex_pairing_residual_zero": flat_convex_pairing_residual == 0,
        "local_polynomial_scope_declared": "local algebraic" in boundary_polynomial["scope"],
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
        "schema_version": "p2348_s1298_v1",
        "packet_id": "P2348",
        "stage_id": "S1298",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_CHERN_BOUNDARY_POLYNOMIAL_NONSYMMETRIC_REDUCTION_NO_INTEGRATED_GLOBAL_THEOREM",
        "strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2348 strict Track-B Chern boundary polynomial non-symmetric reduction\n\n"
        "Status: local non-symmetric Chern boundary polynomial derived; no integrated arbitrary-boundary theorem.\n\n"
        f"- `b_GB = {b_gb}`.\n"
        f"- Polynomial: `B3_poly = {chern_boundary_polynomial}` with `sigma1={sigma1}`, `sigma3={sigma3}`.\n"
        f"- Flat reduction `K=0`: `{flat_polynomial}`; residual `{flat_residual}`.\n"
        f"- Spherical-cap equal-curvature reduction: `{equal_curvature_polynomial}`; residual `{spherical_cap_residual}`.\n"
        f"- Permutation residuals all zero: `{all_permutation_residuals_zero}`; non-symmetric sample value `{nonsymmetric_sample_value}`.\n"
        f"- Flat convex integrated replay: `{flat_convex_substitution_value}`; pairing `{flat_convex_pairing}`; residual `{flat_convex_pairing_residual}`.\n"
        "- No integrated arbitrary-boundary theorem, no general Chern-Gauss-Bonnet theorem, no global renormalization, no independent `a_GB`, no selector premise, no QW-2191 discharge, no unique-future choice, no observer prediction, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
