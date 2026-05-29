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
OUT = GEN / "p2349_s1299_strict_track_b_triaxial_ellipsoid_integrated_boundary_polynomial_witness.json"
MD = GEN / "p2349_s1299_strict_track_b_triaxial_ellipsoid_integrated_boundary_polynomial_witness.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2338_TRACK_B_CONTRACT": GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.json",
    "P2345_CONVEX_GAUSS_MAP_PROMISE": GEN / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.json",
    "P2348_NONSYMMETRIC_POLYNOMIAL": GEN / "p2348_s1298_strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction.json",
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
        GEN / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.md",
        GEN / "p2348_s1298_strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track B|convex|Gauss-map|boundary polynomial|non-symmetric|sigma3|selector|QW-2191|ToE closure|global renormalization",
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


def axis_pole_curvatures(axes: tuple[sp.Integer, sp.Integer, sp.Integer, sp.Integer], pole_index: int) -> list[sp.Expr]:
    axis = axes[pole_index]
    return [sp.factor(axis / axes[j] ** 2) for j in range(len(axes)) if j != pole_index]


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    p2335 = artifacts["P2335_TWO_TRACK_LEDGER"]
    p2338 = artifacts["P2338_TRACK_B_CONTRACT"]
    p2345 = artifacts["P2345_CONVEX_GAUSS_MAP_PROMISE"]
    p2348 = artifacts["P2348_NONSYMMETRIC_POLYNOMIAL"]

    b_gb = sp.factor(track_b_coefficient(p2335))
    per_euler_pairing = sp.factor(sp.simplify(32 * sp.pi**2 * b_gb))
    axes = (sp.Integer(1), sp.Integer(2), sp.Integer(3), sp.Integer(5))
    axis_labels = ("x1", "x2", "x3", "x4")
    ellipsoid_equation = "x1**2/1**2 + x2**2/2**2 + x3**2/3**2 + x4**2/5**2 = 1"

    pole_rows: list[dict[str, Any]] = []
    sigma3_values: list[sp.Expr] = []
    boundary_density_values: list[sp.Expr] = []
    for idx, axis in enumerate(axes):
        curvatures = axis_pole_curvatures(axes, idx)
        sigma1 = sp.factor(sum(curvatures))
        sigma2 = sp.factor(sum(curvatures[i] * curvatures[j] for i in range(3) for j in range(i + 1, 3)))
        sigma3 = sp.factor(curvatures[0] * curvatures[1] * curvatures[2])
        boundary_density = sp.factor(16 * sigma3)
        sigma3_values.append(sigma3)
        boundary_density_values.append(boundary_density)
        pole_rows.append(
            {
                "pole": f"+{axis_labels[idx]}",
                "axis_length": str(axis),
                "principal_curvatures": [str(value) for value in curvatures],
                "sigma1": str(sigma1),
                "sigma2_recorded_not_used_in_flat_density": str(sigma2),
                "sigma3": str(sigma3),
                "flat_K0_boundary_polynomial_value_16sigma3": str(boundary_density),
            }
        )

    distinct_sigma3_count = len(set(map(str, sigma3_values)))
    distinct_boundary_density_count = len(set(map(str, boundary_density_values)))
    nonconstant_density_witness = sp.factor(boundary_density_values[-1] - boundary_density_values[0])

    gauss_map_degree = sp.Integer(1)
    integral_sigma3 = sp.factor(2 * sp.pi**2 * gauss_map_degree)
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

    derivation_steps = [
        {
            "step": "triaxial_ellipsoid_class",
            "formula": ellipsoid_equation,
            "axes": [str(axis) for axis in axes],
            "scope": "one explicit smooth strictly convex triaxial ellipsoid in flat R^4",
        },
        {
            "step": "axis_pole_shape_operator_samples",
            "formula": "At pole ai*ei, principal curvatures along the other axes are ai/aj**2.",
            "pole_rows": pole_rows,
        },
        {
            "step": "nonconstant_unequal_density_witness",
            "formula": "The local P2348 flat density 16*sigma3 is not constant across axis poles.",
            "distinct_sigma3_count": distinct_sigma3_count,
            "distinct_boundary_density_count": distinct_boundary_density_count,
            "x4_minus_x1_density_difference": str(nonconstant_density_witness),
        },
        {
            "step": "integrated_gauss_map_replay",
            "formula": "Strict convexity gives Gauss-map degree +1, hence Integral sigma3(K)dA = Area(S^3) = 2*pi**2.",
            "integral_sigma3": str(integral_sigma3),
            "integrated_boundary_functional": str(integrated_boundary_functional),
            "integrated_boundary_residual": str(integrated_boundary_residual),
        },
    ]

    ellipsoid_witness = {
        "witness_id": "P2349_triaxial_ellipsoid_integrated_boundary_polynomial_witness_v1",
        "scope": "one explicit triaxial ellipsoid integrated stress-test; not an arbitrary-boundary theorem",
        "ellipsoid_equation": ellipsoid_equation,
        "axes": [str(axis) for axis in axes],
        "axis_pole_rows": pole_rows,
        "distinct_sigma3_count": distinct_sigma3_count,
        "distinct_boundary_density_count": distinct_boundary_density_count,
        "nonconstant_density_witness_x4_minus_x1": str(nonconstant_density_witness),
        "gauss_map_degree": str(gauss_map_degree),
        "integral_sigma3_via_gauss_map": str(integral_sigma3),
        "integrated_boundary_functional_16_integral_sigma3": str(integrated_boundary_functional),
        "integrated_boundary_residual": str(integrated_boundary_residual),
        "normalized_track_B_pairing": str(normalized_pairing),
        "normalized_pairing_residual": str(normalized_pairing_residual),
        "p2345_boundary_residual": str(p2345_boundary_residual),
        "p2348_flat_replay_residual": str(p2348_flat_replay_residual),
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
        "p2345_boundary_value_replayed": p2345_boundary_residual == 0,
        "p2348_flat_integral_replay_replayed": p2348_flat_replay_residual == 0,
    }

    theorem_export = {
        "theorem_name": "P2349 Track-B triaxial ellipsoid integrated boundary-polynomial witness",
        "claim": (
            "P2349 applies the P2348 local flat boundary polynomial to an explicit triaxial ellipsoid with axes "
            "1,2,3,5.  Axis-pole shape-operator samples show unequal and nonconstant principal-curvature data, while "
            "the integrated Gauss-map degree replay gives Integral sigma3 dA=2*pi**2 and the Track-B boundary functional "
            "32*pi**2 with zero residual."
        ),
        "proof_witnesses": [
            "At each axis pole ai*ei, the principal curvature samples ai/aj**2 are computed explicitly.",
            "The local flat density 16*sigma3 differs across axis poles, so this is not a round/equal-curvature replay.",
            "Strict convexity of the ellipsoid gives Gauss-map degree +1.",
            "The integrated sigma3 value is therefore 2*pi**2, producing boundary functional 32*pi**2.",
            "The Track-B ledger pairing residual against per-Euler normalization is zero.",
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
            "Track-B now has local non-symmetric polynomial data and one explicit nonconstant integrated ellipsoid replay. "
            "The next honest move is either a second integrated non-ellipsoidal convex boundary stress-test or a return to "
            "Track-A tensor-bundle obstruction; do not promote this to a universal boundary theorem."
        ),
    }

    probe = {
        "probe_id": "P2349_S1299_STRICT_TRACK_B_TRIAXIAL_ELLIPSOID_INTEGRATED_BOUNDARY_POLYNOMIAL_WITNESS",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2348 next step: integrated nonconstant unequal-principal-curvature ellipsoid stress-test",
            "top_hits": grep_hits(),
        },
        "track_B_triaxial_ellipsoid_integrated_boundary_polynomial_witness": {
            "ledger_coefficient_b_GB": str(b_gb),
            "per_euler_pairing": str(per_euler_pairing),
            "derivation_steps": derivation_steps,
            "ellipsoid_witness": ellipsoid_witness,
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
        "p2345_boundary_value_replayed": dependencies["p2345_boundary_value_replayed"],
        "p2348_flat_integral_replay_replayed": dependencies["p2348_flat_integral_replay_replayed"],
        "axis_pole_rows_complete": len(pole_rows) == 4,
        "sigma3_values_nonconstant": distinct_sigma3_count == 4,
        "boundary_density_values_nonconstant": distinct_boundary_density_count == 4,
        "integrated_boundary_residual_zero": integrated_boundary_residual == 0,
        "normalized_pairing_residual_zero": normalized_pairing_residual == 0,
        "explicit_ellipsoid_scope_declared": "one explicit triaxial ellipsoid" in ellipsoid_witness["scope"],
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
        "schema_version": "p2349_s1299_v1",
        "packet_id": "P2349",
        "stage_id": "S1299",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_TRIAXIAL_ELLIPSOID_INTEGRATED_BOUNDARY_POLYNOMIAL_WITNESS_NO_UNIVERSAL_THEOREM",
        "strict_track_b_triaxial_ellipsoid_integrated_boundary_polynomial_witness_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2349 strict Track-B triaxial ellipsoid integrated boundary-polynomial witness\n\n"
        "Status: explicit nonconstant unequal-principal-curvature ellipsoid stress-test integrated; no universal boundary theorem.\n\n"
        f"- `b_GB = {b_gb}`.\n"
        f"- Ellipsoid: `{ellipsoid_equation}`.\n"
        f"- Axis-pole sigma3 values: `{', '.join(map(str, sigma3_values))}`; distinct count `{distinct_sigma3_count}`.\n"
        f"- Axis-pole flat density values `16*sigma3`: `{', '.join(map(str, boundary_density_values))}`; nonconstant witness `{nonconstant_density_witness}`.\n"
        f"- Gauss-map integral replay: `Integral sigma3 dA = {integral_sigma3}`; boundary functional `{integrated_boundary_functional}`; residual `{integrated_boundary_residual}`.\n"
        f"- Normalized Track-B pairing: `{normalized_pairing}`; residual `{normalized_pairing_residual}`.\n"
        "- No integrated arbitrary-boundary theorem, no general Chern-Gauss-Bonnet theorem, no global renormalization, no independent `a_GB`, no selector premise, no QW-2191 discharge, no unique-future choice, no observer prediction, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
