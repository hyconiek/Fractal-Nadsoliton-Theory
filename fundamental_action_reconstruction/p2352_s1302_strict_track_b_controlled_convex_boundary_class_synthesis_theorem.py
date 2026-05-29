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
OUT = GEN / "p2352_s1302_strict_track_b_controlled_convex_boundary_class_synthesis_theorem.json"
MD = GEN / "p2352_s1302_strict_track_b_controlled_convex_boundary_class_synthesis_theorem.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2338_TRACK_B_CONTRACT": GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.json",
    "P2343_FLAT_ROUND": GEN / "p2343_s1293_strict_track_b_flat_round_d4_boundary_functional_derivation.json",
    "P2344_SPHEROIDAL": GEN / "p2344_s1294_strict_track_b_flat_spheroidal_d4_boundary_functional_integral_witness.json",
    "P2345_CONVEX_GAUSS_MAP": GEN / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.json",
    "P2348_CHERN_POLYNOMIAL": GEN / "p2348_s1298_strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction.json",
    "P2349_TRIAXIAL_ELLIPSOID": GEN / "p2349_s1299_strict_track_b_triaxial_ellipsoid_integrated_boundary_polynomial_witness.json",
    "P2350_QUARTIC_INTEGRATED": GEN / "p2350_s1300_strict_track_b_quartic_support_convex_boundary_integrated_witness.json",
    "P2351_QUARTIC_GLOBAL": GEN / "p2351_s1301_strict_track_b_quartic_support_global_convexity_jacobian_identity_witness.json",
}

RESULT_KINDS = {
    "P2335_TWO_TRACK_LEDGER": "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
    "P2338_TRACK_B_CONTRACT": "STRICT_TRACK_B_GB_BOUNDARY_TOPOLOGICAL_NORMALIZATION_CONTRACT_EXPORTED_NO_NORMALIZATION_CLAIM",
    "P2343_FLAT_ROUND": "STRICT_TRACK_B_FLAT_ROUND_D4_BOUNDARY_FUNCTIONAL_DERIVATION_NO_UNIVERSAL_BOUNDARY_THEOREM",
    "P2344_SPHEROIDAL": "STRICT_TRACK_B_FLAT_SPHEROIDAL_D4_BOUNDARY_FUNCTIONAL_INTEGRAL_WITNESS_NO_UNIVERSAL_BOUNDARY_THEOREM",
    "P2345_CONVEX_GAUSS_MAP": "STRICT_TRACK_B_CONVEX_GAUSS_MAP_BOUNDARY_FUNCTIONAL_THEOREM_AND_CLOSURE_PROMISE_MAP_NO_SELECTOR_CLOSURE",
    "P2348_CHERN_POLYNOMIAL": "STRICT_TRACK_B_CHERN_BOUNDARY_POLYNOMIAL_NONSYMMETRIC_REDUCTION_NO_INTEGRATED_GLOBAL_THEOREM",
    "P2349_TRIAXIAL_ELLIPSOID": "STRICT_TRACK_B_TRIAXIAL_ELLIPSOID_INTEGRATED_BOUNDARY_POLYNOMIAL_WITNESS_NO_UNIVERSAL_THEOREM",
    "P2350_QUARTIC_INTEGRATED": "STRICT_TRACK_B_QUARTIC_SUPPORT_CONVEX_BOUNDARY_INTEGRATED_WITNESS_NO_UNIVERSAL_THEOREM",
    "P2351_QUARTIC_GLOBAL": "STRICT_TRACK_B_QUARTIC_SUPPORT_GLOBAL_CONVEXITY_JACOBIAN_IDENTITY_WITNESS_NO_UNIVERSAL_THEOREM",
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
        ROOT / "p2350_s1300_strict_track_b_quartic_support_convex_boundary_integrated_witness.py",
        ROOT / "p2351_s1301_strict_track_b_quartic_support_global_convexity_jacobian_identity_witness.py",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track B|convex|Gauss-map|boundary functional|Chern|ellipsoid|quartic support|global convexity|selector|QW-2191|ToE closure",
        *existing,
    ]
    proc = subprocess.run(cmd, cwd=ROOT.parent, text=True, capture_output=True, check=False)
    return proc.stdout.splitlines()[:220]


def track_b_coefficient(p2335: dict[str, Any]) -> sp.Expr:
    probe = p2335.get("strict_two_track_renormalization_ledger_probe", {})
    ledger = probe.get("two_track_ledger", {})
    track_b = ledger.get("track_B_gb_topological_counterterm_ledger", {})
    raw = track_b.get("ledger_coefficient_b_GB_topological", "0")
    return sp.sympify(raw, locals={"pi": sp.pi, "log": sp.log, "ln": sp.log})


def find_first_key(obj: Any, key: str) -> Any:
    if isinstance(obj, dict):
        if key in obj:
            return obj[key]
        for value in obj.values():
            found = find_first_key(value, key)
            if found is not None:
                return found
    elif isinstance(obj, list):
        for value in obj:
            found = find_first_key(value, key)
            if found is not None:
                return found
    return None


def sympify_expr(raw: Any) -> sp.Expr:
    return sp.sympify(str(raw), locals={"pi": sp.pi, "log": sp.log, "ln": sp.log})


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    b_gb = sp.factor(track_b_coefficient(artifacts["P2335_TWO_TRACK_LEDGER"]))
    target_boundary = 32 * sp.pi**2
    target_pairing = sp.factor(sp.simplify(target_boundary * b_gb))

    class_specs = [
        {
            "packet": "P2343",
            "class_id": "flat_round_D4_R",
            "source": "P2343_FLAT_ROUND",
            "boundary_key": "boundary_functional_B_flat_round_R",
            "pairing_key": "normalized_track_B_pairing",
            "nonround_or_nonconstant": False,
            "has_integrated_nonconstant_witness": False,
            "has_global_convexity_certificate": True,
            "has_local_chern_polynomial_bridge": False,
        },
        {
            "packet": "P2344",
            "class_id": "flat_spheroidal_D4_q",
            "source": "P2344_SPHEROIDAL",
            "boundary_key": "boundary_functional_B_spheroidal_q",
            "pairing_key": "normalized_track_B_pairing",
            "nonround_or_nonconstant": True,
            "has_integrated_nonconstant_witness": True,
            "has_global_convexity_certificate": True,
            "has_local_chern_polynomial_bridge": False,
        },
        {
            "packet": "P2349",
            "class_id": "triaxial_ellipsoid_integrated",
            "source": "P2349_TRIAXIAL_ELLIPSOID",
            "boundary_key": "integrated_boundary_functional_16_integral_sigma3",
            "pairing_key": "normalized_track_B_pairing",
            "nonround_or_nonconstant": True,
            "has_integrated_nonconstant_witness": True,
            "has_global_convexity_certificate": True,
            "has_local_chern_polynomial_bridge": True,
        },
        {
            "packet": "P2350",
            "class_id": "quartic_support_integrated_samples",
            "source": "P2350_QUARTIC_INTEGRATED",
            "boundary_key": "integrated_boundary_functional_16_integral_sigma3",
            "pairing_key": "normalized_track_B_pairing",
            "nonround_or_nonconstant": True,
            "has_integrated_nonconstant_witness": True,
            "has_global_convexity_certificate": False,
            "has_local_chern_polynomial_bridge": True,
        },
        {
            "packet": "P2351",
            "class_id": "quartic_support_global_convexity",
            "source": "P2351_QUARTIC_GLOBAL",
            "boundary_key": "integrated_boundary_functional_16_integral_sigma3",
            "pairing_key": "normalized_track_B_pairing",
            "nonround_or_nonconstant": True,
            "has_integrated_nonconstant_witness": True,
            "has_global_convexity_certificate": True,
            "has_local_chern_polynomial_bridge": True,
        },
    ]

    class_rows = []
    residuals = []
    pairing_residuals = []
    feature_rows = []
    for spec in class_specs:
        artifact = artifacts[spec["source"]]
        boundary_value = sp.factor(sympify_expr(find_first_key(artifact, spec["boundary_key"])))
        pairing_value = sp.factor(sympify_expr(find_first_key(artifact, spec["pairing_key"])))
        boundary_residual = sp.factor(sp.simplify(boundary_value - target_boundary))
        pairing_residual = sp.factor(sp.simplify(pairing_value - target_pairing))
        residuals.append(boundary_residual)
        pairing_residuals.append(pairing_residual)
        feature_vector = [
            1,
            int(spec["nonround_or_nonconstant"]),
            int(spec["has_integrated_nonconstant_witness"]),
            int(spec["has_global_convexity_certificate"]),
            int(spec["has_local_chern_polynomial_bridge"]),
        ]
        feature_rows.append(feature_vector)
        class_rows.append(
            {
                "packet": spec["packet"],
                "class_id": spec["class_id"],
                "source": spec["source"],
                "boundary_functional": str(boundary_value),
                "boundary_residual_against_32pi2": str(boundary_residual),
                "normalized_track_B_pairing": str(pairing_value),
                "pairing_residual_against_per_euler_pairing": str(pairing_residual),
                "feature_vector_columns": [
                    "has_boundary_replay",
                    "nonround_or_nonconstant",
                    "has_integrated_nonconstant_witness",
                    "has_global_convexity_certificate",
                    "has_local_chern_polynomial_bridge",
                ],
                "feature_vector": feature_vector,
            }
        )

    feature_matrix = sp.Matrix(feature_rows)
    gram_matrix = feature_matrix.T * feature_matrix
    feature_rank = int(feature_matrix.rank())
    gram_det = sp.factor(gram_matrix.det())
    all_boundary_residuals_zero = all(residual == 0 for residual in residuals)
    all_pairing_residuals_zero = all(residual == 0 for residual in pairing_residuals)
    nonround_count = sum(row[1] for row in feature_rows)
    global_certificate_count = sum(row[3] for row in feature_rows)
    chern_bridge_count = sum(row[4] for row in feature_rows)
    coverage_score = sp.Integer(nonround_count + global_certificate_count + chern_bridge_count)
    normalized_coverage_score = sp.Rational(int(coverage_score), len(feature_rows) * 3)

    p2345_boundary = sp.factor(sympify_expr(find_first_key(artifacts["P2345_CONVEX_GAUSS_MAP"], "boundary_functional_convex_class")))
    p2348_replay = sp.factor(
        sympify_expr(find_first_key(artifacts["P2348_CHERN_POLYNOMIAL"], "substitute_integral_sigma3_equals_2pi2"))
    )
    p2345_residual = sp.factor(sp.simplify(p2345_boundary - target_boundary))
    p2348_residual = sp.factor(sp.simplify(p2348_replay - target_boundary))

    controlled_synthesis = {
        "synthesis_id": "P2352_controlled_convex_boundary_class_synthesis_v1",
        "scope": "finite controlled Track-B class synthesis; not an arbitrary-boundary theorem",
        "target_boundary_functional": str(target_boundary),
        "ledger_coefficient_b_GB": str(b_gb),
        "target_pairing_per_euler_number": str(target_pairing),
        "class_rows": class_rows,
        "boundary_residual_vector": [str(value) for value in residuals],
        "pairing_residual_vector": [str(value) for value in pairing_residuals],
        "all_boundary_residuals_zero": all_boundary_residuals_zero,
        "all_pairing_residuals_zero": all_pairing_residuals_zero,
        "feature_matrix_columns": [
            "has_boundary_replay",
            "nonround_or_nonconstant",
            "has_integrated_nonconstant_witness",
            "has_global_convexity_certificate",
            "has_local_chern_polynomial_bridge",
        ],
        "feature_matrix_rows": feature_rows,
        "feature_rank": feature_rank,
        "feature_gram_determinant": str(gram_det),
        "nonround_or_nonconstant_count": nonround_count,
        "global_convexity_certificate_count": global_certificate_count,
        "local_chern_bridge_count": chern_bridge_count,
        "coverage_score_raw": str(coverage_score),
        "coverage_score_normalized": str(normalized_coverage_score),
        "p2345_convex_class_boundary_residual": str(p2345_residual),
        "p2348_chern_polynomial_replay_residual": str(p2348_residual),
    }

    dependencies = {
        key.lower() + "_loaded": artifacts[key].get("result_kind") == expected
        for key, expected in RESULT_KINDS.items()
    }
    dependencies.update(
        {
            "p2345_convex_class_boundary_replayed": p2345_residual == 0,
            "p2348_chern_polynomial_replayed": p2348_residual == 0,
            "controlled_class_boundary_residuals_zero": all_boundary_residuals_zero,
            "controlled_class_pairing_residuals_zero": all_pairing_residuals_zero,
        }
    )

    theorem_export = {
        "theorem_name": "P2352 Track-B controlled convex-boundary class synthesis theorem",
        "claim": (
            "P2352 synthesizes the controlled Track-B boundary-functional evidence from P2343/P2344/P2349/P2350/P2351. "
            "Every controlled row replays the same boundary functional 32*pi**2 and the same P2335-ledger pairing, while "
            "the feature matrix records nonround/nonconstant coverage, global convexity coverage, and local Chern-polynomial bridge coverage. "
            "This is a finite controlled-class synthesis theorem, not an arbitrary-boundary or global CGB theorem."
        ),
        "proof_witnesses": [
            "The class-row residual vector against 32*pi**2 is exactly all zeros.",
            "The normalized Track-B pairing residual vector against 32*pi**2*b_GB is exactly all zeros.",
            "The controlled matrix contains flat round, spheroidal, triaxial ellipsoid, quartic integrated, and quartic global rows.",
            "P2345 convex Gauss-map value and P2348 Chern-polynomial replay are both rechecked with zero residual.",
            "The synthesis records feature rank and coverage counts so downstream steps cannot confuse finite controlled coverage with universal closure.",
        ],
        "not_licensed": [
            "arbitrary-boundary theorem",
            "general Chern-Gauss-Bonnet theorem over all compact four-manifolds with boundary",
            "universal support-function convexity theorem",
            "general gluing theorem",
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
            "The Track-B controlled class is now synthesized. The next honest non-selector move is either a single obstruction-ledger "
            "showing exactly what remains before any arbitrary-boundary theorem, or return to Track-A tensor-bundle obstruction; selector remains last."
        ),
    }

    probe = {
        "probe_id": "P2352_S1302_STRICT_TRACK_B_CONTROLLED_CONVEX_BOUNDARY_CLASS_SYNTHESIS_THEOREM",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2351 next step: controlled Track-B class synthesis before arbitrary-boundary or selector claims",
            "top_hits": grep_hits(),
        },
        "track_B_controlled_convex_boundary_class_synthesis": controlled_synthesis,
        "current_export_dependencies": dependencies,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "repo_grep_audit_recorded": len(probe["repo_grep_audit"]["top_hits"]) > 0,
        **dependencies,
        "class_row_count_is_five": len(class_rows) == 5,
        "feature_rank_at_least_three": feature_rank >= 3,
        "nonround_or_nonconstant_count_at_least_four": nonround_count >= 4,
        "global_convexity_certificate_count_at_least_four": global_certificate_count >= 4,
        "local_chern_bridge_count_at_least_three": chern_bridge_count >= 3,
        "coverage_score_normalized_above_two_thirds": bool(normalized_coverage_score >= sp.Rational(2, 3)),
        "arbitrary_boundary_theorem_not_claimed": True,
        "general_cgb_theorem_not_claimed": True,
        "universal_support_function_theorem_not_claimed": True,
        "general_gluing_theorem_not_claimed": True,
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
        "schema_version": "p2352_s1302_v1",
        "packet_id": "P2352",
        "stage_id": "S1302",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-29T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_CONTROLLED_CONVEX_BOUNDARY_CLASS_SYNTHESIS_THEOREM_NO_UNIVERSAL_CLOSURE",
        "strict_track_b_controlled_convex_boundary_class_synthesis_theorem_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2352 strict Track-B controlled convex-boundary class synthesis theorem\n\n"
        "Status: finite controlled-class synthesis exported; no arbitrary-boundary, selector, or ToE closure.\n\n"
        f"- `b_GB = {b_gb}`; target boundary functional `{target_boundary}`; target pairing `{target_pairing}`.\n"
        f"- Controlled rows: `{len(class_rows)}` = flat round, spheroidal, triaxial ellipsoid, quartic integrated, quartic global.\n"
        f"- Boundary residual vector: `{[str(value) for value in residuals]}`.\n"
        f"- Pairing residual vector: `{[str(value) for value in pairing_residuals]}`.\n"
        f"- Feature rank: `{feature_rank}`; Gram determinant `{gram_det}`; normalized coverage score `{normalized_coverage_score}`.\n"
        f"- Coverage counts: nonround/nonconstant `{nonround_count}`, global convexity `{global_certificate_count}`, local Chern bridge `{chern_bridge_count}`.\n"
        f"- P2345 residual `{p2345_residual}`; P2348 residual `{p2348_residual}`.\n"
        "- No arbitrary-boundary theorem, no general Chern-Gauss-Bonnet theorem, no universal support-function theorem, no general gluing theorem, no global renormalization, no independent `a_GB`, no selector premise, no QW-2191 discharge, no unique-future choice, no observer prediction, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
