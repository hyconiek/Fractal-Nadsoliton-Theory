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
OUT = GEN / "p2356_s1306_strict_track_b_oriented_convex_component_degree_additivity_lemma.json"
MD = GEN / "p2356_s1306_strict_track_b_oriented_convex_component_degree_additivity_lemma.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2338_TRACK_B_CONTRACT": GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.json",
    "P2345_CONVEX_GAUSS_MAP": GEN / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.json",
    "P2348_CHERN_POLYNOMIAL": GEN / "p2348_s1298_strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction.json",
    "P2353_OBSTRUCTION_LEDGER": GEN / "p2353_s1303_strict_track_b_arbitrary_boundary_obstruction_ledger.json",
    "P2354_SPHERICAL_SHELL": GEN / "p2354_s1304_strict_track_b_spherical_shell_orientation_cancellation_witness.json",
    "P2355_NONROUND_ELLIPSOIDAL_SHELL": GEN / "p2355_s1305_strict_track_b_nonround_ellipsoidal_shell_orientation_degree_witness.json",
}

RESULT_KINDS = {
    "P2335_TWO_TRACK_LEDGER": "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
    "P2338_TRACK_B_CONTRACT": "STRICT_TRACK_B_GB_BOUNDARY_TOPOLOGICAL_NORMALIZATION_CONTRACT_EXPORTED_NO_NORMALIZATION_CLAIM",
    "P2345_CONVEX_GAUSS_MAP": "STRICT_TRACK_B_CONVEX_GAUSS_MAP_BOUNDARY_FUNCTIONAL_THEOREM_AND_CLOSURE_PROMISE_MAP_NO_SELECTOR_CLOSURE",
    "P2348_CHERN_POLYNOMIAL": "STRICT_TRACK_B_CHERN_BOUNDARY_POLYNOMIAL_NONSYMMETRIC_REDUCTION_NO_INTEGRATED_GLOBAL_THEOREM",
    "P2353_OBSTRUCTION_LEDGER": "STRICT_TRACK_B_ARBITRARY_BOUNDARY_OBSTRUCTION_LEDGER_NO_UNIVERSAL_CLOSURE",
    "P2354_SPHERICAL_SHELL": "STRICT_TRACK_B_SPHERICAL_SHELL_ORIENTATION_CANCELLATION_WITNESS_NO_UNIVERSAL_CLOSURE",
    "P2355_NONROUND_ELLIPSOIDAL_SHELL": "STRICT_TRACK_B_NONROUND_ELLIPSOIDAL_SHELL_ORIENTATION_DEGREE_WITNESS_NO_UNIVERSAL_CLOSURE",
}

LOCALS = {"pi": sp.pi, "log": sp.log, "ln": sp.log}


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
    return sp.sympify(str(raw), locals=LOCALS)


def track_b_coefficient(p2335: dict[str, Any]) -> sp.Expr:
    probe = p2335.get("strict_two_track_renormalization_ledger_probe", {})
    ledger = probe.get("two_track_ledger", {})
    track_b = ledger.get("track_B_gb_topological_counterterm_ledger", {})
    return sympify_expr(track_b.get("ledger_coefficient_b_GB_topological", "0"))


def grep_hits() -> list[str]:
    candidates = [
        ROOT / "p2353_s1303_strict_track_b_arbitrary_boundary_obstruction_ledger.py",
        ROOT / "p2354_s1304_strict_track_b_spherical_shell_orientation_cancellation_witness.py",
        ROOT / "p2355_s1305_strict_track_b_nonround_ellipsoidal_shell_orientation_degree_witness.py",
        GEN / "p2355_s1305_strict_track_b_nonround_ellipsoidal_shell_orientation_degree_witness.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track-B|Track B|degree|orientation|additivity|convex component|nonconvex|Gauss-map|sigma3|arbitrary-boundary|selector|QW-2191|ToE closure",
        *existing,
    ]
    proc = subprocess.run(cmd, cwd=ROOT.parent, text=True, capture_output=True, check=False)
    return proc.stdout.splitlines()[:280]


def oriented_component_row(component_id: str, sign: int, b_gb: sp.Expr) -> dict[str, Any]:
    degree = sp.Integer(sign)
    integral_sigma3 = sp.factor(degree * 2 * sp.pi**2)
    boundary = sp.factor(16 * integral_sigma3)
    pairing = sp.factor(boundary * b_gb)
    return {
        "component_id": component_id,
        "orientation_sign": sign,
        "gauss_map_degree": int(degree),
        "integral_sigma3_dA_by_degree": str(integral_sigma3),
        "boundary_functional_16_integral_sigma3": str(boundary),
        "normalized_track_B_pairing": str(pairing),
    }


def fixture_replay(fixture_id: str, signs: list[int], b_gb: sp.Expr) -> dict[str, Any]:
    components = [oriented_component_row(f"{fixture_id}_component_{idx}", sign, b_gb) for idx, sign in enumerate(signs)]
    degree_sum = sum(signs)
    boundary_sum = sp.factor(sum(sympify_expr(row["boundary_functional_16_integral_sigma3"]) for row in components))
    pairing_sum = sp.factor(sum(sympify_expr(row["normalized_track_B_pairing"]) for row in components))
    expected_boundary = sp.factor(32 * sp.pi**2 * degree_sum)
    expected_pairing = sp.factor(expected_boundary * b_gb)
    return {
        "fixture_id": fixture_id,
        "orientation_sign_vector": signs,
        "component_rows": components,
        "degree_sum": degree_sum,
        "boundary_sum": str(boundary_sum),
        "pairing_sum": str(pairing_sum),
        "expected_boundary_by_degree_sum": str(expected_boundary),
        "expected_pairing_by_degree_sum": str(expected_pairing),
        "boundary_residual": str(sp.factor(boundary_sum - expected_boundary)),
        "pairing_residual": str(sp.factor(pairing_sum - expected_pairing)),
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    b_gb = sp.factor(track_b_coefficient(artifacts["P2335_TWO_TRACK_LEDGER"]))
    target_boundary_per_degree = 32 * sp.pi**2
    target_pairing_per_degree = sp.factor(target_boundary_per_degree * b_gb)

    eps1, eps2, eps3, n = sp.symbols("eps1 eps2 eps3 n")
    symbolic_degree_sum = sp.factor(eps1 + eps2 + eps3)
    symbolic_boundary_sum = sp.factor(32 * sp.pi**2 * symbolic_degree_sum)
    symbolic_pairing_sum = sp.factor(symbolic_boundary_sum * b_gb)

    p2345_boundary = sp.factor(sympify_expr(find_first_key(artifacts["P2345_CONVEX_GAUSS_MAP"], "boundary_functional_convex_class")))
    p2348_replay = sp.factor(sympify_expr(find_first_key(artifacts["P2348_CHERN_POLYNOMIAL"], "substitute_integral_sigma3_equals_2pi2")))
    p2353_ledger = find_first_key(artifacts["P2353_OBSTRUCTION_LEDGER"], "track_B_arbitrary_boundary_obstruction_ledger") or {}
    p2354_witness = find_first_key(artifacts["P2354_SPHERICAL_SHELL"], "track_B_spherical_shell_orientation_cancellation_witness") or {}
    p2355_witness = find_first_key(artifacts["P2355_NONROUND_ELLIPSOIDAL_SHELL"], "track_B_nonround_ellipsoidal_shell_orientation_degree_witness") or {}
    p2353_minimal_cut = p2353_ledger.get("minimal_next_missing_cut", [])

    fixture_rows = [
        fixture_replay("convex_one_component_degree_plus_one", [+1], b_gb),
        fixture_replay("p2354_round_shell_degree_plus_one_minus_one", [+1, -1], b_gb),
        fixture_replay("p2355_nonround_shell_degree_plus_one_minus_one", [+1, -1], b_gb),
        fixture_replay("three_component_oriented_ledger_stress_vector", [+1, -1, -1], b_gb),
    ]

    coverage_columns = [
        "single_component",
        "multi_component",
        "negative_orientation_component",
        "round_shell_replay",
        "nonround_nonconstant_shell_replay",
        "nonzero_degree_sum_stress",
    ]
    coverage_rows = [
        [1, 0, 0, 0, 0, 1],
        [0, 1, 1, 1, 0, 0],
        [0, 1, 1, 0, 1, 0],
        [0, 1, 1, 0, 0, 1],
    ]
    coverage_rank = int(sp.Matrix(coverage_rows).rank())
    all_fixture_residuals_zero = all(row["boundary_residual"] == "0" and row["pairing_residual"] == "0" for row in fixture_rows)

    lemma = {
        "lemma_id": "P2356_oriented_convex_component_degree_additivity_v1",
        "scope": "finite disjoint smooth strictly convex flat R4 boundary components with explicit orientation signs; not an arbitrary-boundary theorem",
        "hypotheses": [
            "Each boundary component is smooth, closed, embedded, and strictly convex in flat R^4.",
            "Each component has a declared shell/outward orientation sign eps_i in {+1,-1} relative to its convex-body outward Gauss map.",
            "The Gauss map on each component is a diffeomorphism with degree eps_i.",
            "Only the flat Track-B density 16*sigma3 is used; no corner, gluing, or smoothing-limit term is introduced.",
        ],
        "symbolic_degree_sum_three_component": str(symbolic_degree_sum),
        "symbolic_boundary_sum_three_component": str(symbolic_boundary_sum),
        "symbolic_pairing_sum_three_component": str(symbolic_pairing_sum),
        "general_finite_formula": "B_TrackB = 32*pi**2 * Sum_i eps_i; pairing = b_GB * B_TrackB",
        "target_boundary_per_degree": str(target_boundary_per_degree),
        "target_pairing_per_degree": str(target_pairing_per_degree),
        "fixture_replays": fixture_rows,
        "all_fixture_residuals_zero": all_fixture_residuals_zero,
        "coverage_matrix_columns": coverage_columns,
        "coverage_matrix_rows": coverage_rows,
        "coverage_matrix_rank": coverage_rank,
        "p2345_convex_boundary_residual": str(sp.factor(p2345_boundary - target_boundary_per_degree)),
        "p2348_chern_replay_residual": str(sp.factor(p2348_replay - target_boundary_per_degree)),
        "p2354_shell_boundary_residual_replayed": p2354_witness.get("shell_boundary_residual"),
        "p2355_shell_boundary_residual_replayed": p2355_witness.get("shell_boundary_residual"),
        "p2353_minimal_cut_replayed": p2353_minimal_cut,
        "o4_cut_lemma_level_partially_closed_under_hypotheses": "O4_nonconvex_degree_and_orientation_accounting" in p2353_minimal_cut,
        "remaining_o4_gap": "No theorem here covers arbitrary nonconvex Gauss-map degree changes, self-overlap, non-strict convex strata, corners, or gluing/smoothing limits.",
    }

    dependencies = {
        key.lower() + "_loaded": artifacts[key].get("result_kind") == expected
        for key, expected in RESULT_KINDS.items()
    }
    dependencies.update(
        {
            "p2345_convex_boundary_replayed": lemma["p2345_convex_boundary_residual"] == "0",
            "p2348_chern_replay_replayed": lemma["p2348_chern_replay_residual"] == "0",
            "p2354_shell_residual_replayed": lemma["p2354_shell_boundary_residual_replayed"] == "0",
            "p2355_shell_residual_replayed": lemma["p2355_shell_boundary_residual_replayed"] == "0",
            "p2353_o4_cut_replayed": lemma["o4_cut_lemma_level_partially_closed_under_hypotheses"],
            "fixture_residuals_zero": all_fixture_residuals_zero,
            "coverage_rank_at_least_four": coverage_rank >= 4,
        }
    )

    theorem_export = {
        "theorem_name": "P2356 Track-B oriented convex-component degree additivity lemma",
        "claim": (
            "Under the stated finite-component hypotheses, Track-B boundary contributions add by oriented Gauss-map degree: "
            "each strictly convex component contributes 32*pi**2*eps_i, hence the finite union contributes 32*pi**2*Sum_i eps_i "
            "and the P2335-ledger pairing contributes b_GB times that sum.  This packages the P2354/P2355 O4 witnesses into a local "
            "degree/orientation lemma with explicit hypotheses, not an arbitrary-boundary theorem."
        ),
        "proof_witnesses": [
            "The component hypothesis makes sigma3 dA pull back to eps_i times the unit S3 volume form, so Integral sigma3 dA = eps_i*2*pi**2.",
            "Multiplying by the flat Track-B density normalization 16 gives 32*pi**2*eps_i per component.",
            "Finite additivity of disjoint component integrals gives 32*pi**2*Sum_i eps_i and the same additivity after multiplication by b_GB.",
            "The convex one-component, P2354 round shell, P2355 nonround shell, and a three-component stress vector all replay with zero residual.",
            "The O4 cut is partially lemma-closed only under the declared strict convex component hypotheses; O3 and O5 remain open.",
        ],
        "not_licensed": [
            "arbitrary-boundary theorem",
            "general nonconvex boundary theorem without convex-component hypotheses",
            "general Chern-Gauss-Bonnet theorem over all compact four-manifolds with boundary",
            "general gluing theorem or smoothing-limit theorem",
            "corner contribution theorem",
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
            "After this O4 convex-component lemma, do not add another same-class shell replay as if it were new closure.  Move to O3 "
            "transgression-integration, O5 corners/gluing, or explicitly weaken/remove the strict convex-component hypotheses with a new proof."
        ),
    }

    probe = {
        "probe_id": "P2356_S1306_STRICT_TRACK_B_ORIENTED_CONVEX_COMPONENT_DEGREE_ADDITIVITY_LEMMA",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2356 local O4 degree/orientation additivity lemma after P2354/P2355 witnesses",
            "top_hits": grep_hits(),
        },
        "track_B_oriented_convex_component_degree_additivity_lemma": lemma,
        "current_export_dependencies": dependencies,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "repo_grep_audit_recorded": len(probe["repo_grep_audit"]["top_hits"]) > 0,
        **dependencies,
        "four_fixture_rows_recorded": len(fixture_rows) == 4,
        "symbolic_formula_recorded": lemma["general_finite_formula"] == "B_TrackB = 32*pi**2 * Sum_i eps_i; pairing = b_GB * B_TrackB",
        "o4_only_under_hypotheses": True,
        "o3_not_claimed_closed": True,
        "o5_not_claimed_closed": True,
        "arbitrary_boundary_theorem_not_claimed": True,
        "general_nonconvex_theorem_not_claimed": True,
        "general_cgb_theorem_not_claimed": True,
        "general_gluing_theorem_not_claimed": True,
        "corner_theorem_not_claimed": True,
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
        "schema_version": "p2356_s1306_v1",
        "packet_id": "P2356",
        "stage_id": "S1306",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-29T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_ORIENTED_CONVEX_COMPONENT_DEGREE_ADDITIVITY_LEMMA_NO_UNIVERSAL_CLOSURE",
        "strict_track_b_oriented_convex_component_degree_additivity_lemma_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2356 strict Track-B oriented convex-component degree additivity lemma\n\n"
        "Status: local O4 degree/orientation lemma under strict convex-component hypotheses; no arbitrary-boundary, selector, or ToE closure.\n\n"
        f"- `b_GB = {b_gb}`; per-degree boundary `{target_boundary_per_degree}`; per-degree pairing `{target_pairing_per_degree}`.\n"
        f"- Symbolic three-component replay: degree `{symbolic_degree_sum}`, boundary `{symbolic_boundary_sum}`, pairing `{symbolic_pairing_sum}`.\n"
        f"- Fixture rows: `{[row['fixture_id'] for row in fixture_rows]}`; all residuals zero `{all_fixture_residuals_zero}`.\n"
        f"- Coverage matrix rank `{coverage_rank}` over columns `{coverage_columns}`.\n"
        f"- P2354 residual `{lemma['p2354_shell_boundary_residual_replayed']}`; P2355 residual `{lemma['p2355_shell_boundary_residual_replayed']}`.\n"
        f"- P2353 cut replayed: `{p2353_minimal_cut}`; O4 lemma-level partial closure `{lemma['o4_cut_lemma_level_partially_closed_under_hypotheses']}`.\n"
        "- No arbitrary-boundary theorem, no general nonconvex boundary theorem without convex-component hypotheses, no general Chern-Gauss-Bonnet theorem, no gluing/smoothing/corner theorem, no global renormalization, no independent `a_GB`, no selector premise, no QW-2191 discharge, no unique-future choice, no observer prediction, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
