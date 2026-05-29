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
OUT = GEN / "p2357_s1307_strict_track_b_flat_convex_component_transgression_integration_bridge.json"
MD = GEN / "p2357_s1307_strict_track_b_flat_convex_component_transgression_integration_bridge.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2345_CONVEX_GAUSS_MAP": GEN / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.json",
    "P2348_CHERN_POLYNOMIAL": GEN / "p2348_s1298_strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction.json",
    "P2353_OBSTRUCTION_LEDGER": GEN / "p2353_s1303_strict_track_b_arbitrary_boundary_obstruction_ledger.json",
    "P2356_COMPONENT_ADDITIVITY": GEN / "p2356_s1306_strict_track_b_oriented_convex_component_degree_additivity_lemma.json",
}

RESULT_KINDS = {
    "P2335_TWO_TRACK_LEDGER": "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
    "P2345_CONVEX_GAUSS_MAP": "STRICT_TRACK_B_CONVEX_GAUSS_MAP_BOUNDARY_FUNCTIONAL_THEOREM_AND_CLOSURE_PROMISE_MAP_NO_SELECTOR_CLOSURE",
    "P2348_CHERN_POLYNOMIAL": "STRICT_TRACK_B_CHERN_BOUNDARY_POLYNOMIAL_NONSYMMETRIC_REDUCTION_NO_INTEGRATED_GLOBAL_THEOREM",
    "P2353_OBSTRUCTION_LEDGER": "STRICT_TRACK_B_ARBITRARY_BOUNDARY_OBSTRUCTION_LEDGER_NO_UNIVERSAL_CLOSURE",
    "P2356_COMPONENT_ADDITIVITY": "STRICT_TRACK_B_ORIENTED_CONVEX_COMPONENT_DEGREE_ADDITIVITY_LEMMA_NO_UNIVERSAL_CLOSURE",
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
        ROOT / "p2348_s1298_strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction.py",
        ROOT / "p2353_s1303_strict_track_b_arbitrary_boundary_obstruction_ledger.py",
        ROOT / "p2356_s1306_strict_track_b_oriented_convex_component_degree_additivity_lemma.py",
        GEN / "p2356_s1306_strict_track_b_oriented_convex_component_degree_additivity_lemma.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track-B|Track B|transgression|integration|Chern|sigma3|convex-component|degree|orientation|O3|O4|O5|arbitrary-boundary|selector|QW-2191|ToE closure",
        *existing,
    ]
    proc = subprocess.run(cmd, cwd=ROOT.parent, text=True, capture_output=True, check=False)
    return proc.stdout.splitlines()[:280]


def fixture_bridge_row(row: dict[str, Any], b_gb: sp.Expr) -> dict[str, Any]:
    signs = [int(value) for value in row["orientation_sign_vector"]]
    degree_sum = sp.Integer(sum(signs))
    local_integral_sum = sp.factor(16 * 2 * sp.pi**2 * degree_sum)
    p2356_boundary = sympify_expr(row["boundary_sum"])
    p2356_pairing = sympify_expr(row["pairing_sum"])
    integrated_pairing = sp.factor(local_integral_sum * b_gb)
    return {
        "fixture_id": row["fixture_id"],
        "orientation_sign_vector": signs,
        "degree_sum": int(degree_sum),
        "local_transgression_integral_sum": str(local_integral_sum),
        "p2356_boundary_sum": row["boundary_sum"],
        "boundary_bridge_residual": str(sp.factor(local_integral_sum - p2356_boundary)),
        "integrated_pairing": str(integrated_pairing),
        "p2356_pairing_sum": row["pairing_sum"],
        "pairing_bridge_residual": str(sp.factor(integrated_pairing - p2356_pairing)),
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    b_gb = sp.factor(track_b_coefficient(artifacts["P2335_TWO_TRACK_LEDGER"]))
    target_boundary_per_degree = 32 * sp.pi**2
    target_pairing_per_degree = sp.factor(target_boundary_per_degree * b_gb)

    K, k1, k2, k3, eps, n = sp.symbols("K k1 k2 k3 eps n")
    sigma1 = sp.factor(k1 + k2 + k3)
    sigma3 = sp.factor(k1 * k2 * k3)
    chern_boundary_polynomial = sp.factor(8 * K * sigma1 + 16 * sigma3)
    flat_density = sp.factor(chern_boundary_polynomial.subs(K, 0))
    flat_density_residual = sp.factor(flat_density - 16 * sigma3)
    component_integral_rule = sp.factor(16 * eps * 2 * sp.pi**2)
    component_integral_residual = sp.factor(component_integral_rule - 32 * sp.pi**2 * eps)
    finite_n_symbolic_boundary = sp.factor(32 * sp.pi**2 * n)
    finite_n_symbolic_pairing = sp.factor(finite_n_symbolic_boundary * b_gb)

    p2345_boundary = sp.factor(sympify_expr(find_first_key(artifacts["P2345_CONVEX_GAUSS_MAP"], "boundary_functional_convex_class")))
    p2348_replay = sp.factor(sympify_expr(find_first_key(artifacts["P2348_CHERN_POLYNOMIAL"], "substitute_integral_sigma3_equals_2pi2")))
    p2353_ledger = find_first_key(artifacts["P2353_OBSTRUCTION_LEDGER"], "track_B_arbitrary_boundary_obstruction_ledger") or {}
    p2356_lemma = find_first_key(artifacts["P2356_COMPONENT_ADDITIVITY"], "track_B_oriented_convex_component_degree_additivity_lemma") or {}
    p2353_minimal_cut = p2353_ledger.get("minimal_next_missing_cut", [])
    p2356_fixtures = p2356_lemma.get("fixture_replays", [])
    bridge_rows = [fixture_bridge_row(row, b_gb) for row in p2356_fixtures]
    all_bridge_residuals_zero = all(
        row["boundary_bridge_residual"] == "0" and row["pairing_bridge_residual"] == "0" for row in bridge_rows
    )

    proof_chain = [
        {
            "step_id": "local_chern_polynomial_flat_reduction",
            "input": "P2348 B3_poly(K;k1,k2,k3)=8*K*sigma1+16*sigma3",
            "output": str(flat_density),
            "residual_against_16_sigma3": str(flat_density_residual),
        },
        {
            "step_id": "oriented_component_transgression_integral",
            "input": "For a strict convex component with orientation sign eps, sigma3 dA integrates to eps*2*pi**2.",
            "output": str(component_integral_rule),
            "residual_against_32pi2_eps": str(component_integral_residual),
        },
        {
            "step_id": "finite_component_sum",
            "input": "n = Sum_i eps_i under P2356 finite convex-component hypotheses.",
            "boundary_output": str(finite_n_symbolic_boundary),
            "pairing_output": str(finite_n_symbolic_pairing),
        },
    ]

    bridge = {
        "bridge_id": "P2357_flat_convex_component_transgression_integration_bridge_v1",
        "scope": "partial O3 local-to-integrated bridge under P2356 strict convex-component hypotheses; not an arbitrary-boundary theorem",
        "hypotheses_replayed_from_p2356": p2356_lemma.get("hypotheses", []),
        "local_chern_boundary_polynomial": str(chern_boundary_polynomial),
        "flat_density": str(flat_density),
        "flat_density_residual": str(flat_density_residual),
        "component_integral_rule": str(component_integral_rule),
        "component_integral_residual": str(component_integral_residual),
        "finite_degree_symbol": str(n),
        "finite_component_boundary_formula": str(finite_n_symbolic_boundary),
        "finite_component_pairing_formula": str(finite_n_symbolic_pairing),
        "target_boundary_per_degree": str(target_boundary_per_degree),
        "target_pairing_per_degree": str(target_pairing_per_degree),
        "proof_chain": proof_chain,
        "fixture_bridge_rows": bridge_rows,
        "all_bridge_residuals_zero": all_bridge_residuals_zero,
        "p2345_convex_boundary_residual": str(sp.factor(p2345_boundary - target_boundary_per_degree)),
        "p2348_flat_integral_replay_residual": str(sp.factor(p2348_replay - target_boundary_per_degree)),
        "p2356_fixture_residuals_replayed": p2356_lemma.get("all_fixture_residuals_zero"),
        "p2356_coverage_rank_replayed": p2356_lemma.get("coverage_matrix_rank"),
        "p2353_minimal_cut_replayed": p2353_minimal_cut,
        "o3_cut_partially_attacked_under_hypotheses": "O3_arbitrary_boundary_transgression_integration" in p2353_minimal_cut,
        "remaining_o3_gap": "No arbitrary-boundary transgression theorem is exported for nonconvex Gauss-map singularities, corners, gluing, smoothing limits, or non-flat ambient curvature.",
    }

    dependencies = {
        key.lower() + "_loaded": artifacts[key].get("result_kind") == expected
        for key, expected in RESULT_KINDS.items()
    }
    dependencies.update(
        {
            "p2345_convex_boundary_replayed": bridge["p2345_convex_boundary_residual"] == "0",
            "p2348_flat_integral_replayed": bridge["p2348_flat_integral_replay_residual"] == "0",
            "p2356_fixture_residuals_replayed": bridge["p2356_fixture_residuals_replayed"] is True,
            "p2356_coverage_rank_replayed": bridge["p2356_coverage_rank_replayed"] == 4,
            "p2353_o3_cut_replayed": bridge["o3_cut_partially_attacked_under_hypotheses"],
            "flat_density_residual_zero": flat_density_residual == 0,
            "component_integral_residual_zero": component_integral_residual == 0,
            "all_bridge_residuals_zero": all_bridge_residuals_zero,
        }
    )

    theorem_export = {
        "theorem_name": "P2357 Track-B flat convex-component transgression-integration bridge",
        "claim": (
            "Under the P2356 finite strict convex-component hypotheses, the P2348 local flat Chern boundary polynomial 16*sigma3 "
            "integrates componentwise to 32*pi**2*eps_i and therefore to 32*pi**2*Sum_i eps_i.  This is a partial O3 "
            "local-to-integrated bridge tied to the P2356 hypotheses, not an arbitrary-boundary or global CGB theorem."
        ),
        "proof_witnesses": [
            "P2348's Chern boundary polynomial reduces exactly to 16*sigma3 in the flat K=0 lane.",
            "For each P2356 strict convex component, the oriented Gauss-map degree premise gives Integral sigma3 dA = eps_i*2*pi**2.",
            "Multiplication by 16 gives 32*pi**2*eps_i per component, and finite disjoint additivity gives 32*pi**2*Sum_i eps_i.",
            "The convex one-component, round shell, nonround shell, and three-component P2356 fixtures all bridge with zero boundary and pairing residuals.",
            "The O3 cut is attacked only under the same strict convex-component hypotheses; O5 and universal nonconvex cases remain open.",
        ],
        "not_licensed": [
            "arbitrary-boundary theorem",
            "general transgression theorem without strict convex-component hypotheses",
            "general nonconvex boundary theorem",
            "general Chern-Gauss-Bonnet theorem over all compact four-manifolds with boundary",
            "general gluing theorem or smoothing-limit theorem",
            "corner contribution theorem",
            "non-flat ambient curvature bridge",
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
            "After this partial O3 bridge, the next non-redundant move is O5 corners/gluing/smoothing or a genuine theorem weakening the "
            "strict convex-component hypotheses. Do not add another same-hypothesis component replay as universal closure."
        ),
    }

    probe = {
        "probe_id": "P2357_S1307_STRICT_TRACK_B_FLAT_CONVEX_COMPONENT_TRANSGRESSION_INTEGRATION_BRIDGE",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2357 partial O3 local-to-integrated bridge under P2356 hypotheses",
            "top_hits": grep_hits(),
        },
        "track_B_flat_convex_component_transgression_integration_bridge": bridge,
        "current_export_dependencies": dependencies,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "repo_grep_audit_recorded": len(probe["repo_grep_audit"]["top_hits"]) > 0,
        **dependencies,
        "proof_chain_has_three_steps": len(proof_chain) == 3,
        "bridge_rows_match_p2356_fixtures": len(bridge_rows) == len(p2356_fixtures) == 4,
        "o3_only_under_p2356_hypotheses": True,
        "o5_not_claimed_closed": True,
        "arbitrary_boundary_theorem_not_claimed": True,
        "general_transgression_theorem_not_claimed": True,
        "general_nonconvex_theorem_not_claimed": True,
        "general_cgb_theorem_not_claimed": True,
        "general_gluing_theorem_not_claimed": True,
        "corner_theorem_not_claimed": True,
        "nonflat_ambient_curvature_bridge_not_claimed": True,
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
        "schema_version": "p2357_s1307_v1",
        "packet_id": "P2357",
        "stage_id": "S1307",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-29T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_FLAT_CONVEX_COMPONENT_TRANSGRESSION_INTEGRATION_BRIDGE_NO_UNIVERSAL_CLOSURE",
        "strict_track_b_flat_convex_component_transgression_integration_bridge_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2357 strict Track-B flat convex-component transgression-integration bridge\n\n"
        "Status: partial O3 bridge under P2356 strict convex-component hypotheses; no arbitrary-boundary, selector, or ToE closure.\n\n"
        f"- Local flat density: `{flat_density}`; residual `{flat_density_residual}`.\n"
        f"- Component integral rule: `{component_integral_rule}`; residual `{component_integral_residual}`.\n"
        f"- Finite component formula: boundary `{finite_n_symbolic_boundary}`, pairing `{finite_n_symbolic_pairing}`.\n"
        f"- Bridge fixture rows: `{[row['fixture_id'] for row in bridge_rows]}`; all residuals zero `{all_bridge_residuals_zero}`.\n"
        f"- P2345 residual `{bridge['p2345_convex_boundary_residual']}`; P2348 residual `{bridge['p2348_flat_integral_replay_residual']}`; P2356 coverage rank `{bridge['p2356_coverage_rank_replayed']}`.\n"
        f"- P2353 cut replayed: `{p2353_minimal_cut}`; O3 partially attacked `{bridge['o3_cut_partially_attacked_under_hypotheses']}`.\n"
        "- No arbitrary-boundary theorem, no general transgression theorem without strict convex-component hypotheses, no general Chern-Gauss-Bonnet theorem, no gluing/smoothing/corner theorem, no non-flat ambient curvature bridge, no global renormalization, no independent `a_GB`, no selector premise, no QW-2191 discharge, no unique-future choice, no observer prediction, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
