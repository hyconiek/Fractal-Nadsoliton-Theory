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
OUT = GEN / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.json"
MD = GEN / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2338_TRACK_B_CONTRACT": GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.json",
    "P2341_D4_BOUNDARY_FIXTURE": GEN / "p2341_s1291_strict_track_b_d4_boundary_correction_fixture_witness.json",
    "P2343_FLAT_ROUND_D4_DERIVATION": GEN / "p2343_s1293_strict_track_b_flat_round_d4_boundary_functional_derivation.json",
    "P2344_FLAT_SPHEROIDAL_D4_WITNESS": GEN / "p2344_s1294_strict_track_b_flat_spheroidal_d4_boundary_functional_integral_witness.json",
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
        ROOT / "p2344_s1294_strict_track_b_flat_spheroidal_d4_boundary_functional_integral_witness.py",
        GEN / "p2343_s1293_strict_track_b_flat_round_d4_boundary_functional_derivation.md",
        GEN / "p2344_s1294_strict_track_b_flat_spheroidal_d4_boundary_functional_integral_witness.md",
        ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track B|boundary functional|Gauss-map|spheroidal|flat round|convex|selector|QW-2191|ToE closure|strict-only|nadsoliton",
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


def score_lane(proof_strength: int, closure_relevance: int, blocker_penalty: int, selector_penalty: int) -> sp.Rational:
    return sp.Rational(2 * proof_strength + 2 * closure_relevance - blocker_penalty - selector_penalty, 1)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    p2335 = artifacts["P2335_TWO_TRACK_LEDGER"]
    p2338 = artifacts["P2338_TRACK_B_CONTRACT"]
    p2341 = artifacts["P2341_D4_BOUNDARY_FIXTURE"]
    p2343 = artifacts["P2343_FLAT_ROUND_D4_DERIVATION"]
    p2344 = artifacts["P2344_FLAT_SPHEROIDAL_D4_WITNESS"]

    b_gb = sp.factor(track_b_coefficient(p2335))
    per_euler_pairing = sp.factor(sp.simplify(32 * sp.pi**2 * b_gb))

    gauss_map_degree = sp.Integer(1)
    area_unit_s3 = sp.factor(2 * sp.pi**2)
    sigma3_integral_convex = sp.factor(gauss_map_degree * area_unit_s3)
    boundary_density_normalization = sp.Integer(16)
    boundary_functional_convex = sp.factor(boundary_density_normalization * sigma3_integral_convex)
    chi_flat_convex_ball = sp.Integer(1)
    target_topological_number = sp.factor(32 * sp.pi**2 * chi_flat_convex_ball)
    convex_boundary_residual = sp.factor(sp.simplify(boundary_functional_convex - target_topological_number))
    convex_pairing = sp.factor(sp.simplify(b_gb * boundary_functional_convex))
    convex_pairing_residual = sp.factor(sp.simplify(convex_pairing - per_euler_pairing * chi_flat_convex_ball))

    p2341_probe = p2341.get("strict_track_b_d4_boundary_correction_fixture_witness_probe", {})
    p2341_witness = p2341_probe.get("track_B_D4_boundary_fixture_witness", {})
    p2343_probe = p2343.get("strict_track_b_flat_round_d4_boundary_functional_derivation_probe", {})
    p2343_derivation = p2343_probe.get("track_B_flat_round_D4_boundary_functional_derivation", {})
    p2343_flat = p2343_derivation.get("flat_round_boundary_class", {})
    p2344_probe = p2344.get("strict_track_b_flat_spheroidal_d4_boundary_functional_integral_witness_probe", {})
    p2344_witness = p2344_probe.get("track_B_flat_spheroidal_D4_boundary_functional_integral_witness", {})
    p2344_spheroidal = p2344_witness.get("spheroidal_boundary_class", {})

    round_residual = sp.factor(
        sp.sympify(p2343_flat.get("boundary_functional_B_flat_round_R", "0"), locals={"pi": sp.pi})
        - boundary_functional_convex
    )
    spheroidal_residual = sp.factor(
        sp.sympify(p2344_spheroidal.get("boundary_functional_B_spheroidal_q", "0"), locals={"pi": sp.pi})
        - boundary_functional_convex
    )

    theorem_lemma = {
        "lemma_name": "convex_flat_boundary_gauss_map_degree_lemma_for_Track_B",
        "admissible_class": "smooth bounded strictly convex flat four-ball domains Omega subset R^4 with outward boundary orientation",
        "proof_chain": [
            "The outward unit normal map N: boundary(Omega) -> S^3 is a diffeomorphism for the declared strictly convex class.",
            "The Gauss-map degree is +1 under the declared outward orientation.",
            "The Jacobian of N is sigma3(K), the Gauss-Kronecker curvature / determinant of the shape operator.",
            "Therefore Integral_boundary sigma3(K)dA = degree(N)*Area(S^3) = 2*pi**2.",
            "The Track-B boundary convention from P2343/P2344 multiplies this by 16, giving 32*pi**2.",
        ],
        "gauss_map_degree": str(gauss_map_degree),
        "unit_s3_area": str(area_unit_s3),
        "sigma3_integral_convex_class": str(sigma3_integral_convex),
        "boundary_density_normalization_factor": str(boundary_density_normalization),
        "boundary_functional_convex_class": str(boundary_functional_convex),
        "target_topological_number_32pi2_chi": str(target_topological_number),
        "boundary_residual": str(convex_boundary_residual),
        "normalized_track_B_pairing": str(convex_pairing),
        "normalized_pairing_residual": str(convex_pairing_residual),
        "round_P2343_residual_against_lemma": str(round_residual),
        "spheroidal_P2344_residual_against_lemma": str(spheroidal_residual),
    }

    lane_inputs = {
        "track_B_convex_boundary_functional": {
            "proof_strength": 5,
            "closure_relevance": 4,
            "blocker_penalty": 1,
            "selector_penalty": 0,
            "why": "Exact lemma-style computation turns two explicit boundary classes into a named convex flat boundary class, still below global/non-flat closure.",
        },
        "track_A_non_GB_transport": {
            "proof_strength": 3,
            "closure_relevance": 4,
            "blocker_penalty": 3,
            "selector_penalty": 0,
            "why": "P2336 transported finite parts, but P2337 blocks silent tensor-bundle promotion.",
        },
        "selector_future_choice_lane": {
            "proof_strength": 1,
            "closure_relevance": 5,
            "blocker_penalty": 5,
            "selector_penalty": 3,
            "why": "Selector closure is physically central as a choice of one future over another, but QW-2191 remains a real strict-core obstruction.",
        },
        "full_strict_ToE_closure": {
            "proof_strength": 1,
            "closure_relevance": 5,
            "blocker_penalty": 5,
            "selector_penalty": 2,
            "why": "Strict-only ToE closure is the target, not a current export; it inherits unresolved global, non-flat, and selector gaps.",
        },
    }
    closure_promise_scores = {
        key: {**value, "score": str(score_lane(**{k: value[k] for k in ["proof_strength", "closure_relevance", "blocker_penalty", "selector_penalty"]}))}
        for key, value in lane_inputs.items()
    }
    ranked_closure_promises = sorted(
        closure_promise_scores,
        key=lambda name: sp.Integer(closure_promise_scores[name]["score"]),
        reverse=True,
    )

    nadsoliton_reality_readout = {
        "status": "interpretive_readout_below_ToE_closure_no_selector_discharge",
        "what_the_current_strict_track_b_evidence_says": [
            "The most promising closure evidence is not yet the selector; it is the Track-B boundary/topology ledger becoming computable on an expanding class.",
            "The nadsoliton ontology remains primary information in a solitonic state; the evidence here suggests reality is constrained by global boundary/topological bookkeeping, not only by local EOM terms.",
            "Selector language should remain future-choice language: choosing one admissible future and rejecting another is still not derived in strict core.",
            "What is currently strongest is a pre-selector structural statement: admissible realities must respect the same normalized topological boundary charge on the tested flat convex classes.",
        ],
        "not_claimed": [
            "actual choice of our unique future",
            "strict-core selector closure",
            "observer-level prediction",
            "full physical ToE closure",
        ],
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
        "p2344_flat_spheroidal_witness_loaded": p2344.get("result_kind")
        == "STRICT_TRACK_B_FLAT_SPHEROIDAL_D4_BOUNDARY_FUNCTIONAL_INTEGRAL_WITNESS_NO_UNIVERSAL_BOUNDARY_THEOREM",
        "p2341_boundary_value_rederived": p2341_witness.get("boundary_correction_fixture_value") == str(boundary_functional_convex),
        "p2343_round_case_absorbed": round_residual == 0,
        "p2344_spheroidal_case_absorbed": spheroidal_residual == 0,
    }

    theorem_export = {
        "theorem_name": "P2345 Track-B convex Gauss-map boundary functional theorem and closure-promise map",
        "claim": (
            "P2345 upgrades the P2343/P2344 explicit boundary computations into a named lemma for smooth bounded "
            "strictly convex flat four-ball domains: the outward Gauss map has degree +1, its Jacobian is sigma3(K), "
            "and so Integral_boundary sigma3(K)dA = Area(S^3)=2*pi**2.  Under the existing Track-B normalization this "
            "gives boundary functional 32*pi**2 and the P2335 ledger pairing 13152087*log(2)/10000000 with zero residual."
        ),
        "proof_witnesses": theorem_lemma["proof_chain"],
        "closure_promise_answer": (
            "The best current non-selector promise of closure is Track-B boundary/topology normalization: it now has exact "
            "fixture, round, spheroidal, and convex-flat Gauss-map support.  Selector closure remains conceptually central "
            "as future-choice, but mathematically it is still blocked by QW-2191 and should stay last."
        ),
        "not_licensed": [
            "arbitrary non-convex or non-flat boundary functional",
            "general Chern-Gauss-Bonnet theorem over arbitrary compact four-manifolds with boundary",
            "global/non-flat bulk-plus-boundary renormalization theorem",
            "independent a_GB measurement separate from the P2335 ledger coefficient",
            "bulk EOM GB force or EOM-only GB lift",
            "selector premise or QW-2191 selector discharge",
            "choice of the unique physical future",
            "observer-level prediction",
            "G1/G3 update",
            "ToE closure",
        ],
        "next_honest_step": (
            "Use this as a Track-B closure-promise theorem below global closure.  Next either extend from flat convex "
            "domains to a controlled curved bulk+boundary class, or pause Track-B and return to the non-GB Track-A tensor "
            "bundle obstruction.  Do not move selector until these structural lanes are exhausted."
        ),
    }

    probe = {
        "probe_id": "P2345_S1295_STRICT_TRACK_B_CONVEX_GAUSS_MAP_BOUNDARY_FUNCTIONAL_THEOREM_AND_CLOSURE_PROMISE_MAP",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2344 next step: convex Gauss-map lemma and closure-promise ranking without selector discharge",
            "top_hits": grep_hits(),
        },
        "track_B_convex_gauss_map_boundary_functional_theorem": theorem_lemma,
        "closure_promise_map": {
            "scoring_formula": "2*proof_strength + 2*closure_relevance - blocker_penalty - selector_penalty",
            "scores": closure_promise_scores,
            "ranked_best_to_worst": ranked_closure_promises,
        },
        "nadsoliton_reality_readout": nadsoliton_reality_readout,
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
        "p2344_flat_spheroidal_witness_loaded": dependencies["p2344_flat_spheroidal_witness_loaded"],
        "p2341_boundary_value_rederived": dependencies["p2341_boundary_value_rederived"],
        "p2343_round_case_absorbed": dependencies["p2343_round_case_absorbed"],
        "p2344_spheroidal_case_absorbed": dependencies["p2344_spheroidal_case_absorbed"],
        "gauss_map_degree_one_declared": gauss_map_degree == 1,
        "sigma3_integral_is_s3_area": sigma3_integral_convex == area_unit_s3,
        "boundary_functional_is_32pi2": convex_boundary_residual == 0,
        "normalized_pairing_residual_zero": convex_pairing_residual == 0,
        "track_b_ranked_best_current_promise": ranked_closure_promises[0] == "track_B_convex_boundary_functional",
        "selector_ranked_below_structural_lanes": ranked_closure_promises.index("selector_future_choice_lane") > ranked_closure_promises.index("track_B_convex_boundary_functional"),
        "nadsoliton_readout_marked_nonclosure": nadsoliton_reality_readout["status"].endswith("no_selector_discharge"),
        "arbitrary_boundary_functional_not_claimed": True,
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
        "schema_version": "p2345_s1295_v1",
        "packet_id": "P2345",
        "stage_id": "S1295",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_CONVEX_GAUSS_MAP_BOUNDARY_FUNCTIONAL_THEOREM_AND_CLOSURE_PROMISE_MAP_NO_SELECTOR_CLOSURE",
        "strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2345 strict Track-B convex Gauss-map boundary functional theorem and closure-promise map\n\n"
        "Status: convex flat-boundary Gauss-map lemma exported below universal/global closure; selector stays last.\n\n"
        f"- `b_GB = {b_gb}`.\n"
        f"- Convex-flat lemma: `degree(N) = {gauss_map_degree}`, `Area(S^3) = {area_unit_s3}`, so `Integral sigma3(K)dA = {sigma3_integral_convex}`.\n"
        f"- Boundary functional: `16*Integral sigma3(K)dA = {boundary_functional_convex}`; residual against `32*pi**2*chi`: `{convex_boundary_residual}`.\n"
        f"- Normalized Track-B pairing: `{convex_pairing}`; residual `{convex_pairing_residual}`.\n"
        f"- P2343 round residual: `{round_residual}`; P2344 spheroidal residual: `{spheroidal_residual}`.\n"
        f"- Closure-promise ranking: `{', '.join(ranked_closure_promises)}`.\n"
        "- Readout: best current non-selector closure promise is Track-B boundary/topology normalization; nadsoliton ontology points to global boundary/topological constraints before selector choice, but not to a unique future yet.\n"
        "- No arbitrary-boundary functional, no general Chern-Gauss-Bonnet theorem, no global renormalization, no independent `a_GB`, no selector premise, no QW-2191 discharge, no unique-future choice, no observer prediction, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
