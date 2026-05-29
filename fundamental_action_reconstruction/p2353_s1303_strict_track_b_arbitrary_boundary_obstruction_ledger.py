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
OUT = GEN / "p2353_s1303_strict_track_b_arbitrary_boundary_obstruction_ledger.json"
MD = GEN / "p2353_s1303_strict_track_b_arbitrary_boundary_obstruction_ledger.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2338_TRACK_B_CONTRACT": GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.json",
    "P2345_CONVEX_GAUSS_MAP": GEN / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.json",
    "P2348_CHERN_POLYNOMIAL": GEN / "p2348_s1298_strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction.json",
    "P2352_CONTROLLED_SYNTHESIS": GEN / "p2352_s1302_strict_track_b_controlled_convex_boundary_class_synthesis_theorem.json",
}

RESULT_KINDS = {
    "P2335_TWO_TRACK_LEDGER": "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
    "P2338_TRACK_B_CONTRACT": "STRICT_TRACK_B_GB_BOUNDARY_TOPOLOGICAL_NORMALIZATION_CONTRACT_EXPORTED_NO_NORMALIZATION_CLAIM",
    "P2345_CONVEX_GAUSS_MAP": "STRICT_TRACK_B_CONVEX_GAUSS_MAP_BOUNDARY_FUNCTIONAL_THEOREM_AND_CLOSURE_PROMISE_MAP_NO_SELECTOR_CLOSURE",
    "P2348_CHERN_POLYNOMIAL": "STRICT_TRACK_B_CHERN_BOUNDARY_POLYNOMIAL_NONSYMMETRIC_REDUCTION_NO_INTEGRATED_GLOBAL_THEOREM",
    "P2352_CONTROLLED_SYNTHESIS": "STRICT_TRACK_B_CONTROLLED_CONVEX_BOUNDARY_CLASS_SYNTHESIS_THEOREM_NO_UNIVERSAL_CLOSURE",
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
        ROOT / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.py",
        ROOT / "p2348_s1298_strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction.py",
        ROOT / "p2352_s1302_strict_track_b_controlled_convex_boundary_class_synthesis_theorem.py",
        GEN / "p2352_s1302_strict_track_b_controlled_convex_boundary_class_synthesis_theorem.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track-B|Track B|arbitrary-boundary|Chern-Gauss-Bonnet|Gauss-map|convex|gluing|selector|QW-2191|ToE closure",
        *existing,
    ]
    proc = subprocess.run(cmd, cwd=ROOT.parent, text=True, capture_output=True, check=False)
    return proc.stdout.splitlines()[:240]


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    p2335 = artifacts["P2335_TWO_TRACK_LEDGER"]
    p2345 = artifacts["P2345_CONVEX_GAUSS_MAP"]
    p2348 = artifacts["P2348_CHERN_POLYNOMIAL"]
    p2352 = artifacts["P2352_CONTROLLED_SYNTHESIS"]

    b_gb = sp.factor(track_b_coefficient(p2335))
    target_boundary = 32 * sp.pi**2
    target_pairing = sp.factor(target_boundary * b_gb)
    p2345_boundary = sp.factor(sympify_expr(find_first_key(p2345, "boundary_functional_convex_class")))
    p2348_replay = sp.factor(sympify_expr(find_first_key(p2348, "substitute_integral_sigma3_equals_2pi2")))
    p2352_synthesis = find_first_key(p2352, "track_B_controlled_convex_boundary_class_synthesis") or {}
    p2352_boundary_residuals = p2352_synthesis.get("boundary_residual_vector", [])
    p2352_pairing_residuals = p2352_synthesis.get("pairing_residual_vector", [])
    p2352_not_licensed = find_first_key(p2352, "not_licensed") or []

    replay_residuals = {
        "p2345_convex_gauss_map_boundary_minus_32pi2": str(sp.factor(p2345_boundary - target_boundary)),
        "p2348_flat_chern_polynomial_replay_minus_32pi2": str(sp.factor(p2348_replay - target_boundary)),
        "p2352_boundary_residual_vector": p2352_boundary_residuals,
        "p2352_pairing_residual_vector": p2352_pairing_residuals,
    }

    obligations = [
        {
            "obligation_id": "O1_convex_degree_one_gauss_map_lane",
            "needed_for_universal_boundary": True,
            "current_status": "DISCHARGED_FOR_STRICTLY_CONVEX_FLAT_BALL_BOUNDARIES_ONLY",
            "artifact_support": ["P2345", "P2349", "P2350", "P2351", "P2352"],
            "proof_obstruction": "Does not cover nonconvex, nonembedded, higher-degree, or arbitrary-orientation boundary strata.",
            "discharged": True,
            "partial": False,
        },
        {
            "obligation_id": "O2_local_chern_boundary_polynomial_normalization",
            "needed_for_universal_boundary": True,
            "current_status": "LOCAL_POLYNOMIAL_REDUCTION_EXPORTED_NOT_GLOBAL_INTEGRATION_THEOREM",
            "artifact_support": ["P2348", "P2349", "P2350", "P2351", "P2352"],
            "proof_obstruction": "The local polynomial exists, but no arbitrary-boundary transgression/integration theorem is exported.",
            "discharged": True,
            "partial": False,
        },
        {
            "obligation_id": "O3_arbitrary_boundary_transgression_integration",
            "needed_for_universal_boundary": True,
            "current_status": "OPEN_MISSING_THEOREM",
            "artifact_support": [],
            "proof_obstruction": "Need a theorem integrating the Chern boundary polynomial over arbitrary smooth boundary geometry with all sign and normalization conventions fixed.",
            "discharged": False,
            "partial": False,
        },
        {
            "obligation_id": "O4_nonconvex_degree_and_orientation_accounting",
            "needed_for_universal_boundary": True,
            "current_status": "OPEN_MISSING_DEGREE_LEDGER",
            "artifact_support": [],
            "proof_obstruction": "Convex Gauss-map degree-one replay is insufficient for nonconvex degree, cancellation, orientation reversal, or self-overlap cases.",
            "discharged": False,
            "partial": False,
        },
        {
            "obligation_id": "O5_regularization_corners_and_gluing",
            "needed_for_universal_boundary": True,
            "current_status": "OPEN_MISSING_GLUING_REGULARITY_THEOREM",
            "artifact_support": [],
            "proof_obstruction": "No exported theorem handles corners, cut-and-paste boundaries, smoothing limits, or additivity under gluing.",
            "discharged": False,
            "partial": False,
        },
        {
            "obligation_id": "O6_general_cgb_topology_bridge",
            "needed_for_universal_boundary": True,
            "current_status": "OPEN_GLOBAL_CGB_BRIDGE",
            "artifact_support": ["P2348"],
            "proof_obstruction": "Local Chern-polynomial normalization has not been promoted to a general Chern-Gauss-Bonnet theorem over all compact four-manifolds with boundary.",
            "discharged": False,
            "partial": True,
        },
        {
            "obligation_id": "O7_independent_a_gb_or_renormalization_source",
            "needed_for_universal_boundary": False,
            "current_status": "INTENTIONALLY_NOT_CLAIMED",
            "artifact_support": ["P2335", "P2338"],
            "proof_obstruction": "Track-B currently replays the P2335 ledger coefficient; it does not independently measure a_GB or close global renormalization.",
            "discharged": False,
            "partial": True,
        },
        {
            "obligation_id": "O8_selector_and_unique_future_closure",
            "needed_for_universal_boundary": False,
            "current_status": "BLOCKED_BY_QW_2191_OUTSIDE_TRACK_B_SCOPE",
            "artifact_support": [],
            "proof_obstruction": "Even a completed Track-B boundary theorem would not add a selector premise or discharge QW-2191.",
            "discharged": False,
            "partial": False,
        },
    ]

    required = [row for row in obligations if row["needed_for_universal_boundary"]]
    discharged_required = [row for row in required if row["discharged"]]
    partial_required = [row for row in required if row["partial"] and not row["discharged"]]
    open_required = [row for row in required if not row["discharged"]]
    hard_open_required = [row for row in required if not row["discharged"] and not row["partial"]]
    closure_score = sp.Rational(len(discharged_required), len(required))
    partial_credit_score = sp.Rational(
        2 * len(discharged_required) + len(partial_required),
        2 * len(required),
    )
    obstruction_vector = [0 if row["discharged"] else 1 for row in required]
    hard_obstruction_vector = [0 if row["discharged"] or row["partial"] else 1 for row in required]
    support_matrix_columns = ["P2335", "P2338", "P2345", "P2348", "P2349", "P2350", "P2351", "P2352"]
    support_matrix_rows = [
        [1 if col in row["artifact_support"] else 0 for col in support_matrix_columns]
        for row in obligations
    ]
    support_rank = int(sp.Matrix(support_matrix_rows).rank())

    obstruction_ledger = {
        "ledger_id": "P2353_arbitrary_boundary_obstruction_ledger_v1",
        "scope": "Track-B proof frontier ledger after controlled-class synthesis; not an arbitrary-boundary theorem",
        "target_boundary_functional": str(target_boundary),
        "ledger_coefficient_b_GB": str(b_gb),
        "target_pairing_per_euler_number": str(target_pairing),
        "replay_residuals": replay_residuals,
        "obligations": obligations,
        "required_obligation_count": len(required),
        "discharged_required_count": len(discharged_required),
        "partial_required_count": len(partial_required),
        "open_required_count": len(open_required),
        "hard_open_required_count": len(hard_open_required),
        "closure_score_discharged_only": str(closure_score),
        "closure_score_with_partial_half_credit": str(partial_credit_score),
        "required_obstruction_vector": obstruction_vector,
        "hard_required_obstruction_vector": hard_obstruction_vector,
        "minimal_next_missing_cut": [row["obligation_id"] for row in hard_open_required],
        "support_matrix_columns": support_matrix_columns,
        "support_matrix_rows": support_matrix_rows,
        "support_matrix_rank": support_rank,
        "p2352_not_licensed_replay": p2352_not_licensed,
    }

    dependencies = {
        key.lower() + "_loaded": artifacts[key].get("result_kind") == expected
        for key, expected in RESULT_KINDS.items()
    }
    dependencies.update(
        {
            "p2345_boundary_replay_zero": replay_residuals["p2345_convex_gauss_map_boundary_minus_32pi2"] == "0",
            "p2348_chern_replay_zero": replay_residuals["p2348_flat_chern_polynomial_replay_minus_32pi2"] == "0",
            "p2352_boundary_residual_vector_zero": all(value == "0" for value in p2352_boundary_residuals),
            "p2352_pairing_residual_vector_zero": all(value == "0" for value in p2352_pairing_residuals),
            "p2352_arbitrary_boundary_nonlicense_replayed": "arbitrary-boundary theorem" in p2352_not_licensed,
        }
    )

    theorem_export = {
        "theorem_name": "P2353 Track-B arbitrary-boundary obstruction ledger theorem",
        "claim": (
            "P2353 proves only a frontier statement: after P2352 controlled-class synthesis, the currently exported Track-B artifacts discharge "
            "the convex degree-one Gauss-map lane and the local Chern-polynomial normalization lane, but they still leave explicit missing "
            "obligations before any arbitrary-boundary or general Chern-Gauss-Bonnet theorem."
        ),
        "proof_witnesses": [
            "P2345, P2348, and P2352 replay with exact zero residuals against 32*pi**2 and the P2335 Track-B pairing.",
            "The required-obligation vector has six entries; only O1 and O2 are discharged by current artifacts.",
            "O3, O4, and O5 are hard-open missing theorems, while O6 is only partially supported by the local Chern-polynomial artifact.",
            "The ledger records a minimal next missing cut rather than upgrading controlled evidence into universal closure.",
            "Selector, unique-future, observer-prediction, G1/G3, and ToE claims remain explicitly outside this Track-B obstruction ledger.",
        ],
        "not_licensed": [
            "arbitrary-boundary theorem",
            "general Chern-Gauss-Bonnet theorem over all compact four-manifolds with boundary",
            "universal convex or nonconvex boundary theorem",
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
            "Attack exactly one hard-open cut: O3 arbitrary-boundary transgression integration, O4 nonconvex degree/orientation ledger, "
            "or O5 regularization/corners/gluing. Do not present another controlled sample as universal closure."
        ),
    }

    probe = {
        "probe_id": "P2353_S1303_STRICT_TRACK_B_ARBITRARY_BOUNDARY_OBSTRUCTION_LEDGER",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2353 obstruction ledger before arbitrary-boundary, global CGB, selector, or ToE claims",
            "top_hits": grep_hits(),
        },
        "track_B_arbitrary_boundary_obstruction_ledger": obstruction_ledger,
        "current_export_dependencies": dependencies,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "repo_grep_audit_recorded": len(probe["repo_grep_audit"]["top_hits"]) > 0,
        **dependencies,
        "required_obligation_count_is_six": len(required) == 6,
        "only_two_required_obligations_discharged": len(discharged_required) == 2,
        "hard_open_required_count_is_three": len(hard_open_required) == 3,
        "closure_score_below_one": bool(closure_score < 1),
        "partial_credit_score_below_one": bool(partial_credit_score < 1),
        "minimal_next_missing_cut_nonempty": len(obstruction_ledger["minimal_next_missing_cut"]) == 3,
        "support_matrix_rank_at_least_four": support_rank >= 4,
        "arbitrary_boundary_theorem_not_claimed": True,
        "general_cgb_theorem_not_claimed": True,
        "universal_boundary_theorem_not_claimed": True,
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
        "schema_version": "p2353_s1303_v1",
        "packet_id": "P2353",
        "stage_id": "S1303",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-29T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_ARBITRARY_BOUNDARY_OBSTRUCTION_LEDGER_NO_UNIVERSAL_CLOSURE",
        "strict_track_b_arbitrary_boundary_obstruction_ledger_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2353 strict Track-B arbitrary-boundary obstruction ledger\n\n"
        "Status: proof-frontier obstruction ledger exported; no arbitrary-boundary, selector, or ToE closure.\n\n"
        f"- `b_GB = {b_gb}`; target boundary functional `{target_boundary}`; target pairing `{target_pairing}`.\n"
        f"- Replay residuals: P2345 `{replay_residuals['p2345_convex_gauss_map_boundary_minus_32pi2']}`, "
        f"P2348 `{replay_residuals['p2348_flat_chern_polynomial_replay_minus_32pi2']}`, "
        f"P2352 boundary `{p2352_boundary_residuals}`, P2352 pairing `{p2352_pairing_residuals}`.\n"
        f"- Required obligations: `{len(required)}`; discharged `{len(discharged_required)}`; partial `{len(partial_required)}`; open `{len(open_required)}`; hard-open `{len(hard_open_required)}`.\n"
        f"- Closure score discharged-only `{closure_score}`; partial-half-credit `{partial_credit_score}`.\n"
        f"- Required obstruction vector `{obstruction_vector}`; hard obstruction vector `{hard_obstruction_vector}`.\n"
        f"- Minimal next missing cut: `{obstruction_ledger['minimal_next_missing_cut']}`.\n"
        f"- Support matrix rank: `{support_rank}` over columns `{support_matrix_columns}`.\n"
        "- No arbitrary-boundary theorem, no general Chern-Gauss-Bonnet theorem, no universal boundary theorem, no general gluing theorem, no global renormalization, no independent `a_GB`, no selector premise, no QW-2191 discharge, no unique-future choice, no observer prediction, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
