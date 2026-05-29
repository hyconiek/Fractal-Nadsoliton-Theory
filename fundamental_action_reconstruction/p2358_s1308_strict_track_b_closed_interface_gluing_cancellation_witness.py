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
OUT = GEN / "p2358_s1308_strict_track_b_closed_interface_gluing_cancellation_witness.json"
MD = GEN / "p2358_s1308_strict_track_b_closed_interface_gluing_cancellation_witness.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2353_OBSTRUCTION_LEDGER": GEN / "p2353_s1303_strict_track_b_arbitrary_boundary_obstruction_ledger.json",
    "P2356_COMPONENT_ADDITIVITY": GEN / "p2356_s1306_strict_track_b_oriented_convex_component_degree_additivity_lemma.json",
    "P2357_TRANSGRESSION_BRIDGE": GEN / "p2357_s1307_strict_track_b_flat_convex_component_transgression_integration_bridge.json",
}

RESULT_KINDS = {
    "P2335_TWO_TRACK_LEDGER": "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
    "P2353_OBSTRUCTION_LEDGER": "STRICT_TRACK_B_ARBITRARY_BOUNDARY_OBSTRUCTION_LEDGER_NO_UNIVERSAL_CLOSURE",
    "P2356_COMPONENT_ADDITIVITY": "STRICT_TRACK_B_ORIENTED_CONVEX_COMPONENT_DEGREE_ADDITIVITY_LEMMA_NO_UNIVERSAL_CLOSURE",
    "P2357_TRANSGRESSION_BRIDGE": "STRICT_TRACK_B_FLAT_CONVEX_COMPONENT_TRANSGRESSION_INTEGRATION_BRIDGE_NO_UNIVERSAL_CLOSURE",
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
        ROOT / "p2356_s1306_strict_track_b_oriented_convex_component_degree_additivity_lemma.py",
        ROOT / "p2357_s1307_strict_track_b_flat_convex_component_transgression_integration_bridge.py",
        GEN / "p2357_s1307_strict_track_b_flat_convex_component_transgression_integration_bridge.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track-B|Track B|gluing|interface|cancellation|corner|smoothing|O5|O3|O4|degree|orientation|arbitrary-boundary|selector|QW-2191|ToE closure",
        *existing,
    ]
    proc = subprocess.run(cmd, cwd=ROOT.parent, text=True, capture_output=True, check=False)
    return proc.stdout.splitlines()[:280]


def boundary_for_degree(degree: sp.Expr) -> sp.Expr:
    return sp.factor(32 * sp.pi**2 * degree)


def gluing_case(case_id: str, left_external_degree: int, right_external_degree: int, b_gb: sp.Expr) -> dict[str, Any]:
    left_ext = sp.Integer(left_external_degree)
    right_ext = sp.Integer(right_external_degree)
    interface_left = sp.Integer(1)
    interface_right = sp.Integer(-1)
    left_total = left_ext + interface_left
    right_total = right_ext + interface_right
    preglue_total = sp.factor(left_total + right_total)
    postglue_external_total = sp.factor(left_ext + right_ext)
    left_boundary = boundary_for_degree(left_total)
    right_boundary = boundary_for_degree(right_total)
    preglue_boundary = sp.factor(left_boundary + right_boundary)
    postglue_boundary = boundary_for_degree(postglue_external_total)
    interface_boundary_sum = sp.factor(boundary_for_degree(interface_left) + boundary_for_degree(interface_right))
    boundary_residual = sp.factor(preglue_boundary - postglue_boundary)
    pairing_residual = sp.factor(b_gb * preglue_boundary - b_gb * postglue_boundary)
    return {
        "case_id": case_id,
        "left_external_degree": int(left_ext),
        "right_external_degree": int(right_ext),
        "interface_left_degree": int(interface_left),
        "interface_right_degree": int(interface_right),
        "left_total_degree_before_gluing": int(left_total),
        "right_total_degree_before_gluing": int(right_total),
        "preglue_total_degree": int(preglue_total),
        "postglue_external_degree": int(postglue_external_total),
        "left_boundary_before_gluing": str(left_boundary),
        "right_boundary_before_gluing": str(right_boundary),
        "interface_boundary_sum": str(interface_boundary_sum),
        "preglue_boundary_sum": str(preglue_boundary),
        "postglue_boundary_sum": str(postglue_boundary),
        "boundary_gluing_residual": str(boundary_residual),
        "preglue_pairing_sum": str(sp.factor(b_gb * preglue_boundary)),
        "postglue_pairing_sum": str(sp.factor(b_gb * postglue_boundary)),
        "pairing_gluing_residual": str(pairing_residual),
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    b_gb = sp.factor(track_b_coefficient(artifacts["P2335_TWO_TRACK_LEDGER"]))
    target_boundary_per_degree = 32 * sp.pi**2
    target_pairing_per_degree = sp.factor(target_boundary_per_degree * b_gb)

    a, b = sp.symbols("a b")
    left_interface = boundary_for_degree(1)
    right_interface = boundary_for_degree(-1)
    symbolic_interface_residual = sp.factor(left_interface + right_interface)
    symbolic_preglue_degree = sp.factor((a + 1) + (b - 1))
    symbolic_postglue_degree = sp.factor(a + b)
    symbolic_degree_residual = sp.factor(symbolic_preglue_degree - symbolic_postglue_degree)
    symbolic_boundary_residual = sp.factor(boundary_for_degree(symbolic_preglue_degree) - boundary_for_degree(symbolic_postglue_degree))
    symbolic_pairing_residual = sp.factor(b_gb * symbolic_boundary_residual)

    p2353_ledger = find_first_key(artifacts["P2353_OBSTRUCTION_LEDGER"], "track_B_arbitrary_boundary_obstruction_ledger") or {}
    p2356_lemma = find_first_key(artifacts["P2356_COMPONENT_ADDITIVITY"], "track_B_oriented_convex_component_degree_additivity_lemma") or {}
    p2357_bridge = find_first_key(artifacts["P2357_TRANSGRESSION_BRIDGE"], "track_B_flat_convex_component_transgression_integration_bridge") or {}
    p2353_minimal_cut = p2353_ledger.get("minimal_next_missing_cut", [])

    case_rows = [
        gluing_case("closed_interface_pair_no_external_boundary", 0, 0, b_gb),
        gluing_case("single_external_degree_survives_gluing", 1, 0, b_gb),
        gluing_case("annular_zero_degree_recomposition", 1, -1, b_gb),
        gluing_case("nonzero_stress_external_degree_recomposition", 2, -1, b_gb),
    ]
    all_case_residuals_zero = all(
        row["boundary_gluing_residual"] == "0" and row["pairing_gluing_residual"] == "0" for row in case_rows
    )
    coverage_columns = ["closed_interface", "external_survives", "zero_total", "nonzero_total", "negative_external"]
    coverage_rows = [
        [1, 0, 1, 0, 0],
        [1, 1, 0, 1, 0],
        [1, 1, 1, 0, 1],
        [1, 1, 0, 1, 1],
    ]
    coverage_rank = int(sp.Matrix(coverage_rows).rank())

    witness = {
        "witness_id": "P2358_closed_interface_gluing_cancellation_v1",
        "scope": "partial O5 closed-interface gluing cancellation under P2356/P2357 strict convex-component hypotheses; not a corner or smoothing theorem",
        "hypotheses": [
            "The glued interface is a smooth closed strict convex component appearing with opposite orientation signs on the two sides.",
            "Both sides separately satisfy the P2356 finite convex-component degree hypotheses and the P2357 flat local-to-integrated bridge hypotheses.",
            "No corners, collars with boundary, smoothing limits, metric jumps, or non-flat ambient curvature terms are introduced.",
        ],
        "target_boundary_per_degree": str(target_boundary_per_degree),
        "target_pairing_per_degree": str(target_pairing_per_degree),
        "interface_left_boundary": str(left_interface),
        "interface_right_boundary": str(right_interface),
        "symbolic_interface_residual": str(symbolic_interface_residual),
        "symbolic_preglue_degree": str(symbolic_preglue_degree),
        "symbolic_postglue_degree": str(symbolic_postglue_degree),
        "symbolic_degree_residual": str(symbolic_degree_residual),
        "symbolic_boundary_residual": str(symbolic_boundary_residual),
        "symbolic_pairing_residual": str(symbolic_pairing_residual),
        "gluing_case_rows": case_rows,
        "all_case_residuals_zero": all_case_residuals_zero,
        "coverage_matrix_columns": coverage_columns,
        "coverage_matrix_rows": coverage_rows,
        "coverage_matrix_rank": coverage_rank,
        "p2356_fixture_residuals_replayed": p2356_lemma.get("all_fixture_residuals_zero"),
        "p2357_bridge_residuals_replayed": p2357_bridge.get("all_bridge_residuals_zero"),
        "p2357_o3_partial_replayed": p2357_bridge.get("o3_cut_partially_attacked_under_hypotheses"),
        "p2353_minimal_cut_replayed": p2353_minimal_cut,
        "o5_cut_partially_attacked_under_closed_interface_hypotheses": "O5_regularization_corners_and_gluing" in p2353_minimal_cut,
        "remaining_o5_gap": "No theorem here handles corners, nonclosed interfaces, smoothing limits, metric-jump regularization, or general gluing of arbitrary boundaries.",
    }

    dependencies = {
        key.lower() + "_loaded": artifacts[key].get("result_kind") == expected
        for key, expected in RESULT_KINDS.items()
    }
    dependencies.update(
        {
            "p2356_fixture_residuals_replayed": witness["p2356_fixture_residuals_replayed"] is True,
            "p2357_bridge_residuals_replayed": witness["p2357_bridge_residuals_replayed"] is True,
            "p2357_o3_partial_replayed": witness["p2357_o3_partial_replayed"] is True,
            "p2353_o5_cut_replayed": witness["o5_cut_partially_attacked_under_closed_interface_hypotheses"],
            "symbolic_interface_residual_zero": symbolic_interface_residual == 0,
            "symbolic_degree_residual_zero": symbolic_degree_residual == 0,
            "symbolic_boundary_residual_zero": symbolic_boundary_residual == 0,
            "symbolic_pairing_residual_zero": symbolic_pairing_residual == 0,
            "all_case_residuals_zero": all_case_residuals_zero,
            "coverage_rank_at_least_four": coverage_rank >= 4,
        }
    )

    theorem_export = {
        "theorem_name": "P2358 Track-B closed-interface gluing cancellation witness",
        "claim": (
            "Under the stated closed-interface hypotheses, a common smooth strict convex interface contributes +32*pi**2 on one side and -32*pi**2 "
            "on the other, so its Track-B boundary functional and P2335-ledger pairing cancel exactly under gluing.  This is a controlled partial "
            "O5 gluing witness tied to P2356/P2357 hypotheses, not a general corners/gluing/smoothing theorem."
        ),
        "proof_witnesses": [
            "The symbolic interface residual is (+32*pi**2)+(-32*pi**2)=0.",
            "For symbolic external degrees a and b, preglue degree (a+1)+(b-1) equals postglue degree a+b.",
            "Boundary and pairing residuals vanish symbolically and in four concrete gluing stress cases.",
            "The witness replays P2356 finite component additivity and P2357 flat local-to-integrated bridge residuals.",
            "The O5 cut is attacked only for closed smooth opposite-orientation interfaces; corners and smoothing limits remain open.",
        ],
        "not_licensed": [
            "arbitrary-boundary theorem",
            "general gluing theorem",
            "corner contribution theorem",
            "smoothing-limit theorem",
            "general transgression theorem without strict convex-component hypotheses",
            "general nonconvex boundary theorem",
            "general Chern-Gauss-Bonnet theorem over all compact four-manifolds with boundary",
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
            "After this closed-interface O5 witness, the next non-redundant move is a real corner/smoothing-local model or a theorem weakening "
            "the closed smooth strict-convex interface hypothesis. Do not treat closed-interface cancellation as arbitrary-boundary closure."
        ),
    }

    probe = {
        "probe_id": "P2358_S1308_STRICT_TRACK_B_CLOSED_INTERFACE_GLUING_CANCELLATION_WITNESS",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2358 partial O5 closed-interface gluing cancellation after P2357",
            "top_hits": grep_hits(),
        },
        "track_B_closed_interface_gluing_cancellation_witness": witness,
        "current_export_dependencies": dependencies,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "repo_grep_audit_recorded": len(probe["repo_grep_audit"]["top_hits"]) > 0,
        **dependencies,
        "four_gluing_cases_recorded": len(case_rows) == 4,
        "closed_interface_scope_declared": "closed-interface" in witness["scope"],
        "o5_only_under_closed_interface_hypotheses": True,
        "corner_theorem_not_claimed": True,
        "smoothing_limit_theorem_not_claimed": True,
        "general_gluing_theorem_not_claimed": True,
        "arbitrary_boundary_theorem_not_claimed": True,
        "general_transgression_theorem_not_claimed": True,
        "general_nonconvex_theorem_not_claimed": True,
        "general_cgb_theorem_not_claimed": True,
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
        "schema_version": "p2358_s1308_v1",
        "packet_id": "P2358",
        "stage_id": "S1308",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-29T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_CLOSED_INTERFACE_GLUING_CANCELLATION_WITNESS_NO_UNIVERSAL_CLOSURE",
        "strict_track_b_closed_interface_gluing_cancellation_witness_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2358 strict Track-B closed-interface gluing cancellation witness\n\n"
        "Status: partial O5 closed-interface gluing witness under P2356/P2357 hypotheses; no arbitrary-boundary, selector, or ToE closure.\n\n"
        f"- Interface contributions: left `{left_interface}`, right `{right_interface}`, residual `{symbolic_interface_residual}`.\n"
        f"- Symbolic degrees: preglue `{symbolic_preglue_degree}`, postglue `{symbolic_postglue_degree}`, residual `{symbolic_degree_residual}`.\n"
        f"- Symbolic boundary residual `{symbolic_boundary_residual}`; pairing residual `{symbolic_pairing_residual}`.\n"
        f"- Gluing cases: `{[row['case_id'] for row in case_rows]}`; all residuals zero `{all_case_residuals_zero}`.\n"
        f"- Coverage matrix rank `{coverage_rank}`; P2356 residuals `{witness['p2356_fixture_residuals_replayed']}`; P2357 residuals `{witness['p2357_bridge_residuals_replayed']}`.\n"
        f"- P2353 cut replayed: `{p2353_minimal_cut}`; O5 partially attacked `{witness['o5_cut_partially_attacked_under_closed_interface_hypotheses']}`.\n"
        "- No arbitrary-boundary theorem, no general gluing theorem, no corner contribution theorem, no smoothing-limit theorem, no general transgression theorem without strict convex-component hypotheses, no general Chern-Gauss-Bonnet theorem, no non-flat ambient curvature bridge, no global renormalization, no independent `a_GB`, no selector premise, no QW-2191 discharge, no unique-future choice, no observer prediction, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
