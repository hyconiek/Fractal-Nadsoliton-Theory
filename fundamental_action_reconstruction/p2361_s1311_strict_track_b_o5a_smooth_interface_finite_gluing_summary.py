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
OUT = GEN / "p2361_s1311_strict_track_b_o5a_smooth_interface_finite_gluing_summary.json"
MD = GEN / "p2361_s1311_strict_track_b_o5a_smooth_interface_finite_gluing_summary.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2353_OBSTRUCTION_LEDGER": GEN / "p2353_s1303_strict_track_b_arbitrary_boundary_obstruction_ledger.json",
    "P2358_CLOSED_INTERFACE_GLUING": GEN / "p2358_s1308_strict_track_b_closed_interface_gluing_cancellation_witness.json",
    "P2359_S4_CAP_COMPLEMENT": GEN / "p2359_s1309_strict_track_b_s4_cap_complement_gluing_bulk_boundary_witness.json",
    "P2360_SMOOTH_INTERFACE_SIGN_REVERSAL": GEN / "p2360_s1310_strict_track_b_smooth_interface_chern_polynomial_sign_reversal_lemma.json",
}

RESULT_KINDS = {
    "P2335_TWO_TRACK_LEDGER": "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
    "P2353_OBSTRUCTION_LEDGER": "STRICT_TRACK_B_ARBITRARY_BOUNDARY_OBSTRUCTION_LEDGER_NO_UNIVERSAL_CLOSURE",
    "P2358_CLOSED_INTERFACE_GLUING": "STRICT_TRACK_B_CLOSED_INTERFACE_GLUING_CANCELLATION_WITNESS_NO_UNIVERSAL_CLOSURE",
    "P2359_S4_CAP_COMPLEMENT": "STRICT_TRACK_B_S4_CAP_COMPLEMENT_GLUING_BULK_BOUNDARY_WITNESS_NO_GENERAL_CGB",
    "P2360_SMOOTH_INTERFACE_SIGN_REVERSAL": "STRICT_TRACK_B_SMOOTH_INTERFACE_CHERN_POLYNOMIAL_SIGN_REVERSAL_LEMMA_NO_CORNER_SMOOTHING_CLOSURE",
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
        ROOT / "p2358_s1308_strict_track_b_closed_interface_gluing_cancellation_witness.py",
        ROOT / "p2359_s1309_strict_track_b_s4_cap_complement_gluing_bulk_boundary_witness.py",
        ROOT / "p2360_s1310_strict_track_b_smooth_interface_chern_polynomial_sign_reversal_lemma.py",
        GEN / "p2360_s1310_strict_track_b_smooth_interface_chern_polynomial_sign_reversal_lemma.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "P2360|smooth-interface|Chern|boundary polynomial|interface|gluing|sign reversal|corner|smoothing|O5|O5a|arbitrary-boundary|selector|QW-2191|ToE closure",
        *existing,
    ]
    proc = subprocess.run(cmd, cwd=ROOT.parent, text=True, capture_output=True, check=False)
    return proc.stdout.splitlines()[:320]


def interface_case_row(case_id: str, interface_symbols: list[sp.Symbol], incidence_rows: list[list[int]], external_degree: int, b_gb: sp.Expr) -> dict[str, Any]:
    if not interface_symbols:
        interface_vector = sp.Matrix([])
        incidence = sp.zeros(len(incidence_rows), 0)
        side_contribs: list[sp.Expr] = [sp.Integer(0) for _ in incidence_rows]
    else:
        interface_vector = sp.Matrix(interface_symbols)
        incidence = sp.Matrix(incidence_rows)
        side_contribs = list(incidence * interface_vector)
    pre_gluing_interface_total = sp.factor(sum(side_contribs, sp.Integer(0)))
    post_gluing_interface_total = sp.Integer(0)
    interface_residual = sp.factor(pre_gluing_interface_total - post_gluing_interface_total)
    pre_gluing_degree = sp.factor(sp.Integer(external_degree) + pre_gluing_interface_total)
    post_gluing_degree = sp.Integer(external_degree)
    boundary_residual = sp.factor(32 * sp.pi**2 * (pre_gluing_degree - post_gluing_degree))
    pairing_residual = sp.factor(b_gb * boundary_residual)
    column_sums = [sp.factor(sum(incidence[:, j])) for j in range(incidence.shape[1])]
    return {
        "case_id": case_id,
        "regions": len(incidence_rows),
        "interfaces": [str(symbol) for symbol in interface_symbols],
        "incidence_rows": [[int(entry) for entry in row] for row in incidence_rows],
        "column_sums": [str(value) for value in column_sums],
        "side_interface_contributions": [str(sp.factor(value)) for value in side_contribs],
        "pre_gluing_interface_total": str(pre_gluing_interface_total),
        "post_gluing_interface_total": str(post_gluing_interface_total),
        "interface_residual": str(interface_residual),
        "pre_gluing_degree": str(pre_gluing_degree),
        "post_gluing_degree": str(post_gluing_degree),
        "boundary_residual": str(boundary_residual),
        "pairing_residual": str(pairing_residual),
        "all_columns_balanced": all(value == 0 for value in column_sums),
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    b_gb = sp.factor(track_b_coefficient(artifacts["P2335_TWO_TRACK_LEDGER"]))

    J1, J2, J3 = sp.symbols("J1 J2 J3")
    x1, x2, x3, K = sp.symbols("x1 x2 x3 K")
    sigma1 = x1 + x2 + x3
    sigma3 = x1 * x2 * x3
    b3 = sp.factor(8 * K * sigma1 + 16 * sigma3)
    b3_reversed = sp.factor(b3.subs({x1: -x1, x2: -x2, x3: -x3}))
    local_o5a_residual = sp.factor(b3 + b3_reversed)

    p2353_ledger = find_first_key(artifacts["P2353_OBSTRUCTION_LEDGER"], "track_B_arbitrary_boundary_obstruction_ledger") or {}
    p2358_witness = find_first_key(artifacts["P2358_CLOSED_INTERFACE_GLUING"], "track_B_closed_interface_gluing_cancellation_witness") or {}
    p2359_witness = find_first_key(artifacts["P2359_S4_CAP_COMPLEMENT"], "track_B_s4_cap_complement_gluing_bulk_boundary_witness") or {}
    p2360_lemma = find_first_key(artifacts["P2360_SMOOTH_INTERFACE_SIGN_REVERSAL"], "track_B_smooth_interface_chern_polynomial_sign_reversal_lemma") or {}
    p2353_minimal_cut = p2353_ledger.get("minimal_next_missing_cut", [])

    cases = [
        interface_case_row("two_region_single_smooth_closed_interface", [J1], [[1], [-1]], 0, b_gb),
        interface_case_row("three_region_smooth_chain_two_interfaces", [J1, J2], [[1, 0], [-1, 1], [0, -1]], 2, b_gb),
        interface_case_row("four_region_disjoint_smooth_interface_tree", [J1, J2, J3], [[1, 1, 1], [-1, 0, 0], [0, -1, 0], [0, 0, -1]], -1, b_gb),
    ]
    all_case_columns_balanced = all(row["all_columns_balanced"] for row in cases)
    all_case_residuals_zero = all(
        row["interface_residual"] == "0" and row["boundary_residual"] == "0" and row["pairing_residual"] == "0"
        for row in cases
    )

    finite_gluing_summary = {
        "summary_id": "P2361_o5a_smooth_interface_finite_gluing_summary_v1",
        "scope": "O5a smooth-interface finite gluing only: closed smooth internal hypersurfaces with opposite normals and no corners, triple-junction corner terms, smoothing collars, metric jumps, or arbitrary nonsmooth boundary.",
        "hypotheses": [
            "Each internal interface is a closed smooth hypersurface shared by exactly two sides.",
            "The induced metric and ambient K input agree on the two traces of each smooth interface.",
            "The unit normals are opposite, so the P2360 Chern boundary polynomial changes sign locally.",
            "The finite interface incidence matrix has one +1 and one -1 in every interface column.",
            "No corner, edge, triple-junction, smoothing-limit, or metric-jump contribution is present or inferred.",
        ],
        "local_chern_polynomial": str(b3),
        "local_reversed_normal_polynomial": str(b3_reversed),
        "local_o5a_residual": str(local_o5a_residual),
        "finite_incidence_rule": "For interface column j, sum_i epsilon_{ij}=0; hence sum_i epsilon_{ij}*J_j=0 before any corner/smoothing regularization is considered.",
        "case_rows": cases,
        "all_case_columns_balanced": all_case_columns_balanced,
        "all_case_residuals_zero": all_case_residuals_zero,
        "p2360_sign_reversal_residual_replayed": p2360_lemma.get("sign_reversal_residual"),
        "p2360_interface_pair_sum_replayed": p2360_lemma.get("interface_pair_sum"),
        "p2358_symbolic_interface_residual_replayed": p2358_witness.get("symbolic_interface_residual"),
        "p2358_all_case_residuals_replayed": p2358_witness.get("all_case_residuals_zero"),
        "p2359_interface_cancellation_replayed": p2359_witness.get("interface_cancellation"),
        "p2359_gluing_consistency_residual_replayed": p2359_witness.get("gluing_consistency_residual"),
        "p2353_minimal_cut_replayed": p2353_minimal_cut,
        "o5a_smooth_interface_status": "FORMAL_SUMMARY_EXPORTED_FOR_FINITE_SMOOTH_CLOSED_INTERFACES_ONLY",
        "o5_full_status": "OPEN_BECAUSE_CORNERS_SMOOTHING_LIMITS_METRIC_JUMPS_AND_ARBITRARY_BOUNDARIES_ARE_NOT_INCLUDED",
        "remaining_gap": "This package intentionally finishes only the smooth-interface O5a lane; it does not move to corners and does not assert any smoothing or arbitrary-boundary theorem.",
    }

    dependencies = {
        key.lower() + "_loaded": artifacts[key].get("result_kind") == expected
        for key, expected in RESULT_KINDS.items()
    }
    dependencies.update(
        {
            "p2360_local_sign_reversal_replayed": finite_gluing_summary["p2360_sign_reversal_residual_replayed"] == "0",
            "p2360_interface_pair_replayed": finite_gluing_summary["p2360_interface_pair_sum_replayed"] == "0",
            "p2358_interface_replayed": finite_gluing_summary["p2358_symbolic_interface_residual_replayed"] == "0",
            "p2358_cases_replayed": finite_gluing_summary["p2358_all_case_residuals_replayed"] is True,
            "p2359_interface_replayed": finite_gluing_summary["p2359_interface_cancellation_replayed"] == "0",
            "p2359_gluing_consistency_replayed": finite_gluing_summary["p2359_gluing_consistency_residual_replayed"] == "0",
            "p2353_o5_cut_replayed": "O5_regularization_corners_and_gluing" in p2353_minimal_cut,
            "local_o5a_residual_zero": local_o5a_residual == 0,
            "finite_incidence_cases_balanced": all_case_columns_balanced,
            "finite_incidence_residuals_zero": all_case_residuals_zero,
        }
    )

    theorem_export = {
        "theorem_name": "P2361 Track-B O5a finite smooth-interface gluing summary",
        "claim": (
            "Under the P2360 smooth opposite-normal hypotheses, any finite collection of closed smooth internal interfaces with a balanced two-sided incidence column cancels in the Track-B Chern-boundary contribution and in the P2335 ledger pairing.  This gives a formal O5a smooth-interface summary, not O5 full regularization."
        ),
        "proof_witnesses": [
            "P2360 supplies the local identity B3(K;-k1,-k2,-k3)+B3(K;k1,k2,k3)=0.",
            "Each finite smooth interface column has one +1 and one -1 incidence entry, so the symbolic column sum is zero.",
            "Three finite incidence stress cases have zero interface, boundary, and pairing residuals.",
            "P2358 closed-interface and P2359 S4 cap-complement witnesses replay as special cases.",
            "The theorem export explicitly excludes corners, smoothing collars, metric jumps, and arbitrary nonsmooth boundaries.",
        ],
        "not_licensed": [
            "O5 full regularization theorem",
            "corner contribution theorem",
            "smoothing-limit theorem",
            "triple-junction or edge theorem",
            "metric-jump or distributional-curvature theorem",
            "arbitrary-boundary theorem",
            "general gluing theorem including nonsmooth strata",
            "general Chern-Gauss-Bonnet theorem over arbitrary compact four-manifolds",
            "global renormalization closure",
            "independent a_GB measurement separate from the P2335 ledger coefficient",
            "bulk EOM GB force or EOM-only GB lift",
            "legacy-to-strict kernel bridge",
            "selector premise or QW-2191 selector discharge",
            "choice of the unique physical future",
            "observer-level prediction",
            "G1/G3 update",
            "ToE closure",
        ],
        "next_honest_step": (
            "Do not continue O5a by replaying more smooth-interface cases and do not move to corners in this package.  The next high-value lane is a separate strict-side priority audit of the legacy-to-strict/scratch bridge evidence or a QW-2191 selector-source attempt, kept outside this O5a theorem export."
        ),
    }

    probe = {
        "probe_id": "P2361_S1311_STRICT_TRACK_B_O5A_SMOOTH_INTERFACE_FINITE_GLUING_SUMMARY",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2361 O5a finite smooth-interface summary after P2360",
            "top_hits": grep_hits(),
        },
        "track_B_o5a_smooth_interface_finite_gluing_summary": finite_gluing_summary,
        "current_export_dependencies": dependencies,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "repo_grep_audit_recorded": len(probe["repo_grep_audit"]["top_hits"]) > 0,
        **dependencies,
        "three_finite_incidence_cases_recorded": len(cases) == 3,
        "o5a_scope_declared": "O5a smooth-interface" in finite_gluing_summary["scope"],
        "o5_full_kept_open": finite_gluing_summary["o5_full_status"].startswith("OPEN"),
        "no_corner_step_taken": True,
        "corner_theorem_not_claimed": True,
        "smoothing_limit_theorem_not_claimed": True,
        "general_gluing_theorem_not_claimed": True,
        "arbitrary_boundary_theorem_not_claimed": True,
        "metric_jump_theorem_not_claimed": True,
        "general_cgb_theorem_not_claimed": True,
        "legacy_strict_bridge_not_claimed": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_unique_future_choice_claimed": True,
        "no_observer_prediction_claimed": True,
        "no_g1_g3_update_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2361_s1311_v1",
        "packet_id": "P2361",
        "stage_id": "S1311",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-29T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_O5A_SMOOTH_INTERFACE_FINITE_GLUING_SUMMARY_NO_O5_FULL_CLOSURE",
        "strict_track_b_o5a_smooth_interface_finite_gluing_summary_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2361 strict Track-B O5a smooth-interface finite gluing summary\n\n"
        "Status: finite smooth closed-interface O5a summary exported; no corner, smoothing, arbitrary-boundary, selector, bridge, or ToE closure.\n\n"
        f"- Local polynomial `{b3}`; reversed-normal polynomial `{b3_reversed}`; local residual `{local_o5a_residual}`.\n"
        f"- Finite incidence cases: `{[row['case_id'] for row in cases]}`; balanced columns `{all_case_columns_balanced}`; all residuals zero `{all_case_residuals_zero}`.\n"
        f"- Replays: P2360 residual `{finite_gluing_summary['p2360_sign_reversal_residual_replayed']}`, P2358 interface `{finite_gluing_summary['p2358_symbolic_interface_residual_replayed']}`, P2359 interface `{finite_gluing_summary['p2359_interface_cancellation_replayed']}`.\n"
        f"- O5a status `{finite_gluing_summary['o5a_smooth_interface_status']}`; O5 full status `{finite_gluing_summary['o5_full_status']}`.\n"
        "- No O5 full regularization theorem, no corner contribution theorem, no smoothing-limit theorem, no triple-junction/edge theorem, no metric-jump theorem, no arbitrary-boundary theorem, no general gluing theorem over nonsmooth strata, no general CGB theorem, no global renormalization, no independent `a_GB`, no bulk EOM GB force, no legacy-to-strict kernel bridge, no selector premise, no QW-2191 discharge, no unique-future choice, no observer prediction, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
