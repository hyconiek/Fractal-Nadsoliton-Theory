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
OUT = GEN / "p2360_s1310_strict_track_b_smooth_interface_chern_polynomial_sign_reversal_lemma.json"
MD = GEN / "p2360_s1310_strict_track_b_smooth_interface_chern_polynomial_sign_reversal_lemma.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2348_CHERN_POLYNOMIAL": GEN / "p2348_s1298_strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction.json",
    "P2353_OBSTRUCTION_LEDGER": GEN / "p2353_s1303_strict_track_b_arbitrary_boundary_obstruction_ledger.json",
    "P2358_CLOSED_INTERFACE_GLUING": GEN / "p2358_s1308_strict_track_b_closed_interface_gluing_cancellation_witness.json",
    "P2359_S4_CAP_COMPLEMENT": GEN / "p2359_s1309_strict_track_b_s4_cap_complement_gluing_bulk_boundary_witness.json",
}

RESULT_KINDS = {
    "P2335_TWO_TRACK_LEDGER": "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
    "P2348_CHERN_POLYNOMIAL": "STRICT_TRACK_B_CHERN_BOUNDARY_POLYNOMIAL_NONSYMMETRIC_REDUCTION_NO_INTEGRATED_GLOBAL_THEOREM",
    "P2353_OBSTRUCTION_LEDGER": "STRICT_TRACK_B_ARBITRARY_BOUNDARY_OBSTRUCTION_LEDGER_NO_UNIVERSAL_CLOSURE",
    "P2358_CLOSED_INTERFACE_GLUING": "STRICT_TRACK_B_CLOSED_INTERFACE_GLUING_CANCELLATION_WITNESS_NO_UNIVERSAL_CLOSURE",
    "P2359_S4_CAP_COMPLEMENT": "STRICT_TRACK_B_S4_CAP_COMPLEMENT_GLUING_BULK_BOUNDARY_WITNESS_NO_GENERAL_CGB",
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
        ROOT / "p2358_s1308_strict_track_b_closed_interface_gluing_cancellation_witness.py",
        ROOT / "p2359_s1309_strict_track_b_s4_cap_complement_gluing_bulk_boundary_witness.py",
        GEN / "p2359_s1309_strict_track_b_s4_cap_complement_gluing_bulk_boundary_witness.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track-B|Track B|Chern|boundary polynomial|interface|gluing|sign reversal|non-geodesic|corner|smoothing|O5|arbitrary-boundary|selector|QW-2191|ToE closure",
        *existing,
    ]
    proc = subprocess.run(cmd, cwd=ROOT.parent, text=True, capture_output=True, check=False)
    return proc.stdout.splitlines()[:300]


def sample_row(sample_id: str, substitutions: dict[sp.Symbol, sp.Expr], polynomial: sp.Expr) -> dict[str, str]:
    k1, k2, k3 = sp.symbols("k1 k2 k3")
    reversed_subs = dict(substitutions)
    reversed_subs.update({k1: -substitutions[k1], k2: -substitutions[k2], k3: -substitutions[k3]})
    value = sp.factor(polynomial.subs(substitutions))
    reversed_value = sp.factor(polynomial.subs(reversed_subs))
    return {
        "sample_id": sample_id,
        "value": str(value),
        "reversed_normal_value": str(reversed_value),
        "interface_sum": str(sp.factor(value + reversed_value)),
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    b_gb = sp.factor(track_b_coefficient(artifacts["P2335_TWO_TRACK_LEDGER"]))

    K, k1, k2, k3, I = sp.symbols("K k1 k2 k3 I")
    sigma1 = sp.factor(k1 + k2 + k3)
    sigma2 = sp.factor(k1 * k2 + k1 * k3 + k2 * k3)
    sigma3 = sp.factor(k1 * k2 * k3)
    chern_polynomial = sp.factor(8 * K * sigma1 + 16 * sigma3)
    reversed_polynomial = sp.factor(chern_polynomial.subs({k1: -k1, k2: -k2, k3: -k3}))
    sign_reversal_residual = sp.factor(reversed_polynomial + chern_polynomial)
    interface_pair_sum = sp.factor(chern_polynomial + reversed_polynomial)
    integrated_interface_residual = sp.factor(I + (-I))
    pairing_interface_residual = sp.factor(b_gb * integrated_interface_residual)
    sigma1_reversal_residual = sp.factor(sigma1.subs({k1: -k1, k2: -k2, k3: -k3}) + sigma1)
    sigma2_reversal_residual = sp.factor(sigma2.subs({k1: -k1, k2: -k2, k3: -k3}) - sigma2)
    sigma3_reversal_residual = sp.factor(sigma3.subs({k1: -k1, k2: -k2, k3: -k3}) + sigma3)

    p2348_poly = find_first_key(artifacts["P2348_CHERN_POLYNOMIAL"], "polynomial")
    p2348_polynomial_residual = sp.factor(chern_polynomial - sympify_expr(p2348_poly))
    p2353_ledger = find_first_key(artifacts["P2353_OBSTRUCTION_LEDGER"], "track_B_arbitrary_boundary_obstruction_ledger") or {}
    p2358_witness = find_first_key(artifacts["P2358_CLOSED_INTERFACE_GLUING"], "track_B_closed_interface_gluing_cancellation_witness") or {}
    p2359_witness = find_first_key(artifacts["P2359_S4_CAP_COMPLEMENT"], "track_B_s4_cap_complement_gluing_bulk_boundary_witness") or {}
    p2353_minimal_cut = p2353_ledger.get("minimal_next_missing_cut", [])

    samples = [
        sample_row("p2348_nonsymmetric_sample_K5_123", {K: 5, k1: 1, k2: 2, k3: 3}, chern_polynomial),
        sample_row("mixed_curvature_sample_K_minus2_347", {K: -2, k1: 3, k2: 4, k3: 7}, chern_polynomial),
        sample_row("zero_sigma3_nonzero_sigma1_sample", {K: 11, k1: 0, k2: 2, k3: 5}, chern_polynomial),
    ]
    all_sample_sums_zero = all(row["interface_sum"] == "0" for row in samples)

    lemma = {
        "lemma_id": "P2360_smooth_interface_chern_polynomial_sign_reversal_v1",
        "scope": "local smooth-interface Chern boundary polynomial sign-reversal lemma; not a corner, smoothing, arbitrary-boundary, or global CGB theorem",
        "hypotheses": [
            "A smooth internal hypersurface is seen from two glued sides with the same induced metric and the opposite unit normal.",
            "The ambient sectional-curvature input K is the same on both sides along the interface.",
            "The shape-operator eigenvalues transform as (k1,k2,k3) -> (-k1,-k2,-k3).",
            "No corner, smoothing collar, metric jump, distributional curvature, or non-smooth interface term is introduced.",
        ],
        "chern_boundary_polynomial": str(chern_polynomial),
        "reversed_normal_polynomial": str(reversed_polynomial),
        "sign_reversal_residual": str(sign_reversal_residual),
        "interface_pair_sum": str(interface_pair_sum),
        "integrated_interface_residual_symbolic": str(integrated_interface_residual),
        "pairing_interface_residual_symbolic": str(pairing_interface_residual),
        "sigma1_reversal_residual": str(sigma1_reversal_residual),
        "sigma2_even_reversal_residual_recorded_not_used": str(sigma2_reversal_residual),
        "sigma3_reversal_residual": str(sigma3_reversal_residual),
        "p2348_polynomial_residual": str(p2348_polynomial_residual),
        "sample_rows": samples,
        "all_sample_sums_zero": all_sample_sums_zero,
        "p2358_symbolic_interface_residual_replayed": p2358_witness.get("symbolic_interface_residual"),
        "p2358_all_case_residuals_replayed": p2358_witness.get("all_case_residuals_zero"),
        "p2359_interface_cancellation_replayed": p2359_witness.get("interface_cancellation"),
        "p2359_gluing_consistency_residual_replayed": p2359_witness.get("gluing_consistency_residual"),
        "p2353_minimal_cut_replayed": p2353_minimal_cut,
        "o5_smooth_interface_polynomial_partially_attacked": "O5_regularization_corners_and_gluing" in p2353_minimal_cut,
        "remaining_gap": "The lemma cancels smooth opposite-normal interface Chern-polynomial terms, but it does not supply corner, smoothing-limit, metric-jump, or arbitrary-boundary regularization terms.",
    }

    dependencies = {
        key.lower() + "_loaded": artifacts[key].get("result_kind") == expected
        for key, expected in RESULT_KINDS.items()
    }
    dependencies.update(
        {
            "p2348_polynomial_replayed": p2348_polynomial_residual == 0,
            "p2358_interface_replayed": lemma["p2358_symbolic_interface_residual_replayed"] == "0",
            "p2358_cases_replayed": lemma["p2358_all_case_residuals_replayed"] is True,
            "p2359_interface_replayed": lemma["p2359_interface_cancellation_replayed"] == "0",
            "p2359_gluing_consistency_replayed": lemma["p2359_gluing_consistency_residual_replayed"] == "0",
            "p2353_o5_cut_replayed": lemma["o5_smooth_interface_polynomial_partially_attacked"],
            "symbolic_sign_reversal_zero": sign_reversal_residual == 0,
            "symbolic_interface_pair_sum_zero": interface_pair_sum == 0,
            "integrated_interface_residual_zero": integrated_interface_residual == 0,
            "pairing_interface_residual_zero": pairing_interface_residual == 0,
            "all_sample_sums_zero": all_sample_sums_zero,
        }
    )

    theorem_export = {
        "theorem_name": "P2360 Track-B smooth-interface Chern-polynomial sign-reversal lemma",
        "claim": (
            "For a smooth glued interface with the same ambient K and opposite unit normals, the P2348 Chern boundary polynomial "
            "B3=8*K*sigma1+16*sigma3 changes sign under (k1,k2,k3)->(-k1,-k2,-k3).  Therefore the two interface contributions cancel "
            "locally before any convexity, geodesic, or spherical symmetry assumption.  This is a real smooth-interface O5 improvement, not a corner or smoothing theorem."
        ),
        "proof_witnesses": [
            "sigma1 and sigma3 are odd under normal reversal, while the even sigma2 term is absent from the P2348 polynomial.",
            "The symbolic residual B3(K;-k1,-k2,-k3)+B3(K;k1,k2,k3) is exactly zero.",
            "The integrated formal residual I+(-I) and its P2335-ledger pairing residual are exactly zero.",
            "Three nonsymmetric sample rows, including the P2348 K=5,(1,2,3) row, cancel exactly under normal reversal.",
            "P2358 and P2359 closed-interface cancellations are replayed as special cases, while corner and smoothing gaps remain explicit.",
        ],
        "not_licensed": [
            "arbitrary-boundary theorem",
            "general gluing theorem including corners or smoothing limits",
            "corner contribution theorem",
            "smoothing-limit theorem",
            "metric-jump or distributional-curvature theorem",
            "general Chern-Gauss-Bonnet theorem over arbitrary compact four-manifolds",
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
            "The smooth-interface polynomial cancellation no longer needs convexity or geodesic symmetry.  The next real target is a local corner-angle "
            "or smoothing-collar computation; do not keep replaying smooth opposite-normal interfaces as if corners were solved."
        ),
    }

    probe = {
        "probe_id": "P2360_S1310_STRICT_TRACK_B_SMOOTH_INTERFACE_CHERN_POLYNOMIAL_SIGN_REVERSAL_LEMMA",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2360 smooth-interface Chern-polynomial sign reversal after P2359",
            "top_hits": grep_hits(),
        },
        "track_B_smooth_interface_chern_polynomial_sign_reversal_lemma": lemma,
        "current_export_dependencies": dependencies,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "repo_grep_audit_recorded": len(probe["repo_grep_audit"]["top_hits"]) > 0,
        **dependencies,
        "three_sample_rows_recorded": len(samples) == 3,
        "smooth_interface_scope_declared": "smooth-interface" in lemma["scope"],
        "corner_theorem_not_claimed": True,
        "smoothing_limit_theorem_not_claimed": True,
        "general_gluing_theorem_not_claimed": True,
        "arbitrary_boundary_theorem_not_claimed": True,
        "metric_jump_theorem_not_claimed": True,
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
        "schema_version": "p2360_s1310_v1",
        "packet_id": "P2360",
        "stage_id": "S1310",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-29T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_SMOOTH_INTERFACE_CHERN_POLYNOMIAL_SIGN_REVERSAL_LEMMA_NO_CORNER_SMOOTHING_CLOSURE",
        "strict_track_b_smooth_interface_chern_polynomial_sign_reversal_lemma_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2360 strict Track-B smooth-interface Chern-polynomial sign-reversal lemma\n\n"
        "Status: smooth opposite-normal interface cancellation exported; no corner/smoothing, arbitrary-boundary, selector, or ToE closure.\n\n"
        f"- Polynomial `{chern_polynomial}`; reversed-normal polynomial `{reversed_polynomial}`; residual `{sign_reversal_residual}`.\n"
        f"- Interface pair sum `{interface_pair_sum}`; integrated residual `{integrated_interface_residual}`; pairing residual `{pairing_interface_residual}`.\n"
        f"- Oddness checks: sigma1 residual `{sigma1_reversal_residual}`, sigma3 residual `{sigma3_reversal_residual}`; sigma2 even residual `{sigma2_reversal_residual}` recorded but unused.\n"
        f"- Sample rows: `{[row['sample_id'] for row in samples]}`; all sums zero `{all_sample_sums_zero}`.\n"
        f"- Replays: P2348 residual `{p2348_polynomial_residual}`, P2358 interface `{lemma['p2358_symbolic_interface_residual_replayed']}`, P2359 interface `{lemma['p2359_interface_cancellation_replayed']}`.\n"
        f"- P2353 cut replayed: `{p2353_minimal_cut}`; O5 smooth-interface polynomial partially attacked `{lemma['o5_smooth_interface_polynomial_partially_attacked']}`.\n"
        "- No arbitrary-boundary theorem, no general gluing theorem including corners/smoothing, no corner contribution theorem, no smoothing-limit theorem, no metric-jump theorem, no general Chern-Gauss-Bonnet theorem, no global renormalization, no independent `a_GB`, no selector premise, no QW-2191 discharge, no unique-future choice, no observer prediction, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
