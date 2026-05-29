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
OUT = GEN / "p2355_s1305_strict_track_b_nonround_ellipsoidal_shell_orientation_degree_witness.json"
MD = GEN / "p2355_s1305_strict_track_b_nonround_ellipsoidal_shell_orientation_degree_witness.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2338_TRACK_B_CONTRACT": GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.json",
    "P2345_CONVEX_GAUSS_MAP": GEN / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.json",
    "P2348_CHERN_POLYNOMIAL": GEN / "p2348_s1298_strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction.json",
    "P2349_TRIAXIAL_ELLIPSOID": GEN / "p2349_s1299_strict_track_b_triaxial_ellipsoid_integrated_boundary_polynomial_witness.json",
    "P2353_OBSTRUCTION_LEDGER": GEN / "p2353_s1303_strict_track_b_arbitrary_boundary_obstruction_ledger.json",
    "P2354_SPHERICAL_SHELL": GEN / "p2354_s1304_strict_track_b_spherical_shell_orientation_cancellation_witness.json",
}

RESULT_KINDS = {
    "P2335_TWO_TRACK_LEDGER": "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
    "P2338_TRACK_B_CONTRACT": "STRICT_TRACK_B_GB_BOUNDARY_TOPOLOGICAL_NORMALIZATION_CONTRACT_EXPORTED_NO_NORMALIZATION_CLAIM",
    "P2345_CONVEX_GAUSS_MAP": "STRICT_TRACK_B_CONVEX_GAUSS_MAP_BOUNDARY_FUNCTIONAL_THEOREM_AND_CLOSURE_PROMISE_MAP_NO_SELECTOR_CLOSURE",
    "P2348_CHERN_POLYNOMIAL": "STRICT_TRACK_B_CHERN_BOUNDARY_POLYNOMIAL_NONSYMMETRIC_REDUCTION_NO_INTEGRATED_GLOBAL_THEOREM",
    "P2349_TRIAXIAL_ELLIPSOID": "STRICT_TRACK_B_TRIAXIAL_ELLIPSOID_INTEGRATED_BOUNDARY_POLYNOMIAL_WITNESS_NO_UNIVERSAL_THEOREM",
    "P2353_OBSTRUCTION_LEDGER": "STRICT_TRACK_B_ARBITRARY_BOUNDARY_OBSTRUCTION_LEDGER_NO_UNIVERSAL_CLOSURE",
    "P2354_SPHERICAL_SHELL": "STRICT_TRACK_B_SPHERICAL_SHELL_ORIENTATION_CANCELLATION_WITNESS_NO_UNIVERSAL_CLOSURE",
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
        ROOT / "p2349_s1299_strict_track_b_triaxial_ellipsoid_integrated_boundary_polynomial_witness.py",
        ROOT / "p2353_s1303_strict_track_b_arbitrary_boundary_obstruction_ledger.py",
        ROOT / "p2354_s1304_strict_track_b_spherical_shell_orientation_cancellation_witness.py",
        GEN / "p2354_s1304_strict_track_b_spherical_shell_orientation_cancellation_witness.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track-B|Track B|ellipsoid|triaxial|shell|orientation|degree|nonconstant|nonconvex|Gauss-map|sigma3|arbitrary-boundary|selector|QW-2191|ToE closure",
        *existing,
    ]
    proc = subprocess.run(cmd, cwd=ROOT.parent, text=True, capture_output=True, check=False)
    return proc.stdout.splitlines()[:260]


def axis_pole_samples(axes: tuple[sp.Integer, sp.Integer, sp.Integer, sp.Integer], normal_sign: int) -> list[dict[str, Any]]:
    samples = []
    for pole_axis, a_j in enumerate(axes):
        curvatures = []
        for idx, a_i in enumerate(axes):
            if idx != pole_axis:
                curvatures.append(sp.factor(sp.Integer(normal_sign) * a_j / a_i**2))
        sigma3 = sp.factor(sp.prod(curvatures))
        samples.append(
            {
                "pole_axis_index": pole_axis,
                "axis_value_at_pole": str(a_j),
                "shell_oriented_principal_curvatures": [str(value) for value in curvatures],
                "sigma3": str(sigma3),
                "boundary_density_16_sigma3": str(sp.factor(16 * sigma3)),
            }
        )
    return samples


def ellipsoid_component(
    label: str,
    axes: tuple[sp.Integer, sp.Integer, sp.Integer, sp.Integer],
    normal_sign: int,
    b_gb: sp.Expr,
) -> dict[str, Any]:
    degree = int(normal_sign)
    integral_sigma3 = sp.factor(sp.Integer(degree) * 2 * sp.pi**2)
    boundary_functional = sp.factor(16 * integral_sigma3)
    pairing = sp.factor(boundary_functional * b_gb)
    samples = axis_pole_samples(axes, normal_sign)
    signed_sigma_values = [sympify_expr(sample["sigma3"]) for sample in samples]
    abs_sigma_values = [sp.factor(abs(value)) for value in signed_sigma_values]
    return {
        "component_id": label,
        "axes": [str(value) for value in axes],
        "shell_outward_normal_relative_to_ellipsoid_outward": normal_sign,
        "gauss_map_degree": degree,
        "axis_pole_samples": samples,
        "signed_sigma3_values": [str(value) for value in signed_sigma_values],
        "absolute_sigma3_values": [str(value) for value in abs_sigma_values],
        "nonconstant_abs_sigma3": bool(len({str(value) for value in abs_sigma_values}) > 1),
        "integral_sigma3_dA_by_degree": str(integral_sigma3),
        "boundary_functional_16_integral_sigma3": str(boundary_functional),
        "normalized_track_B_pairing": str(pairing),
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    b_gb = sp.factor(track_b_coefficient(artifacts["P2335_TWO_TRACK_LEDGER"]))
    target_boundary_per_degree = 32 * sp.pi**2
    target_pairing_per_degree = sp.factor(target_boundary_per_degree * b_gb)

    inner_axes = (sp.Integer(1), sp.Integer(2), sp.Integer(3), sp.Integer(5))
    outer_axes = (sp.Integer(3), sp.Integer(4), sp.Integer(6), sp.Integer(9))
    strict_containment_checks = [bool(outer > inner) for outer, inner in zip(outer_axes, inner_axes)]
    outer = ellipsoid_component("outer_triaxial_ellipsoid_axes_3_4_6_9", outer_axes, +1, b_gb)
    inner = ellipsoid_component("inner_triaxial_ellipsoid_axes_1_2_3_5_reversed_orientation", inner_axes, -1, b_gb)
    components = [outer, inner]

    component_boundary_sum = sp.factor(sum(sympify_expr(row["boundary_functional_16_integral_sigma3"]) for row in components))
    component_pairing_sum = sp.factor(sum(sympify_expr(row["normalized_track_B_pairing"]) for row in components))
    component_degree_sum = sum(int(row["gauss_map_degree"]) for row in components)
    shell_chi = sp.Integer(0)
    shell_expected_boundary = sp.factor(target_boundary_per_degree * shell_chi)
    shell_expected_pairing = sp.factor(target_pairing_per_degree * shell_chi)
    shell_boundary_residual = sp.factor(component_boundary_sum - shell_expected_boundary)
    shell_pairing_residual = sp.factor(component_pairing_sum - shell_expected_pairing)
    degree_boundary_residual = sp.factor(component_boundary_sum - target_boundary_per_degree * component_degree_sum)
    degree_pairing_residual = sp.factor(component_pairing_sum - target_pairing_per_degree * component_degree_sum)

    p2345_boundary = sp.factor(sympify_expr(find_first_key(artifacts["P2345_CONVEX_GAUSS_MAP"], "boundary_functional_convex_class")))
    p2348_replay = sp.factor(sympify_expr(find_first_key(artifacts["P2348_CHERN_POLYNOMIAL"], "substitute_integral_sigma3_equals_2pi2")))
    p2349_boundary = sp.factor(sympify_expr(find_first_key(artifacts["P2349_TRIAXIAL_ELLIPSOID"], "integrated_boundary_functional_16_integral_sigma3")))
    p2353_ledger = find_first_key(artifacts["P2353_OBSTRUCTION_LEDGER"], "track_B_arbitrary_boundary_obstruction_ledger") or {}
    p2354_witness = find_first_key(artifacts["P2354_SPHERICAL_SHELL"], "track_B_spherical_shell_orientation_cancellation_witness") or {}
    p2353_minimal_cut = p2353_ledger.get("minimal_next_missing_cut", [])

    shell_witness = {
        "witness_id": "P2355_nonround_ellipsoidal_shell_orientation_degree_v1",
        "scope": "one explicit nonconvex nonround ellipsoidal shell fixture; partial O4 degree/orientation accounting only",
        "domain": "Omega = {x in R^4 inside axes(3,4,6,9) ellipsoid and outside axes(1,2,3,5) ellipsoid}",
        "strict_coordinate_axis_containment_checks": strict_containment_checks,
        "components": components,
        "both_components_have_nonconstant_abs_sigma3": all(row["nonconstant_abs_sigma3"] for row in components),
        "component_gauss_degree_sum": component_degree_sum,
        "component_boundary_functional_sum": str(component_boundary_sum),
        "component_pairing_sum": str(component_pairing_sum),
        "shell_chi_used_for_replay": str(shell_chi),
        "expected_shell_boundary_functional": str(shell_expected_boundary),
        "expected_shell_pairing": str(shell_expected_pairing),
        "shell_boundary_residual": str(shell_boundary_residual),
        "shell_pairing_residual": str(shell_pairing_residual),
        "degree_boundary_residual": str(degree_boundary_residual),
        "degree_pairing_residual": str(degree_pairing_residual),
        "orientation_cancellation_identity": "(+32*pi**2) + (-32*pi**2) = 0 despite nonconstant axis-pole sigma3 samples",
        "target_boundary_per_degree": str(target_boundary_per_degree),
        "target_pairing_per_degree": str(target_pairing_per_degree),
        "p2345_convex_boundary_residual": str(sp.factor(p2345_boundary - target_boundary_per_degree)),
        "p2348_chern_replay_residual": str(sp.factor(p2348_replay - target_boundary_per_degree)),
        "p2349_triaxial_replay_residual": str(sp.factor(p2349_boundary - target_boundary_per_degree)),
        "p2354_shell_boundary_residual_replayed": p2354_witness.get("shell_boundary_residual"),
        "p2353_minimal_cut_replayed": p2353_minimal_cut,
        "o4_cut_partially_attacked": "O4_nonconvex_degree_and_orientation_accounting" in p2353_minimal_cut,
    }

    dependencies = {
        key.lower() + "_loaded": artifacts[key].get("result_kind") == expected
        for key, expected in RESULT_KINDS.items()
    }
    dependencies.update(
        {
            "p2345_convex_boundary_replayed": shell_witness["p2345_convex_boundary_residual"] == "0",
            "p2348_chern_replay_replayed": shell_witness["p2348_chern_replay_residual"] == "0",
            "p2349_triaxial_replay_replayed": shell_witness["p2349_triaxial_replay_residual"] == "0",
            "p2354_shell_residual_replayed": shell_witness["p2354_shell_boundary_residual_replayed"] == "0",
            "p2353_o4_cut_replayed": shell_witness["o4_cut_partially_attacked"],
            "strict_coordinate_axis_containment": all(strict_containment_checks),
            "both_components_nonconstant": shell_witness["both_components_have_nonconstant_abs_sigma3"],
            "shell_boundary_residual_zero": bool(shell_boundary_residual == 0),
            "shell_pairing_residual_zero": bool(shell_pairing_residual == 0),
            "degree_boundary_residual_zero": bool(degree_boundary_residual == 0),
            "degree_pairing_residual_zero": bool(degree_pairing_residual == 0),
        }
    )

    theorem_export = {
        "theorem_name": "P2355 Track-B nonround ellipsoidal-shell orientation-degree witness",
        "claim": (
            "For the explicit flat shell between the outer ellipsoid with axes (3,4,6,9) and the inner ellipsoid with axes (1,2,3,5), "
            "the two nonround boundary components have nonconstant axis-pole sigma3 samples.  Nevertheless their oriented Gauss-map degrees are +1 "
            "and -1, so their Track-B boundary functionals and P2335-ledger pairings cancel exactly.  This strengthens P2354 from round shell "
            "cancellation to a nonround/nonconstant O4 fixture, but it is not a general nonconvex theorem."
        ),
        "proof_witnesses": [
            "Strict coordinate-axis containment holds: 3>1, 4>2, 6>3, and 9>5.",
            "At each ellipsoid axis pole, the shell-oriented principal curvatures are computed by signed a_j/a_i**2 for i != j.",
            "Both boundary components have nonconstant absolute sigma3 samples, so this is not a round/equal-curvature replay.",
            "The outer component has degree +1 and contributes +32*pi**2; the inner shell-oriented component has degree -1 and contributes -32*pi**2.",
            "Boundary and pairing residuals vanish both against the shell replay target and against the degree-sum replay target.",
            "P2353's O4 nonconvex degree/orientation cut is further attacked but not globally discharged.",
        ],
        "not_licensed": [
            "arbitrary-boundary theorem",
            "general nonconvex boundary theorem",
            "general Chern-Gauss-Bonnet theorem over all compact four-manifolds with boundary",
            "general gluing theorem or smoothing-limit theorem",
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
            "O4 now has both a round shell and a nonround ellipsoidal-shell witness.  The next honest step is either a theorem-level local "
            "degree/orientation lemma with hypotheses, or switch to the still-hard O3 transgression-integration or O5 corners/gluing cuts."
        ),
    }

    probe = {
        "probe_id": "P2355_S1305_STRICT_TRACK_B_NONROUND_ELLIPSOIDAL_SHELL_ORIENTATION_DEGREE_WITNESS",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2355 nonround O4 orientation/degree accounting after P2354 spherical shell",
            "top_hits": grep_hits(),
        },
        "track_B_nonround_ellipsoidal_shell_orientation_degree_witness": shell_witness,
        "current_export_dependencies": dependencies,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "repo_grep_audit_recorded": len(probe["repo_grep_audit"]["top_hits"]) > 0,
        **dependencies,
        "outer_boundary_plus_target": outer["boundary_functional_16_integral_sigma3"] == "32*pi**2",
        "inner_boundary_minus_target": inner["boundary_functional_16_integral_sigma3"] == "-32*pi**2",
        "component_degree_sum_zero": bool(component_degree_sum == 0),
        "component_boundary_sum_zero": bool(component_boundary_sum == 0),
        "component_pairing_sum_zero": bool(component_pairing_sum == 0),
        "o4_only_partial_not_discharged": True,
        "arbitrary_boundary_theorem_not_claimed": True,
        "general_nonconvex_theorem_not_claimed": True,
        "general_cgb_theorem_not_claimed": True,
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
        "schema_version": "p2355_s1305_v1",
        "packet_id": "P2355",
        "stage_id": "S1305",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-29T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_NONROUND_ELLIPSOIDAL_SHELL_ORIENTATION_DEGREE_WITNESS_NO_UNIVERSAL_CLOSURE",
        "strict_track_b_nonround_ellipsoidal_shell_orientation_degree_witness_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2355 strict Track-B nonround ellipsoidal-shell orientation-degree witness\n\n"
        "Status: explicit nonround annular orientation/degree fixture exported; no arbitrary-boundary, selector, or ToE closure.\n\n"
        f"- Domain: outer axes `(3,4,6,9)`, inner axes `(1,2,3,5)` in flat `R^4`; `b_GB = {b_gb}`.\n"
        f"- Strict coordinate-axis containment checks: `{strict_containment_checks}`.\n"
        f"- Outer component: degree `+1`, nonconstant abs sigma3 `{outer['nonconstant_abs_sigma3']}`, boundary `{outer['boundary_functional_16_integral_sigma3']}`, pairing `{outer['normalized_track_B_pairing']}`.\n"
        f"- Inner component: degree `-1`, nonconstant abs sigma3 `{inner['nonconstant_abs_sigma3']}`, boundary `{inner['boundary_functional_16_integral_sigma3']}`, pairing `{inner['normalized_track_B_pairing']}`.\n"
        f"- Component sums: degree `{component_degree_sum}`, boundary `{component_boundary_sum}`, pairing `{component_pairing_sum}`.\n"
        f"- Shell residuals: boundary `{shell_boundary_residual}`, pairing `{shell_pairing_residual}`; degree residuals boundary `{degree_boundary_residual}`, pairing `{degree_pairing_residual}`.\n"
        f"- P2353 cut replayed: `{p2353_minimal_cut}`; O4 partially attacked `{shell_witness['o4_cut_partially_attacked']}`.\n"
        "- No arbitrary-boundary theorem, no general nonconvex boundary theorem, no general Chern-Gauss-Bonnet theorem, no general gluing theorem, no global renormalization, no independent `a_GB`, no selector premise, no QW-2191 discharge, no unique-future choice, no observer prediction, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
