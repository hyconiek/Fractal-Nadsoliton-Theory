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
OUT = GEN / "p2354_s1304_strict_track_b_spherical_shell_orientation_cancellation_witness.json"
MD = GEN / "p2354_s1304_strict_track_b_spherical_shell_orientation_cancellation_witness.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2338_TRACK_B_CONTRACT": GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.json",
    "P2345_CONVEX_GAUSS_MAP": GEN / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.json",
    "P2348_CHERN_POLYNOMIAL": GEN / "p2348_s1298_strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction.json",
    "P2352_CONTROLLED_SYNTHESIS": GEN / "p2352_s1302_strict_track_b_controlled_convex_boundary_class_synthesis_theorem.json",
    "P2353_OBSTRUCTION_LEDGER": GEN / "p2353_s1303_strict_track_b_arbitrary_boundary_obstruction_ledger.json",
}

RESULT_KINDS = {
    "P2335_TWO_TRACK_LEDGER": "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
    "P2338_TRACK_B_CONTRACT": "STRICT_TRACK_B_GB_BOUNDARY_TOPOLOGICAL_NORMALIZATION_CONTRACT_EXPORTED_NO_NORMALIZATION_CLAIM",
    "P2345_CONVEX_GAUSS_MAP": "STRICT_TRACK_B_CONVEX_GAUSS_MAP_BOUNDARY_FUNCTIONAL_THEOREM_AND_CLOSURE_PROMISE_MAP_NO_SELECTOR_CLOSURE",
    "P2348_CHERN_POLYNOMIAL": "STRICT_TRACK_B_CHERN_BOUNDARY_POLYNOMIAL_NONSYMMETRIC_REDUCTION_NO_INTEGRATED_GLOBAL_THEOREM",
    "P2352_CONTROLLED_SYNTHESIS": "STRICT_TRACK_B_CONTROLLED_CONVEX_BOUNDARY_CLASS_SYNTHESIS_THEOREM_NO_UNIVERSAL_CLOSURE",
    "P2353_OBSTRUCTION_LEDGER": "STRICT_TRACK_B_ARBITRARY_BOUNDARY_OBSTRUCTION_LEDGER_NO_UNIVERSAL_CLOSURE",
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
        ROOT / "p2353_s1303_strict_track_b_arbitrary_boundary_obstruction_ledger.py",
        GEN / "p2353_s1303_strict_track_b_arbitrary_boundary_obstruction_ledger.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track-B|Track B|spherical shell|orientation|degree|nonconvex|Gauss-map|sigma3|arbitrary-boundary|selector|QW-2191|ToE closure",
        *existing,
    ]
    proc = subprocess.run(cmd, cwd=ROOT.parent, text=True, capture_output=True, check=False)
    return proc.stdout.splitlines()[:240]


def sphere_component(label: str, radius: sp.Rational, normal_sign: int, b_gb: sp.Expr) -> dict[str, str | int]:
    kappa = sp.factor(sp.Integer(normal_sign) / radius)
    sigma3 = sp.factor(kappa**3)
    area = sp.factor(2 * sp.pi**2 * radius**3)
    integral_sigma3 = sp.factor(sigma3 * area)
    boundary_functional = sp.factor(16 * integral_sigma3)
    pairing = sp.factor(boundary_functional * b_gb)
    gauss_degree = int(normal_sign)
    return {
        "component_id": label,
        "radius": str(radius),
        "outward_normal_relative_to_radial": normal_sign,
        "principal_curvatures": [str(kappa), str(kappa), str(kappa)],
        "sigma3": str(sigma3),
        "area": str(area),
        "integral_sigma3_dA": str(integral_sigma3),
        "gauss_map_degree": gauss_degree,
        "boundary_functional_16_integral_sigma3": str(boundary_functional),
        "normalized_track_B_pairing": str(pairing),
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    b_gb = sp.factor(track_b_coefficient(artifacts["P2335_TWO_TRACK_LEDGER"]))
    target_boundary_per_degree = 32 * sp.pi**2
    target_pairing_per_degree = sp.factor(target_boundary_per_degree * b_gb)

    inner_radius = sp.Rational(1, 2)
    outer_radius = sp.Integer(3)
    outer = sphere_component("outer_boundary_S3_R3_outward_radial", outer_radius, +1, b_gb)
    inner = sphere_component("inner_boundary_S3_R1_over_2_shell_outward_inward_radial", inner_radius, -1, b_gb)
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
    p2352_synthesis = find_first_key(artifacts["P2352_CONTROLLED_SYNTHESIS"], "track_B_controlled_convex_boundary_class_synthesis") or {}
    p2353_ledger = find_first_key(artifacts["P2353_OBSTRUCTION_LEDGER"], "track_B_arbitrary_boundary_obstruction_ledger") or {}
    p2353_minimal_cut = p2353_ledger.get("minimal_next_missing_cut", [])

    shell_witness = {
        "witness_id": "P2354_spherical_shell_orientation_cancellation_v1",
        "scope": "one explicit nonconvex annular flat R4 boundary fixture; partial O4 degree/orientation accounting only",
        "domain": "Omega = {x in R^4 : 1/2 <= |x| <= 3}",
        "inner_radius": str(inner_radius),
        "outer_radius": str(outer_radius),
        "components": components,
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
        "orientation_cancellation_identity": "(+32*pi**2) + (-32*pi**2) = 0",
        "target_boundary_per_degree": str(target_boundary_per_degree),
        "target_pairing_per_degree": str(target_pairing_per_degree),
        "p2345_convex_boundary_residual": str(sp.factor(p2345_boundary - target_boundary_per_degree)),
        "p2348_chern_replay_residual": str(sp.factor(p2348_replay - target_boundary_per_degree)),
        "p2352_boundary_residual_vector": p2352_synthesis.get("boundary_residual_vector", []),
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
            "p2352_boundary_residual_vector_zero": all(value == "0" for value in shell_witness["p2352_boundary_residual_vector"]),
            "p2353_o4_cut_replayed": shell_witness["o4_cut_partially_attacked"],
            "shell_boundary_residual_zero": shell_boundary_residual == 0,
            "shell_pairing_residual_zero": shell_pairing_residual == 0,
            "degree_boundary_residual_zero": degree_boundary_residual == 0,
            "degree_pairing_residual_zero": degree_pairing_residual == 0,
        }
    )

    theorem_export = {
        "theorem_name": "P2354 Track-B spherical-shell orientation cancellation witness",
        "claim": (
            "For the explicit flat annular domain Omega={1/2<=|x|<=3} in R^4, the outer boundary contributes +32*pi**2 "
            "and the inner boundary, with the shell outward orientation, contributes -32*pi**2.  Thus the Track-B boundary "
            "functional and P2335-ledger pairing cancel exactly.  This is a concrete O4 degree/orientation accounting witness, not a general nonconvex theorem."
        ),
        "proof_witnesses": [
            "The outer S3 component has kappa=1/3, sigma3=1/27, area=54*pi**2, integral sigma3 dA=2*pi**2, and boundary functional 32*pi**2.",
            "The inner S3 component has shell-outward normal opposite to radial, kappa=-2, sigma3=-8, area=pi**2/4, integral sigma3 dA=-2*pi**2, and boundary functional -32*pi**2.",
            "The component Gauss degrees +1 and -1 sum to 0, matching the shell chi replay used in this fixture.",
            "Boundary and pairing residuals vanish both against the shell replay target and against the degree-sum replay target.",
            "P2353's O4 nonconvex degree/orientation cut is partially attacked but not globally discharged.",
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
            "Continue O4 only by adding a genuine nonconstant-curvature nonconvex degree/orientation fixture, or switch to the still-hard "
            "O3 transgression-integration or O5 corners/gluing cuts. Do not treat the shell cancellation as arbitrary-boundary closure."
        ),
    }

    probe = {
        "probe_id": "P2354_S1304_STRICT_TRACK_B_SPHERICAL_SHELL_ORIENTATION_CANCELLATION_WITNESS",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2354 O4 orientation/degree accounting after P2353 obstruction ledger",
            "top_hits": grep_hits(),
        },
        "track_B_spherical_shell_orientation_cancellation_witness": shell_witness,
        "current_export_dependencies": dependencies,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "repo_grep_audit_recorded": len(probe["repo_grep_audit"]["top_hits"]) > 0,
        **dependencies,
        "outer_boundary_plus_target": outer["boundary_functional_16_integral_sigma3"] == "32*pi**2",
        "inner_boundary_minus_target": inner["boundary_functional_16_integral_sigma3"] == "-32*pi**2",
        "component_degree_sum_zero": component_degree_sum == 0,
        "component_boundary_sum_zero": component_boundary_sum == 0,
        "component_pairing_sum_zero": component_pairing_sum == 0,
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
        "schema_version": "p2354_s1304_v1",
        "packet_id": "P2354",
        "stage_id": "S1304",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-29T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_SPHERICAL_SHELL_ORIENTATION_CANCELLATION_WITNESS_NO_UNIVERSAL_CLOSURE",
        "strict_track_b_spherical_shell_orientation_cancellation_witness_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2354 strict Track-B spherical-shell orientation cancellation witness\n\n"
        "Status: explicit annular orientation/degree fixture exported; no arbitrary-boundary, selector, or ToE closure.\n\n"
        f"- Domain: `1/2 <= |x| <= 3` in flat `R^4`; `b_GB = {b_gb}`.\n"
        f"- Outer component: degree `+1`, boundary functional `{outer['boundary_functional_16_integral_sigma3']}`, pairing `{outer['normalized_track_B_pairing']}`.\n"
        f"- Inner component: degree `-1`, boundary functional `{inner['boundary_functional_16_integral_sigma3']}`, pairing `{inner['normalized_track_B_pairing']}`.\n"
        f"- Component sums: degree `{component_degree_sum}`, boundary `{component_boundary_sum}`, pairing `{component_pairing_sum}`.\n"
        f"- Shell residuals: boundary `{shell_boundary_residual}`, pairing `{shell_pairing_residual}`; degree residuals boundary `{degree_boundary_residual}`, pairing `{degree_pairing_residual}`.\n"
        f"- P2353 cut replayed: `{p2353_minimal_cut}`; O4 partially attacked `{shell_witness['o4_cut_partially_attacked']}`.\n"
        "- No arbitrary-boundary theorem, no general nonconvex boundary theorem, no general Chern-Gauss-Bonnet theorem, no general gluing theorem, no global renormalization, no independent `a_GB`, no selector premise, no QW-2191 discharge, no unique-future choice, no observer prediction, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
