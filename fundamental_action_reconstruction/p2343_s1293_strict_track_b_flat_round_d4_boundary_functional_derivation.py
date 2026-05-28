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
OUT = GEN / "p2343_s1293_strict_track_b_flat_round_d4_boundary_functional_derivation.json"
MD = GEN / "p2343_s1293_strict_track_b_flat_round_d4_boundary_functional_derivation.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2338_TRACK_B_CONTRACT": GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.json",
    "P2341_D4_BOUNDARY_FIXTURE": GEN / "p2341_s1291_strict_track_b_d4_boundary_correction_fixture_witness.json",
    "P2342_D4_TO_S4_GLUING": GEN / "p2342_s1292_strict_track_b_d4_to_s4_boundary_gluing_compatibility_witness.json",
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
        ROOT / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.py",
        ROOT / "p2341_s1291_strict_track_b_d4_boundary_correction_fixture_witness.py",
        ROOT / "p2342_s1292_strict_track_b_d4_to_s4_boundary_gluing_compatibility_witness.py",
        GEN / "p2341_s1291_strict_track_b_d4_boundary_correction_fixture_witness.md",
        GEN / "p2342_s1292_strict_track_b_d4_to_s4_boundary_gluing_compatibility_witness.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track B|D4|D\\^4|boundary correction|boundary functional|flat|round|S3|S\\^3|Euler|chi|a_GB|QW-2191",
        *existing,
    ]
    proc = subprocess.run(cmd, cwd=ROOT.parent, text=True, capture_output=True, check=False)
    return proc.stdout.splitlines()[:160]


def track_b_coefficient(p2335: dict[str, Any]) -> sp.Expr:
    probe = p2335.get("strict_two_track_renormalization_ledger_probe", {})
    ledger = probe.get("two_track_ledger", {})
    track_b = ledger.get("track_B_gb_topological_counterterm_ledger", {})
    raw = track_b.get("ledger_coefficient_b_GB_topological", "0")
    return sp.sympify(raw, locals={"pi": sp.pi, "log": sp.log, "ln": sp.log})


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    p2335 = artifacts["P2335_TWO_TRACK_LEDGER"]
    p2338 = artifacts["P2338_TRACK_B_CONTRACT"]
    p2341 = artifacts["P2341_D4_BOUNDARY_FIXTURE"]
    p2342 = artifacts["P2342_D4_TO_S4_GLUING"]

    radius = sp.symbols("R", positive=True)
    b_gb = sp.factor(track_b_coefficient(p2335))
    per_euler_pairing = sp.factor(sp.simplify(32 * sp.pi**2 * b_gb))

    dimension_boundary = sp.Integer(3)
    chi_flat_round_d4 = sp.Integer(1)
    bulk_gb_flat_value = sp.Integer(0)
    s3_area = sp.factor(2 * sp.pi**2 * radius**3)
    principal_curvature = sp.factor(1 / radius)
    elementary_symmetric_sigma3 = sp.factor(principal_curvature**3)
    boundary_density_normalization = sp.Integer(16)
    local_boundary_density = sp.factor(boundary_density_normalization * elementary_symmetric_sigma3)
    boundary_functional = sp.factor(sp.simplify(local_boundary_density * s3_area))
    target_topological_number = sp.factor(32 * sp.pi**2 * chi_flat_round_d4)
    boundary_functional_residual = sp.factor(sp.simplify(boundary_functional - target_topological_number))
    radius_derivative_residual = sp.factor(sp.diff(boundary_functional, radius))
    unit_radius_value = sp.factor(boundary_functional.subs(radius, 1))

    normalized_flat_round_pairing = sp.factor(sp.simplify(b_gb * (bulk_gb_flat_value + boundary_functional)))
    normalized_pairing_residual = sp.factor(
        sp.simplify(normalized_flat_round_pairing - per_euler_pairing * chi_flat_round_d4)
    )
    unit_radius_pairing = sp.factor(normalized_flat_round_pairing.subs(radius, 1))

    p2341_probe = p2341.get("strict_track_b_d4_boundary_correction_fixture_witness_probe", {})
    p2341_witness = p2341_probe.get("track_B_D4_boundary_fixture_witness", {})
    p2342_probe = p2342.get("strict_track_b_d4_to_s4_boundary_gluing_compatibility_witness_probe", {})
    p2342_witness = p2342_probe.get("track_B_D4_to_S4_boundary_gluing_witness", {})
    p2342_gluing = p2342_witness.get("gluing_model", {})

    two_copy_boundary_functional = sp.factor(sp.simplify(2 * boundary_functional))
    s4_topological_number_from_p2342 = sp.sympify(
        p2342_gluing.get("S4_topological_number", "0"), locals={"pi": sp.pi, "log": sp.log}
    )
    two_copy_residual_against_p2342_s4 = sp.factor(
        sp.simplify(two_copy_boundary_functional - s4_topological_number_from_p2342)
    )

    derivation_steps = [
        {
            "step": "flat_bulk_zero",
            "formula": "For the Euclidean metric on D^4_R, the bulk GB density fixture value is 0.",
            "symbolic_value": str(bulk_gb_flat_value),
        },
        {
            "step": "round_boundary_geometry",
            "formula": "Boundary S^3_R has area 2*pi**2*R**3 and principal curvatures k1=k2=k3=1/R.",
            "s3_area": str(s3_area),
            "sigma3_shape_operator": str(elementary_symmetric_sigma3),
        },
        {
            "step": "normalized_boundary_density",
            "formula": "B_flat_round := 16 * sigma3(K) * dA under the P2338/P2341 convention.",
            "local_density_without_dA": str(local_boundary_density),
        },
        {
            "step": "radius_cancellation",
            "formula": "Integral_{S3_R} 16*(1/R**3)*dA = 16*(1/R**3)*(2*pi**2*R**3).",
            "integrated_boundary_functional": str(boundary_functional),
            "d_dR_boundary_functional": str(radius_derivative_residual),
        },
    ]

    flat_round_class = {
        "class_id": "flat_round_four_ball_radius_R_positive",
        "scope": "one-parameter flat round D^4_R boundary class only; not arbitrary boundaries",
        "radius_symbol": "R > 0",
        "bulk_gb_fixture_value": str(bulk_gb_flat_value),
        "boundary_manifold": "round S^3_R",
        "boundary_area": str(s3_area),
        "shape_operator_principal_curvatures": [str(principal_curvature)] * int(dimension_boundary),
        "sigma3_shape_operator": str(elementary_symmetric_sigma3),
        "boundary_density_normalization_factor": str(boundary_density_normalization),
        "boundary_functional_B_flat_round_R": str(boundary_functional),
        "target_topological_number_32pi2_chi": str(target_topological_number),
        "boundary_functional_residual": str(boundary_functional_residual),
        "radius_derivative_residual": str(radius_derivative_residual),
        "unit_radius_value": str(unit_radius_value),
        "normalized_track_B_pairing": str(normalized_flat_round_pairing),
        "normalized_track_B_pairing_unit_radius": str(unit_radius_pairing),
        "normalized_pairing_residual": str(normalized_pairing_residual),
        "two_copy_boundary_functional": str(two_copy_boundary_functional),
        "two_copy_residual_against_p2342_S4": str(two_copy_residual_against_p2342_s4),
    }

    dependencies = {
        "p2335_two_track_ledger_loaded": p2335.get("result_kind")
        == "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
        "p2338_contract_loaded": p2338.get("result_kind")
        == "STRICT_TRACK_B_GB_BOUNDARY_TOPOLOGICAL_NORMALIZATION_CONTRACT_EXPORTED_NO_NORMALIZATION_CLAIM",
        "p2341_d4_boundary_fixture_loaded": p2341.get("result_kind")
        == "STRICT_TRACK_B_D4_NONZERO_BOUNDARY_FIXTURE_WITNESS_NO_UNIVERSAL_BOUNDARY_THEOREM",
        "p2342_d4_to_s4_gluing_loaded": p2342.get("result_kind")
        == "STRICT_TRACK_B_D4_TO_S4_BOUNDARY_GLUING_COMPATIBILITY_WITNESS_NO_GENERAL_GLUING_THEOREM",
        "p2341_boundary_value_rederived_at_unit_radius": p2341_witness.get("boundary_correction_fixture_value")
        == str(unit_radius_value),
        "p2341_pairing_rederived_at_unit_radius": p2341_witness.get("normalized_D4_pairing") == str(unit_radius_pairing),
        "p2342_two_copy_s4_number_rederived": two_copy_residual_against_p2342_s4 == 0,
    }

    theorem_export = {
        "theorem_name": "P2343 Track-B flat-round D4 boundary functional derivation",
        "claim": (
            "P2343 replaces the P2341 D^4 boundary fixture value with a one-parameter computation on flat round "
            "four-balls D^4_R.  Under the declared Track-B convention B_flat_round=16*sigma3(K)*dA, the round "
            "S^3_R boundary has area 2*pi**2*R**3 and sigma3(K)=1/R**3, so the boundary functional is "
            "32*pi**2 for every R>0.  The bulk GB value remains 0, the radius derivative is zero, and the "
            "Track-B pairing is 13152087*log(2)/10000000 with zero residual."
        ),
        "proof_witnesses": [
            "Flat Euclidean D^4_R has bulk GB fixture value 0 in this boundary-functional class.",
            "The round S^3_R shape operator has three equal principal curvatures 1/R.",
            "The product 16*(1/R**3)*(2*pi**2*R**3) simplifies to 32*pi**2.",
            "The derivative with respect to R is zero, so the value is radius-independent within this class.",
            "At R=1 the value and pairing exactly rederive P2341; two copies rederive the P2342 S4 number.",
        ],
        "not_licensed": [
            "boundary functional for non-round or curved-boundary four-manifolds",
            "general Chern-Gauss-Bonnet theorem over arbitrary compact four-manifolds with boundary",
            "general boundary gluing theorem for arbitrary interfaces",
            "connected-sum theorem beyond previously declared fixture bookkeeping",
            "independent a_GB measurement separate from the P2335 ledger coefficient",
            "bulk EOM GB force or EOM-only GB lift",
            "full/global renormalization closure",
            "selector premise or QW-2191 selector discharge",
            "G1/G3 update",
            "ToE closure",
        ],
        "next_honest_step": (
            "Still leave selector work to the end.  The next proof/computational Track-B step should either extend "
            "the boundary functional from flat round D^4_R to a second explicit boundary class with nonconstant "
            "principal curvatures, or freeze this as a one-parameter class derivation below any universal boundary theorem."
        ),
    }

    probe = {
        "probe_id": "P2343_S1293_STRICT_TRACK_B_FLAT_ROUND_D4_BOUNDARY_FUNCTIONAL_DERIVATION",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2342 next step: more proof/computational Track-B boundary functional derivation",
            "top_hits": grep_hits(),
        },
        "track_B_flat_round_D4_boundary_functional_derivation": {
            "ledger_coefficient_b_GB": str(b_gb),
            "per_euler_pairing": str(per_euler_pairing),
            "derivation_steps": derivation_steps,
            "flat_round_boundary_class": flat_round_class,
        },
        "current_export_dependencies": dependencies,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "repo_grep_audit_recorded": len(probe["repo_grep_audit"]["top_hits"]) > 0,
        "p2335_two_track_ledger_loaded": dependencies["p2335_two_track_ledger_loaded"],
        "p2338_contract_loaded": dependencies["p2338_contract_loaded"],
        "p2341_d4_boundary_fixture_loaded": dependencies["p2341_d4_boundary_fixture_loaded"],
        "p2342_d4_to_s4_gluing_loaded": dependencies["p2342_d4_to_s4_gluing_loaded"],
        "p2341_boundary_value_rederived_at_unit_radius": dependencies["p2341_boundary_value_rederived_at_unit_radius"],
        "p2341_pairing_rederived_at_unit_radius": dependencies["p2341_pairing_rederived_at_unit_radius"],
        "p2342_two_copy_s4_number_rederived": dependencies["p2342_two_copy_s4_number_rederived"],
        "flat_bulk_gb_zero_declared": bulk_gb_flat_value == 0,
        "round_s3_area_formula_recorded": s3_area == 2 * sp.pi**2 * radius**3,
        "sigma3_shape_operator_recorded": elementary_symmetric_sigma3 == radius**-3,
        "boundary_functional_is_32pi2": boundary_functional_residual == 0,
        "radius_derivative_zero": radius_derivative_residual == 0,
        "normalized_pairing_residual_zero": normalized_pairing_residual == 0,
        "one_parameter_class_scope_declared": "one-parameter flat round" in flat_round_class["scope"],
        "nonround_boundary_functional_not_claimed": True,
        "general_cgb_theorem_not_claimed": True,
        "general_boundary_gluing_theorem_not_claimed": True,
        "independent_a_gb_not_claimed": True,
        "bulk_eom_gb_force_not_claimed": True,
        "no_full_global_renormalization_claimed": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_g1_g3_update_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2343_s1293_v1",
        "packet_id": "P2343",
        "stage_id": "S1293",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_FLAT_ROUND_D4_BOUNDARY_FUNCTIONAL_DERIVATION_NO_UNIVERSAL_BOUNDARY_THEOREM",
        "strict_track_b_flat_round_d4_boundary_functional_derivation_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2343 strict Track-B flat-round D4 boundary functional derivation\n\n"
        "Status: one-parameter flat round D4_R boundary functional derived; no universal boundary theorem.\n\n"
        f"- `b_GB = {b_gb}`.\n"
        f"- Boundary class: Euclidean `D^4_R`, `R>0`, boundary `S^3_R`, area `{s3_area}`.\n"
        f"- Shape data: principal curvatures `{principal_curvature}` repeated 3 times, `sigma3(K) = {elementary_symmetric_sigma3}`.\n"
        f"- Boundary functional: `16*sigma3(K)*dA -> {boundary_functional}`; target `32*pi**2*chi = {target_topological_number}`; residual `{boundary_functional_residual}`.\n"
        f"- Radius derivative residual: `{radius_derivative_residual}`; unit-radius value: `{unit_radius_value}`.\n"
        f"- Normalized Track-B pairing: `{normalized_flat_round_pairing}`; residual `{normalized_pairing_residual}`.\n"
        f"- Two-copy check against P2342 S4 topological number: `{two_copy_boundary_functional}` with residual `{two_copy_residual_against_p2342_s4}`.\n"
        "- No non-round boundary functional, no general Chern-Gauss-Bonnet theorem, no independent `a_GB`, no full/global renormalization, no selector premise, no QW-2191 discharge, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
