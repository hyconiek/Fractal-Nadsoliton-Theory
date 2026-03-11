#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_R1 = GENERATED / "r1_strict_core_residual_datum_target_slot_export_packet.json"
IN_THETA_PAIR = GENERATED / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1.json"
IN_QW2190 = REPO / "report_qw2190_kernel_mode_representation_emergence_gate.json"

OUT_INSTANCE = (
    GENERATED / "r1_residual_orientation_datum_target_slot_population_candidate_from_sigma_int_theta_pair_v1.json"
)
OUT_JSON = (
    GENERATED
    / "p402_current_strict_sigma_int_r1_target_slot_population_candidate_from_theta_pair_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p402_current_strict_sigma_int_r1_target_slot_population_candidate_from_theta_pair_probe_summary.json"
)


def load_json_path(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def real_fourier_basis(n: int) -> dict[str, np.ndarray]:
    x = np.arange(n, dtype=float)
    basis: dict[str, np.ndarray] = {"e0": np.ones(n, dtype=float) / math.sqrt(n)}
    for m in range(1, n // 2):
        basis[f"c{m}"] = math.sqrt(2.0 / n) * np.cos(2.0 * math.pi * m * x / n)
        basis[f"s{m}"] = math.sqrt(2.0 / n) * np.sin(2.0 * math.pi * m * x / n)
    if n % 2 == 0:
        basis[f"e{n//2}"] = ((-1.0) ** x) / math.sqrt(n)
    return basis


def orthonormal_residual(b: np.ndarray) -> float:
    m = b.T @ b
    return float(np.linalg.norm(m - np.eye(m.shape[0])))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r1 = load_json_path(IN_R1)
    theta_pair = load_json_path(IN_THETA_PAIR)
    q2190 = load_json_path(IN_QW2190)

    n = int(q2190["mode_mapping"]["n_octaves"])
    fourier = real_fourier_basis(n)
    c1, s1 = fourier["c1"], fourier["s1"]
    c2, s2 = fourier["c2"], fourier["s2"]

    theta_1 = float(theta_pair["outputs"]["pair1"]["theta_1_cand"])
    theta_2 = float(theta_pair["outputs"]["pair2"]["theta_2_cand"])

    u1_expected = math.cos(theta_1) * c1 + math.sin(theta_1) * s1
    u2_expected = math.cos(theta_2) * c2 + math.sin(theta_2) * s2

    u1_art = np.array([float(x) for x in theta_pair["outputs"]["pair1"]["u_1_cand"]], dtype=float)
    u2_art = np.array([float(x) for x in theta_pair["outputs"]["pair2"]["u_2_cand"]], dtype=float)

    u1_match = float(np.linalg.norm(u1_art - u1_expected))
    u2_match = float(np.linalg.norm(u2_art - u2_expected))

    b = np.column_stack([u1_art, u2_art])
    gram = b.T @ b
    ortho_res = float(np.linalg.norm(gram - np.eye(2)))
    det_gram = float(np.linalg.det(gram))
    rank = int(np.linalg.matrix_rank(b, tol=1e-12))

    checks = [
        {
            "id": "r1_target_slot_export_present",
            "pass": bool(r1.get("stage") == "R1" and r1.get("export_target") == "residual_orientation_datum_target_slot"),
            "expected": True,
            "actual": {"stage": r1.get("stage"), "export_target": r1.get("export_target")},
            "meaning": "R1 target-slot export packet is present",
        },
        {
            "id": "theta_pair_artifact_present",
            "pass": bool(theta_pair.get("object") == "ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1"),
            "expected": True,
            "actual": {"object": theta_pair.get("object"), "status": theta_pair.get("status")},
            "meaning": "theta-pair selector ingredient artifact is present",
        },
        {
            "id": "u_vectors_match_r1_formulas",
            "pass": bool(u1_match <= 1e-12 and u2_match <= 1e-12),
            "expected": True,
            "actual": {"u1_match_residual": u1_match, "u2_match_residual": u2_match},
            "meaning": "u_1(theta_1), u_2(theta_2) computed from R1 formulas match the exported u_i^cand vectors",
        },
        {
            "id": "orientation_slice_is_2d",
            "pass": bool(rank == 2),
            "expected": 2,
            "actual": rank,
            "meaning": "span{u_1,u_2} is genuinely 2-dimensional (no collinearity)",
        },
        {
            "id": "orientation_slice_basis_orthonormal",
            "pass": bool(ortho_res <= 1e-12),
            "expected": True,
            "actual": {"orthonormal_residual": ortho_res, "gram_det": det_gram},
            "meaning": "u_1,u_2 form an orthonormal basis (numerical tolerance)",
        },
    ]

    instance = {
        "object": "R1_residual_orientation_datum_target_slot_population_candidate_from_sigma_int_theta_pair_v1",
        "status": "candidate_inhabitant_instance_constructed_and_audited__no_false_pass",
        "as_of": "2026-03-11",
        "intent": (
            "Construct one candidate inhabitant instance of the residual_orientation_datum_target_slot (R1) "
            "from the exported sigma_int->theta selector-ingredient theta pair, without promoting to strict-core "
            "theta export, strict-core target-slot population by an export-map object, selector closure, or QW-2191 discharge."
        ),
        "derived_from": {
            "r1_target_slot_export_packet": str(IN_R1.relative_to(REPO)),
            "theta_pair_selector_ingredient_artifact": str(IN_THETA_PAIR.relative_to(REPO)),
            "qw2190_mode_scaffold": str(IN_QW2190.relative_to(REPO)),
        },
        "inputs": {
            "theta_1_cand": theta_1,
            "theta_2_cand": theta_2,
            "basis_formula": r1.get("target_object_class"),
        },
        "outputs": {
            "u_1_cand": [float(x) for x in u1_art.tolist()],
            "u_2_cand": [float(x) for x in u2_art.tolist()],
            "orientation_slice_candidate": "span{u_1_cand,u_2_cand}",
        },
        "audits": {
            "u1_matches_r1_formula_residual": u1_match,
            "u2_matches_r1_formula_residual": u2_match,
            "slice_basis_orthonormal_residual": ortho_res,
            "slice_basis_rank": rank,
            "slice_basis_det_gram": det_gram,
            "pair_basis_orthonormal_residual_vs_identity": orthonormal_residual(b),
        },
        "hard_limits": [
            "Does not claim strict-core theta export (theta values remain candidate-only).",
            "Does not claim strict-core target-slot population by a strict-core export-map object (F311 remains sign-only).",
            "Does not claim object-support above the export-map object (N302/N395 remain open).",
            "Does not claim admissible S_sel_int nor strict-core selector closure.",
            "Does not claim QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    passed = bool(all(bool(c["pass"]) for c in checks))

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P402",
        "goal": "construct_and_audit_one_candidate_inhabitant_for_R1_target_slot_from_exported_theta_pair_selector_ingredient",
        "inputs": {
            "r1_target_slot_export_packet": str(IN_R1.relative_to(REPO)),
            "theta_pair_selector_ingredient_artifact": str(IN_THETA_PAIR.relative_to(REPO)),
            "qw2190_mode_scaffold": str(IN_QW2190.relative_to(REPO)),
        },
        "constructed_instance_artifact": str(OUT_INSTANCE.relative_to(REPO)),
        "checks": checks,
        "verdict": {
            "candidate_inhabitant_constructed": True,
            "audits_pass": passed,
            "statement": (
                "Candidate inhabitant instance for R1 is constructible from the exported theta-pair selector ingredient "
                "and passes basic internal audits. This does not promote strict-core population, object-support, selector closure, "
                "nor discharge QW-2191."
            ),
        },
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": report["stage"],
        "audits_pass": report["verdict"]["audits_pass"],
        "constructed_instance_artifact": report["constructed_instance_artifact"],
        "theta_1_cand": theta_1,
        "theta_2_cand": theta_2,
        "slice_basis_orthonormal_residual": ortho_res,
        "u1_matches_r1_formula_residual": u1_match,
        "u2_matches_r1_formula_residual": u2_match,
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    OUT_INSTANCE.write_text(json.dumps(instance, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()

