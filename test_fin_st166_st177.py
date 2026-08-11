#!/usr/bin/env python3
"""Regression and live mathematical checks for FIN programs ST166--ST177."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import sympy as sp

import fin_st166_st177_research as research
from fin_st177_operational_validator import ATOMS, digest, run, validate


ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT / "FIN_ST166_ST177_Results.json").read_text(encoding="utf-8"))


def file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_all_programs_and_statuses_present():
    assert all(f"ST{k}" in DATA for k in range(166, 178))
    assert all(DATA[f"ST{k}"]["status"] for k in range(166, 178))


def test_st166_live_stochastic_left_inverses():
    d = research.st166_carrier_invariant()
    assert d["carrier_1_left_inverse_error"] == 0
    assert d["carrier_2_left_inverse_error"] == 0
    assert d["maximum_probability_table_mismatch"] < 1e-15


def test_st166_shannon_entropy_is_not_tomographically_complete():
    c = DATA["ST166"]["same_entropy_counterexample"]
    assert abs(c["entropy_p"] - c["entropy_q"]) < 1e-15
    assert c["vertex_effect_probabilities_p"] != c["vertex_effect_probabilities_q"]


def test_st166_packet_hash():
    d = DATA["ST166"]
    assert file_sha(ROOT / d["packet_file"]) == d["packet_sha256"]


def test_st167_functional_commutant_no_go():
    rows = DATA["ST167"]["rows"]
    assert rows[0]["joint_commutant_dimension"] == 144
    assert all(r["joint_commutant_dimension"] >= 22 for r in rows)
    assert rows[1]["joint_commutant_dimension"] == 22
    assert DATA["ST167"]["strict_commutant_dimension"] == 22


def test_st168_exact_polynomial_and_root_count():
    rho = sp.symbols("rho")
    p = 3*rho**3 + 18*rho**2 - 284*rho + 24
    roots = sp.Poly(p, rho).intervals(eps=sp.Rational(1, 10**14))
    assert len(roots) == 3
    assert sum(1 for (lohi, _) in roots if lohi[0] >= 0) == 2
    assert DATA["ST168"]["positive_rho_intervals"] == [
        [float(lohi[0]), float(lohi[1])] for lohi, _ in roots if lohi[0] >= 0
    ]


def test_st168_constant_mode_is_never_singular_at_stationarity():
    rho = sp.symbols("rho")
    p = 3*rho**3 + 18*rho**2 - 284*rho + 24
    radial_num = 3*(rho**4 + 8*rho**3 + 344*rho**2 - 608*rho + 16)
    assert sp.gcd(p, radial_num) == 1
    assert sp.resultant(p, radial_num, rho) == DATA["ST168"]["constant_mode_exact_resultant"]


def test_st168_fold_count_and_intervals():
    d = DATA["ST168"]
    assert len(d["positive_spectral_classes"]) == 6
    assert len(d["fold_rows"]) == 18
    assert d["geometric_fold_points_null_line_mod_sign"] == 30
    assert d["normalized_augmented_fold_roots_including_v_sign"] == 60
    assert all(r["kappa_interval"][0] <= r["kappa_interval"][1] for r in d["fold_rows"])


def test_st169_identity_limit_and_nonorthogonal_degradation():
    rows = DATA["ST169"]["rows"]
    assert abs(rows[0]["optimal_entanglement_fidelity"] - 1) < 1e-15
    assert rows[-1]["optimal_entanglement_fidelity"] < rows[0]["optimal_entanglement_fidelity"]
    assert all(0 <= r["optimal_entanglement_fidelity"] <= 1 for r in rows)


def test_st169_live_branch_formula():
    details, fidelity = research.binary_phase_recovery(math.pi / 3)
    assert abs(fidelity - sum(r["optimal_contribution"] for r in details)) < 1e-15
    assert all(r["coherence_abs"] <= r["identity_weight"] + r["unitary_weight"] + 1e-15 for r in details)


def test_st170_fisher_cost_is_positive_and_monotone():
    rows = DATA["ST170"]["rows"]
    costs = [r["minimum_geodesic_action_T1"] for r in rows]
    assert all(x > 0 for x in costs)
    assert all(a < b for a, b in zip(costs, costs[1:]))
    d, action = research.fisher_rao_cost(.25, 2)
    assert abs(action - d*d/4) < 1e-15


def test_st170_does_not_claim_strict_selector():
    assert "No strict selector" in DATA["ST170"]["boundary"]


def test_st171_exact_jet_obstruction_scope():
    d = DATA["ST171"]
    assert d["jet_order"] == 11
    assert "twelfth-order" in d["boundary"]


def test_st172_extended_cover_complete():
    d = DATA["ST172"]
    assert d["boxes"] == 125
    assert d["passed_boxes"] == d["boxes"]
    assert len(d["rows"]) == 125
    assert d["minimum_margin"] > 0
    assert d["global_halfwidth"] == 7.5e-5


def test_st172_packet_hash():
    d = DATA["ST172"]
    assert file_sha(ROOT / d["certificate_file"]) == d["certificate_sha256"]


def test_st173_exact_path_parity_formula():
    for eps in (.001, .02, .1):
        q = research.path_parity_error(eps, 3)
        direct = 3*eps*(1-eps)**2 + eps**3
        assert abs(q - direct) < 1e-15


def test_st173_majority_compression_improves_declared_rows():
    for row in DATA["ST173"]["rows"]:
        assert row["five_path_majority_error"] < row["three_edge_path_parity_error"]
        assert row["compression_gain"] > 1


def test_st174_ideal_factor_commutant_and_gap():
    d = DATA["ST174"]
    assert d["ideal_commutant_dimension"] == 4
    assert d["ideal_gap"] == 4
    assert sum(abs(x) < 1e-12 for x in d["ideal_spectrum"]) == 4


def test_st174_analytic_threshold_switches_where_declared():
    rows = DATA["ST174"]["rows"]
    assert all(r["analytic_four_dimensional_cluster_separated"] for r in rows[:-1])
    assert not rows[-1]["analytic_four_dimensional_cluster_separated"]
    assert rows[-2]["analytic_superoperator_perturbation_bound"] < 2
    assert rows[-1]["analytic_superoperator_perturbation_bound"] > 2


def test_st174_live_ideal_laplacian():
    x = np.array([[0, 1], [1, 0]], complex)
    z = np.diag([1, -1]).astype(complex)
    i = np.eye(2)
    vals = np.linalg.eigvalsh(research.commutator_laplacian([np.kron(x, i), np.kron(z, i)]))
    assert np.max(abs(vals[:4])) < 1e-12
    assert abs(vals[4] - 4) < 1e-12


def test_st175_collatz_interval_contains_float_and_sharpens_st162():
    d = DATA["ST175"]
    lo, hi = d["spectral_radius_interval"]
    assert lo <= d["floating_spectral_radius"] <= hi
    assert hi < DATA["ST162_reference_row_sum"]
    assert d["certified_pair_path_error_exponent_lower_bound"] > .32


def test_st175_all_component_intervals_intersect_certificate():
    lo, hi = DATA["ST175"]["spectral_radius_interval"]
    for a, b in DATA["ST175"]["component_ratio_intervals"]:
        assert lo <= a <= b <= hi


def test_st176_dimension_and_depth_tradeoff():
    d = DATA["ST176"]
    assert d["minimum_bath_qubits_from_dimension"] == math.ceil(math.log2(12)) == 4
    assert d["unused_computational_states"] == 4
    assert d["parallel_pairwise_SWAP_gates"] == 4
    assert d["two_qubit_circuit_depth"] == 1
    assert d["energy_conserving_global_SWAP_locality_qubits"] == 8


def test_st177_validator_accepts_structure_but_blocks_physics():
    d = DATA["ST177"]
    assert d["validation"]["structural_packet_valid"]
    assert not d["validation"]["physical_execution_valid"]
    assert d["validation"]["verdict"] == "MATHEMATICAL_PACKET_READY_PHYSICAL_EXECUTION_BLOCKED"
    assert d["all_eleven_deletions_detected"]


def test_st177_live_deletion_audit_and_hashes():
    result = run(write=False)
    assert result["all_eleven_deletions_detected"]
    assert len(result["deletion_audit"]) == len(ATOMS) == 11
    packet = result["packet"]
    record = packet["RECORD"]
    assert digest(record["raw_events"]) == record["raw_events_sha256"]
    assert digest(record["holdout"]) == record["holdout_sha256"]
    assert validate(packet)["structural_packet_valid"]


def test_all_specialized_packet_hashes():
    for key in ("ST166", "ST168", "ST169", "ST173", "ST174", "ST175", "ST176"):
        d = DATA[key]
        assert file_sha(ROOT / d["packet_file"]) == d["packet_sha256"]
    d = DATA["ST177"]
    assert file_sha(ROOT / d["packet_file"]) == d["packet_sha256_file"]
    assert file_sha(ROOT / d["validator_record_file"]) == d["validator_record_sha256"]


def test_recommendations_cover_st178_through_st189():
    recs = DATA["recommended_next_programs"]
    assert [r["id"] for r in recs] == [f"ST{k}" for k in range(178, 190)]
    assert [r["priority"] for r in recs] == list(range(1, 13))


def test_epistemic_boundary_preserves_project_guardrails():
    boundary = DATA["epistemic_boundary"]
    for token in ("QW-2191", "legacy-to-strict", "Standard Model", "gravity", "L_total", "ToE"):
        assert token in boundary
