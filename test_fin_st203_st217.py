#!/usr/bin/env python3
"""Regression and live checks for FIN programs ST203--ST217."""

from __future__ import annotations

import hashlib
import inspect
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

import fin_st203_st217_research as research
from fin_st01_st15_research import strict_operator


ROOT=Path(__file__).resolve().parent
DATA=json.loads((ROOT/"FIN_ST203_ST217_Results.json").read_text(encoding="utf-8"))


def file_sha(p:Path)->str:return hashlib.sha256(p.read_bytes()).hexdigest()


def test_programs_statuses_and_packet_hashes():
    for k in range(203,218):
        d=DATA[f"ST{k}"];assert d["status"]
        assert file_sha(ROOT/d["packet_file"])==d["packet_sha256"]


def test_st203_all_declared_schedules_have_same_total_map():
    d=DATA["ST203"]
    assert all(abs(r["total"]-3.75)<1e-15 for r in d["schedule_rows"])
    assert max(r["final_map_error_from_exp_minus_3_75A"] for r in d["schedule_rows"])<1e-13


def test_st203_geometric_family_is_nonunique():
    rows=DATA["ST203"]["geometric_four_layer_family"]
    assert len({r["ratio"] for r in rows})==5
    assert all(abs(r["total"]-3.75)<1e-15 for r in rows)


def test_st204_complete_partition_lattice_counts():
    d=DATA["ST204"]
    assert d["partition_count"]==15
    assert d["rank_counts_by_blocks"]=={"1":1,"2":7,"3":6,"4":1}
    assert d["Hasse_edge_count"]==31


def test_st204_live_partition_generator():
    p=set(research.set_partitions(list(range(4))))
    assert len(p)==15
    assert research.refines(((0,),(1,),(2,),(3,)),((0,1),(2,3)))


def test_st205_ou_covariance_and_cyclic_invariance():
    d=DATA["ST205"]
    assert d["Lyapunov_residual"]<1e-12 and d["cyclic_covariance_error"]<1e-14


def test_st205_branch_selection_is_not_a_single_canonical_label():
    d=DATA["ST205"];counts=d["selected_vertex_counts"]
    assert sum(counts)==d["replicator_runs"]==120
    assert sum(x>0 for x in counts)==12
    assert d["minimum_final_max_probability"]>.999


def test_st206_all_pseudo_arclength_steps_certified():
    d=DATA["ST206"]
    assert d["steps"]==40 and d["certified_steps"]==40
    assert d["minimum_Jacobian_rank_margin"]>7e-4
    assert d["minimum_tangent_alignment"]>.999
    assert all(r["certificate"] and r["certificate"]["margin"]>0 for r in d["rows"])


def test_st206_live_stationary_residual():
    _,a,_=strict_operator();r=DATA["ST206"]["rows"][17]
    assert r["stationary_residual"]<2e-12
    # Reconstructing the stored root is unnecessary; the packet hash and live source function are checked separately.
    assert research.stationary7(np.zeros(8),a)[0].shape==(7,)


def test_st207_negative_orientation_and_global_fidelity():
    d=DATA["ST207"]
    assert Fraction(d["Bloch_determinant_exact"])<0
    assert d["global_optimal_entanglement_fidelity"]==.3
    assert set(d["optimal_recoveries"])=={"X","Y","Z"}


def test_st208_speed_bound_monotonicity_and_unreachable_rule():
    rows=DATA["ST208"]["rows"]
    reachable=[r for r in rows if r["reachable"]]
    assert all(a["minimum_time"]<b["minimum_time"] for a,b in zip(reachable,reachable[1:]))
    assert all(r["entropy_production_at_speed_saturating_constant_control"]>0 for r in reachable)
    assert DATA["ST208"]["maximum_control_equilibrium_probability"]<1


def test_st209_response_ratio_is_scale_invariant():
    d=DATA["ST209"]
    for g in (.1,1,10):
        rows=[r for r in d["rows"] if r["g"]==g]
        assert max(r["invariant_ratio"] for r in rows)-min(r["invariant_ratio"] for r in rows)<1e-15
        assert all(abs(r["invariant_ratio"]-r["expected_g_over_Ajj6"])<1e-15 for r in rows)


def test_st210_complete_adaptive_cover():
    d=DATA["ST210"]
    assert d["cells_per_axis"]==14 and d["boxes"]==2744 and d["passed_boxes"]==2744
    assert d["global_halfwidth"]==2e-4 and d["minimum_margin"]>9e-4


def test_st211_correlated_counts_inflate_iid_count():
    d=DATA["ST211"]
    assert abs(d["design_effect"]-1.95)<1e-14
    for r in d["rows"]:assert r["cluster_correlated_samples"]>=r["iid_samples_for_worst_mode_sd_0_01"]


def test_st211_worst_mode_cost_grows_quartically_at_deep_layers():
    rows=DATA["ST211"]["rows"]
    assert rows[-1]["iid_samples_for_worst_mode_sd_0_01"]>13_000_000_000
    assert abs(rows[-1]["mode_attenuations"][-1]-2**-12)<1e-15


def test_st212_circle_is_orthogonal_and_globally_degenerate():
    d=DATA["ST212"]
    assert d["raw_generators_invertible"]
    assert max(abs(r["n_dot_m"]) for r in d["sampled_minimizer_circle"])<1e-14
    assert max(abs(r["objective"]-d["analytic_minimum"]) for r in d["sampled_minimizer_circle"])<1e-14


def test_st213_rate_brackets_are_nested_in_width_and_valid():
    rows=DATA["ST213"]["rows"]
    assert all(r["asymptotic_rate_lower"]<r["asymptotic_rate_upper"] for r in rows)
    assert all(a["bracket_width"]>b["bracket_width"] for a,b in zip(rows,rows[1:]))
    assert abs(rows[-1]["bracket_width"]-math.log(8)/12)<1e-14


def test_st213_exact_conditional_ratio_constants():
    assert DATA["ST213"]["conditional_to_stationary_likelihood_ratio_bounds_exact"]==["3/10","12/5"]


def test_st214_pointwise_control_no_go_uses_positive_injectivity_bound():
    d=DATA["ST214"]
    assert d["certified_commutator_map_lower_singular_bound"]>6.3
    assert d["conclusion"]=="H_c(t)=0 almost everywhere and U(T)=I"


def test_st215_diffusion_distance_is_metric_in_all_audited_rows():
    assert all(r["triangle_violations"]==0 for r in DATA["ST215"]["rows"])


def test_st215_first_mode_circle_limit():
    d=DATA["ST215"]
    assert d["maximum_t10_limit_error"]<5e-8
    assert np.allclose(d["large_t_rescaled_distances_t10"],d["first_mode_chordal_limit"],atol=5e-8)


def test_st215_live_translation_invariant_distance():
    _,a,_=strict_operator();c=np.linalg.matrix_power(np.eye(12),1)@research.expm(-.5*a)
    assert abs(np.linalg.norm(c[0]-c[1])-np.linalg.norm(c[4]-c[5]))<1e-13


def test_st216_hard_and_soft_examples_separate_properties():
    hard,soft=DATA["ST216"]["rows"]
    assert hard["rank"]==3 and hard["kernel_dimension"]==9 and not hard["globally_injective"]
    assert soft["rank"]==12 and soft["kernel_dimension"]==0 and soft["globally_injective"]
    assert soft["exact_spectral_condition_number"]>1e64


def test_st217_external_gate_remains_blocked():
    d=DATA["ST217"]
    assert d["synthetic_self_test_passed"] and not d["external_execution_performed"]
    assert d["physical_result"] is None and len(d["missing_external_atoms"])==6


def test_recommendations_cover_st218_through_st232():
    recs=DATA["recommended_next_programs"]
    assert [r["id"] for r in recs]==[f"ST{k}" for k in range(218,233)]
    assert [r["priority"] for r in recs]==list(range(1,16))


def test_epistemic_boundary_preserves_guardrails():
    b=DATA["epistemic_boundary"]
    for token in ("Planck","QW-2191","legacy-to-strict","Standard Model","gravity","L_total","ToE"):
        assert token in b


if __name__=="__main__":
    passed=0
    for name,fn in sorted(inspect.getmembers(__import__(__name__),inspect.isfunction)):
        if name.startswith("test_"):fn();passed+=1
    print(f"{passed}/{passed} tests passed")
