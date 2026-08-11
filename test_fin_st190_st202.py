#!/usr/bin/env python3
"""Regression and live checks for FIN programs ST190--ST202."""

from __future__ import annotations

import hashlib
import inspect
import json
import math
from pathlib import Path

import numpy as np

import fin_st190_st202_research as research
from fin_st01_st15_research import strict_operator


ROOT=Path(__file__).resolve().parent
DATA=json.loads((ROOT/"FIN_ST190_ST202_Results.json").read_text(encoding="utf-8"))


def file_sha(path:Path)->str:return hashlib.sha256(path.read_bytes()).hexdigest()


def test_all_programs_statuses_and_packet_hashes():
    for k in range(190,203):
        d=DATA[f"ST{k}"];assert d["status"]
        assert file_sha(ROOT/d["packet_file"])==d["packet_sha256"]


def test_st190_heat_layers_are_stochastic_and_positive():
    d=DATA["ST190"]
    assert d["semigroup_composition_error"]<2e-13
    assert all(r["minimum_entry"]>=-1e-14 and r["maximum_row_sum_error"]<2e-13 for r in d["rows"])


def test_st190_visibility_matches_strict_fourier_classes():
    d=DATA["ST190"];lam=np.array(d["strict_Fourier_class_eigenvalues_k0_to_k6"])
    for r in d["rows"]:assert np.allclose(r["Fourier_visibility"],np.exp(-r["composite_parameter"]*lam))


def test_st190_base_schedule_is_explicitly_not_canonical():
    assert "conventions" in DATA["ST190"]["boundary"] and DATA["ST190"]["base_parameter_tau0"]==.25


def test_st191_minimal_sufficient_quotient_and_decoder():
    d=DATA["ST191"]
    assert d["minimal_quotient_outcomes"]==4 and d["exact_reconstruction"]
    assert d["minimal_proportionality_classes"]==[[0,1,2],[3,4,5],[6,7,8],[9,10,11]]


def test_st191_overmerge_loses_discrimination():
    d=DATA["ST191"]
    assert d["overmerged_total_variation_after"]<d["overmerged_total_variation_before"]


def test_st191_live_proportional_class_calculation():
    d=research.st191_blackwell_quotient()
    assert d["minimal_quotient_outcomes"]==4 and d["exact_reconstruction"]


def test_st192_uniform_fixed_point_and_instability():
    d=DATA["ST192"]
    assert d["uniform_rhs_norm"]<1e-14
    assert d["minimum_nonconstant_linear_growth_rate"]>0


def test_st192_vertices_are_stable_and_noise_selects_multiple_members():
    d=DATA["ST192"]
    assert d["vertex_local_stability_margin"]>2
    assert len(d["distinct_selected_vertices"])>=4
    assert min(r["final_max_probability"] for r in d["trajectories"])>.99


def test_st193_all_additional_slices_are_certified():
    d=DATA["ST193"]
    assert d["targets"]==120 and d["certified_targets"]==120
    assert all(r["certificate"] and r["certificate"]["margin"]>0 for r in d["rows"])


def test_st193_live_residual_for_one_branch():
    _,a,_=strict_operator();r=DATA["ST193"]["rows"][31]
    q0,_,v=research.uniform_fold_seed(a,r["uniform_amplitude"],r["mode"])
    f,_=research.stationary_slice_float(np.array(r["continued_center"]),a,q0,v,r["slice_epsilon"])
    assert np.linalg.norm(f,np.inf)<2e-12


def test_st194_errors_are_noncommuting_and_certificate_closes():
    d=DATA["ST194"]
    assert min(d["pairwise_commutator_norms"])>.1
    assert d["SVD_orientation"]>.999999
    assert d["primal_upper_gap"]<2e-14
    assert 0<=d["primal_entanglement_fidelity"]<=1


def test_st194_recovery_is_a_rotation():
    r=np.array(DATA["ST194"]["optimal_recovery_rotation"])
    assert np.allclose(r.T@r,np.eye(3)) and np.linalg.det(r)>.999999


def test_st195_dissipation_positive_and_convergent():
    rows=DATA["ST195"]["quasistatic_rows"]
    vals=[r["dissipation"] for r in rows]
    assert all(x>0 for x in vals)
    assert all(x>y for x,y in zip(vals,vals[1:]))
    assert vals[-1]<.005


def test_st195_target_state_normalized():
    d=DATA["ST195"];p0=d["target_preferred_probability"]
    assert 0<p0<1 and abs(p0+11*(p0-d["target_gap"])-1)<1e-12


def test_st196_same_hessian_different_response():
    rows=DATA["ST196"]["rows"]
    assert all(r["Hessian_at_zero_difference_from_A"]==0 for r in rows)
    assert len({r["twelfth_response"] for r in rows})==len(rows)


def test_st197_certified_and_failed_fixed_cover_bracket():
    d=DATA["ST197"]
    assert d["accepted_cover"]["passed"]==343 and d["accepted_cover"]["minimum_margin"]>0
    assert d["first_tested_rejected_cover"]["passed"]<343 and d["first_tested_rejected_cover"]["minimum_margin"]<0
    assert abs(d["method_bracket_width"]-5e-8)<1e-15


def test_st198_leading_sample_scaling_is_quartic_per_layer():
    rows=DATA["ST198"]["rows"]
    ratios=[rows[i+1]["Cramer_Rao_samples_for_sd_0_01"]/rows[i]["Cramer_Rao_samples_for_sd_0_01"] for i in range(5,12)]
    assert all(3.99<r<4.01 for r in ratios)


def test_st198_fisher_information_formula_live():
    d=DATA["ST198"]
    for r in d["rows"]:
        s=r["attenuation"];a=s*d["deep_mode_amplitude"]
        assert abs(r["Fisher_information_per_sample"]-4*s*s/(1-4*a*a))<1e-14


def test_st199_projected_generators_form_exact_factor_numerically():
    for r in DATA["ST199"]["rows"]:
        assert r["generated_algebra_dimension"]==4 and r["commutant_dimension"]==4
        assert r["X_involution_error"]<1e-13 and r["Z_involution_error"]<1e-13 and r["anticommutator_error"]<1e-13
        assert r["X_sign_gap"]>0 and r["Z_odd_sign_gap"]>0


def test_st200_all_finite_intervals_are_valid_and_narrow():
    for r in DATA["ST200"]["rows"]:
        lo,hi=r["Hellinger_coefficient_interval"];rl,rh=r["finite_rate_interval"]
        assert 0<lo<=hi<=1 and rl<=rh
        assert r["interval_width"]<1e-14
        from fractions import Fraction
        assert Fraction(r["exact_rational_width"])>0


def test_st200_exact_rational_enclosures_present():
    for r in DATA["ST200"]["rows"]:
        assert "/" in r["exact_rational_lower"] and "/" in r["exact_rational_upper"]


def test_st201_both_static_spans_have_certified_full_rank():
    d=DATA["ST201"]
    assert d["pairwise_SWAP_interval_lower_singular_bound"]>10
    assert d["expanded_static_generator_count"]==60 and d["expanded_interval_lower_singular_bound"]>6


def test_st201_exact_interval_radius_is_recorded():
    assert DATA["ST201"]["strict_interval_max_entry_radius"]>=0


def test_st202_refuses_local_physical_promotion():
    d=DATA["ST202"]
    assert not d["physical_execution_valid"] and not d["execution_attempted"]
    assert d["missing_external_atoms"]==d["required_external_atoms"]
    assert d["synthetic_self_test"]["test_passed"]


def test_recommendations_cover_st203_through_st217():
    recs=DATA["recommended_next_programs"]
    assert [r["id"] for r in recs]==[f"ST{k}" for k in range(203,218)]
    assert [r["priority"] for r in recs]==list(range(1,16))


def test_epistemic_boundary_preserves_project_guardrails():
    b=DATA["epistemic_boundary"]
    for token in ("Planck","QW-2191","legacy-to-strict","Standard Model","gravity","L_total","ToE"):
        assert token in b


if __name__=="__main__":
    passed=0
    for name,fn in sorted(inspect.getmembers(__import__(__name__),inspect.isfunction)):
        if name.startswith("test_"):fn();passed+=1
    print(f"{passed}/{passed} tests passed")
