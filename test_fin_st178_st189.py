#!/usr/bin/env python3
"""Regression and live checks for FIN programs ST178--ST189."""

from __future__ import annotations

import hashlib
import inspect
import json
import math
from pathlib import Path

import numpy as np

import fin_st178_st189_research as research
from fin_st01_st15_research import strict_operator
from fin_st189_external_record_validator import synthetic_self_test


ROOT=Path(__file__).resolve().parent
DATA=json.loads((ROOT/"FIN_ST178_ST189_Results.json").read_text(encoding="utf-8"))


def file_sha(path:Path)->str:return hashlib.sha256(path.read_bytes()).hexdigest()


def test_programs_and_statuses_present():
    assert all(f"ST{k}" in DATA and DATA[f"ST{k}"]["status"] for k in range(178,190))


def test_all_packet_hashes():
    for k in range(178,190):
        d=DATA[f"ST{k}"]
        assert file_sha(ROOT/d["packet_file"])==d["packet_sha256"]


def test_st178_hard_hierarchy_ranks_and_kernels():
    d=DATA["ST178"]
    assert d["layer_dimensions"]==[12,6,3]
    assert d["composite_ranks"]==[12,6,3]
    assert d["invisible_kernel_dimensions"]==[0,6,9]
    c=np.array(d["composite_map"])
    assert np.linalg.matrix_rank(c)==3
    assert np.allclose(c.sum(axis=0),1)


def test_st178_indistinguishable_states_are_deeply_distinct():
    row=DATA["ST178"]["indistinguishable_deep_states"]
    assert row["coarse_p"]==row["coarse_q"]
    assert row["deep_effect_p"]==1
    assert row["deep_effect_q"]==0


def test_st178_visible_effect_dimension_live():
    d=research.st178_irreversible_carriers()
    c=np.array(d["composite_map"])
    assert np.linalg.matrix_rank(c.T)==3
    assert 12-np.linalg.matrix_rank(c)==9


def test_st179_covariance_and_state_dependence():
    d=DATA["ST179"]
    assert d["C12_covariance_error"]<1e-12
    rows={r["state"]:r for r in d["rows"]}
    assert rows["uniform"]["L_commutant_dimension"]==144
    assert rows["localized_vertex_0"]["L_commutant_dimension"]==22
    assert rows["asymmetric_full_support"]["L_commutant_dimension"]==12
    assert rows["asymmetric_full_support"]["joint_A_L_commutant_dimension"]==1


def test_st179_uniform_state_cannot_select():
    row=DATA["ST179"]["rows"][0]
    assert row["number_distinct_L_values"]==1


def test_st180_all_signed_slices_are_certified():
    d=DATA["ST180"]
    assert d["uniform_folds"]==30
    assert d["signed_slice_targets"]==60
    assert d["certified_continued_states"]==60
    assert all(r["certificate"] and r["certificate"]["margin"]>0 for r in d["rows"])


def test_st180_continued_states_are_nonuniform_and_accurate():
    rows=DATA["ST180"]["rows"]
    assert max(r["residual_inf"] for r in rows)<1e-12
    assert min(r["nonuniform_norm"] for r in rows)>1e-5
    assert {r["slice_epsilon"] for r in rows}=={-.001,.001}


def test_st180_live_slice_residual():
    _,a,_=strict_operator();r=DATA["ST180"]["rows"][17]
    q0,k0,v=research.uniform_fold_seed(a,r["uniform_amplitude"],r["mode"])
    f,_=research.stationary_slice_float(np.array(r["continued_center"]),a,q0,v,r["slice_epsilon"])
    assert np.linalg.norm(f,np.inf)<1e-12


def test_st181_readout_is_stochastic_and_fidelity_normalized():
    d=DATA["ST181"];q=np.array(d["readout_y_given_error"])
    assert np.allclose(q.sum(axis=0),1)
    assert 0<=d["optimal_entanglement_fidelity"]<=1
    assert abs(sum(r["optimal_contribution"] for r in d["rows"])-d["optimal_entanglement_fidelity"])<1e-14


def test_st181_primal_dual_certificates():
    for r in DATA["ST181"]["rows"]:
        assert r["dual_slack_eigenvalues"][0]>-1e-12
        assert abs(r["dual_slack_eigenvalues"][1]-2*r["dual_t"])<1e-12
        assert abs(math.hypot(*r["primal_phase"])-1)<1e-12


def test_st182_stationary_gap_inversion_and_monotonicity():
    for row in DATA["ST182"]["rows"]:
        assert abs(row["dimensionless_field_from_gap"]-row["beta_h"])<1e-12
        gaps=[r["gap"] for r in row["trajectory"]]
        assert all(a<=b for a,b in zip(gaps,gaps[1:]))
        assert gaps[-1]<=row["stationary_gap"]+1e-15


def test_st182_zero_field_has_zero_gap_live():
    row=research.detailed_balance_row(0,1,[0,1])
    assert row["stationary_gap"]==0
    assert all(x["gap"]==0 for x in row["trajectory"])


def test_st183_two_estimators_recover_g():
    for row in DATA["ST183"]["tests"]:
        assert abs(row["g_from_vertex"]-row["g"])<1e-12
        assert abs(row["g_from_first_mode"]-row["g"])<1e-12


def test_st184_complete_343_cell_cover():
    d=DATA["ST184"]
    assert d["subdivisions_per_axis"]==7
    assert d["boxes"]==343
    assert d["passed_boxes"]==343
    assert len(d["rows"])==343
    assert d["minimum_margin"]>0
    assert d["global_halfwidth"]==1e-4


def test_st185_soft_layer_eigenvalues_and_visibility():
    d=DATA["ST185"]
    mu=np.array(d["single_layer_Fourier_eigenvalues"])
    expected=np.array([.75+.25*math.cos(2*math.pi*k/12) for k in range(7)])
    assert np.allclose(mu,expected)
    last=d["visibility_rows"][-1]
    assert abs(last["mode_attenuation"][6]-2**-12)<1e-15
    assert last["inverse_noise_amplification"]==4096


def test_st185_data_processing_monotonicity():
    rows=DATA["ST185"]["data_processing_audit"]
    for key in ("KL","total_variation"):
        vals=[r[key] for r in rows]
        assert all(a>=b-1e-14 for a,b in zip(vals,vals[1:]))
    assert DATA["ST185"]["hard_invisible_dimension"]==9


def test_st185_layer_zero_is_identity_visibility():
    row=DATA["ST185"]["visibility_rows"][0]
    assert row["mode_attenuation"]==[1.0]*7
    assert row["inverse_noise_amplification"]==1


def test_st186_noiseless_factor_is_exact():
    d=DATA["ST186"];r=d["rows"][0]
    assert abs(d["ideal_gap"]-4)<1e-12
    assert r["projector_distance"]==0
    assert r["multiplication_leakage_of_raw_low_basis"]<1e-12


def test_st186_noise_creates_nonclosure():
    rows=DATA["ST186"]["rows"][1:]
    assert all(r["projector_distance"]>0 for r in rows)
    assert all(r["multiplication_leakage_of_raw_low_basis"]>0 for r in rows)


def test_st187_pair_path_bounds_observed_coefficient():
    for r in DATA["ST187"]["rows"]:
        assert r["pair_path_upper_bound"]>=r["observed_Hellinger_coefficient"]
        assert r["observed_finite_rate"]>=r["pair_path_finite_rate"]
        assert r["bound_over_exact_ratio"]>=1


def test_st187_gap_is_large_by_n20():
    row=DATA["ST187"]["rows"][-1]
    assert row["events"]==20
    assert row["bound_over_exact_ratio"]>4000
    assert row["observed_finite_rate"]-row["pair_path_finite_rate"]>.4


def test_st188_global_swap_and_pairwise_ansatz():
    d=DATA["ST188"]
    assert d["global_register_SWAP_commutator_norm"]==0
    assert d["numerical_pairwise_linear_ansatz_nullity"]==0
    assert min(d["Gram_eigenvalues"])>200
    assert d["equal_pairwise_sum_commutator_norm"]>0


def test_st188_live_swap_product_commutes():
    _,a,_=strict_operator();h=np.zeros((16,16));h[:12,:12]=a;i=np.eye(16);ht=np.kron(h,i)+np.kron(i,h)
    swaps=[research.swap_bit_matrix(8,k,4+k) for k in range(4)]
    global_swap=swaps[0]@swaps[1]@swaps[2]@swaps[3]
    assert np.linalg.norm(ht@global_swap-global_swap@ht,np.inf)==0


def test_st189_self_test_blocks_physical_promotion_and_tampering():
    s=DATA["ST189"]["self_test"]
    assert s["test_passed"] and s["synthetic_only"] and s["no_external_record_supplied"]
    assert s["clean"]["structural_record_valid"]
    assert not s["clean"]["physical_execution_valid"]
    assert not s["tampered"]["structural_record_valid"]


def test_st189_live_self_test():
    assert synthetic_self_test()["test_passed"]


def test_recommendations_cover_st190_through_st202():
    recs=DATA["recommended_next_programs"]
    assert [r["id"] for r in recs]==[f"ST{k}" for k in range(190,203)]
    assert [r["priority"] for r in recs]==list(range(1,14))


def test_epistemic_boundary_preserves_guardrails():
    b=DATA["epistemic_boundary"]
    for token in ("Planck","QW-2191","legacy-to-strict","Standard Model","gravity","L_total","ToE"):
        assert token in b


if __name__=="__main__":
    passed=0
    for name,fn in sorted(inspect.getmembers(__import__(__name__),inspect.isfunction)):
        if name.startswith("test_"):fn();passed+=1
    print(f"{passed}/{passed} tests passed")
