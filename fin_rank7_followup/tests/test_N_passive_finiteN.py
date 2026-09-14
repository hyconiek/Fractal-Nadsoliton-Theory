import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
R=json.loads((ROOT/'results/R7P-105_112_passive_finiteN.json').read_text())

def test_hodge_ranks_and_orientation():
    n=R['R7P-105']['numeric']; assert (n['tree_rank'],n['cycle_rank'])==(11,55); assert n['orientation_tree_residual']==0; assert n['orientation_cycle_residual']==0

def test_28_count_is_explicitly_numerical():
    d=R['R7P-106']; assert d['distinct_positive_count_tol_1e-8']==28; assert 'NUMERICAL_ONLY' in d['status']

def test_finiteN_generator_noise_scope():
    d=R['R7P-107']; assert d['state_count']==78; assert d['nonuniform_qv_vs_A_over_6N_norm']>0.7; assert d['mean_zero_orthonormal_mode_variance']=='1/(12N)'

def test_ou_nongaussian_correction():
    rows=R['R7P-108']['rows']; assert all(abs(x['exact_excess_kurtosis'] + 2/x['N'])<1e-15 for x in rows)

def test_visible_memory_degree_five():
    d=R['R7P-109_110']; assert d['sum_visible_residue_ranks']==5; assert d['pole_groups'][-1]['residue_rank']==0

def test_gain_source_firewall():
    assert R['R7P-112']['gain_source_status'].startswith('OPEN')
