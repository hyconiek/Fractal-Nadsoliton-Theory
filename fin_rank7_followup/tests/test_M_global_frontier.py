import json
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def data(): return json.loads((ROOT/'results/R7P-097_104_global_frontier.json').read_text())

def test_rational_upper_witness_is_strict():
    d=data()['R7P-097_099']; assert d['upper_witness']['strict_negative']; assert d['bracket']==['2.8934','3.71835']

def test_g4_uniqueness_remains_unresolved():
    assert data()['R7P-100_101']['unique_global_orbit_at_g4']=='UNRESOLVED'

def test_declared_flow_is_not_physical_promotion():
    d=data()['R7P-103']; assert d['proof_level'].startswith('NUMERICAL'); assert 'theta_dot' in d['law']
    assert {x['classification'] for x in d['branches']}=={'uniform','localized_orbit_rep'}
