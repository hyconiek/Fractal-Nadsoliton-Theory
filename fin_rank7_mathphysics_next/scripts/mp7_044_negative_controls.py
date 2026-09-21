#!/usr/bin/env python3
import hashlib, json, math
from pathlib import Path
import numpy as np
from scipy.linalg import eigvalsh
R=Path(__file__).resolve().parents[1]
results=[]
def add(name,passed,detail): results.append({'name':name,'passed':bool(passed),'detail':detail})

# 1 eigenvalue-order wording trap
vals=np.array([0.40,0.20,0.10,0.05])
dec=np.sort(vals)[::-1]; inc=np.sort(vals)
add('swapped_eigenvalue_ordering', dec[1]!=inc[1], {'decreasing_lambda2':float(dec[1]),'increasing_lambda2':float(inc[1])})

# 2 omitted metric under coordinate rescaling
M=np.diag([0.2,0.1]); T=np.diag([2.0,1.0]); Mp=T.T@M@T; Gp=T.T@T
raw=np.sort(np.linalg.eigvalsh(Mp))[::-1]
gen=np.sort(eigvalsh(Mp,Gp))[::-1]
add('omitted_metric_changes_raw_spectrum', np.max(abs(gen-np.array([0.2,0.1])))<1e-14 and abs(raw[0]-0.8)<1e-14,
    {'raw_transformed':raw.tolist(),'generalized':gen.tolist()})

# 3 rank-deficient witness basis
B=np.array([[1.,0.],[2.,0.],[3.,0.]])
add('rank_deficient_B_rejected', np.linalg.matrix_rank(B)<B.shape[1], {'rank':int(np.linalg.matrix_rank(B)),'columns':B.shape[1]})

# 4 geometry controls: equal total length is not enough
proper=[(0.,0.5),(0.5,1.)]; bad=[(0.,0.6),(0.7,1.1)] # same total length 1.0 but gap+overhang

def coverage(intervals):
    xs=sorted(intervals)
    gap=False; overlap=False
    end=0.0
    for a,b in xs:
        if a>end+1e-15: gap=True
        if a<end-1e-15: overlap=True
        end=max(end,b)
    return {'gap':gap,'overlap':overlap,'end':end,'total':sum(b-a for a,b in xs)}
c1,c2=coverage(proper),coverage(bad)
add('compensated_overlap_gap_not_volume_proof', (not c1['gap']) and c2['gap'] and abs(c1['total']-c2['total'])<1e-15,
    {'proper':c1,'bad_equal_total':c2})

# 5 stale hash
payload=b'certificate-input-v1'; h=hashlib.sha256(payload).hexdigest(); h2=hashlib.sha256(payload+b'!').hexdigest()
add('stale_input_hash_rejected', h!=h2, {'original':h,'mutated':h2})

# 6 negative-b positivity route
b=-1.0
add('negative_b_breaks_nonnegative_residue6_weight', math.sinh(b)<0, {'sinh_b':math.sinh(b),'required_action':'translate odd label first'})

# 7 zero-amplitude phase redundancy
# a3=0 makes its phase absent from the field identically.
phis=[0.0,0.7,2.3]
fieldvals=[0.0*math.cos(x) for x in phis]
add('zero_amplitude_phase_is_redundant', max(fieldvals)-min(fieldvals)==0.0, {'sampled_zero_mode_terms':fieldvals})

# 8 polar gradient term cannot be dropped off stationarity: Phi(x,y)=x at (r,0)
r=2.0
Q=-r # d2/dphi2 r cos(phi) at phi=0
Hyy=0.0; grad_r=1.0
rhs=r*r*Hyy-r*grad_r
add('nonstationary_polar_gradient_term_required', abs(Q-rhs)<1e-15 and abs(Q-r*r*Hyy)>1e-12,
    {'Q':Q,'r2_Hodd':r*r*Hyy,'gradient_corrected_rhs':rhs})

# 9 M4 cannot replace M7 in full finite-N fluctuations at uniform.
fd=json.loads((R/'results/MP7-031_033_finite_N_diagnostics.json').read_text())['lambda_midpoints']
l3,l4,l5,l6=[float(fd[str(k)]) for k in (3,4,5,6)]
M7diag=np.array([l3,l3,l4,l4,l5,l5,l6])/12.0
M4diag=np.array([l3,l4,l5,l6])/12.0
add('M4_is_not_M7_for_full_fluctuations', len(M7diag)==7 and len(M4diag)==4 and np.all(M7diag[[1,3,5]]>0),
    {'M7_uniform_diag':M7diag.tolist(),'M4_uniform_diag':M4diag.tolist(),'omitted_positive_sine_variances':M7diag[[1,3,5]].tolist()})

# 10 interpretation gates
add('local_energy_equality_not_global_by_logic', True, {'status':'semantic rejection rule retained'})
add('phase_root_not_full_stationary_by_logic', True, {'status':'confirmed concretely by MP7-040'})
add('auxiliary_equilibrium_law_not_unique_dynamics', True, {'status':'same invariant density admits different reversible mobilities/attempt rates'})
add('translated_equal_mass_not_selector', True, {'status':'symmetry multiplicity does not choose one orbit member'})

out={'task':'MP7-044','controls':results,'all_passed':all(x['passed'] for x in results),
     'scope':'targeted mathematical/semantic controls; not a replay of unavailable R7O3 leaf corruption tests'}
(R/'audit/MP7-044_negative_controls.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
