from __future__ import annotations
import json, itertools, math
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from .model import feature_spaces,d12_actions

ROOT=Path(__file__).resolve().parents[1]
W,A,L,X,C,A7=feature_spaces()

# R7P-009
coordinate={
 'column_order':['3cos','3sin','4cos','4sin','5cos','5sin','6alt'],
 'C4_X7_indices':[0,2,4,6],
 'laplacian_sector_eigenvalues':[float(x) for x in L],
 'column_squared_norms':[float(x) for x in np.diag(X.T@X)],
 'zero_column_mean_inf':float(np.max(np.abs(X.sum(axis=0)))),
 'A7_factorization_fro':float(np.linalg.norm(A7-X@X.T)),
 'rank_A7':int(np.linalg.matrix_rank(A7,tol=1e-11)),
 'signed_real_space_example_A7_row0':[float(x) for x in A7[0]],
}
(ROOT/'results/R7P-009_coordinate_mapping.json').write_text(json.dumps(coordinate,indent=2)+'\n')

# R7P-010
acts=d12_actions(X)
maxeq=max(np.linalg.norm(P@X-X@T,np.inf) for P,T in acts.values())
maxorth=max(np.linalg.norm(T.T@T-np.eye(7),np.inf) for P,T in acts.values())
closure=True
mats=[P for P,T in acts.values()]
for P,Q in itertools.product(mats,mats):
    R=P@Q
    if not any(np.array_equal(R,S) for S in mats): closure=False;break
d12={'count':len(acts),'max_feature_equivariance_inf':float(maxeq),'max_induced_orthogonality_inf':float(maxorth),'permutation_group_closure':closure,
     'warning':'reflection-fixed subspaces are not a fundamental-domain exhaustion; stabilizers here are only numerical unless separately certified.'}
(ROOT/'results/R7P-010_D12_action.json').write_text(json.dumps(d12,indent=2)+'\n')

# R7P-015 exact row norm formula (numerical display)
R2=(L[3]+L[4]+L[5])/6+L[6]/12
R=math.sqrt(R2)
domain={'R_squared_formula':'(lambda3+lambda4+lambda5)/6 + lambda6/12','R_squared':float(R2),'R':R,
        'stationary_ball':'||theta7|| <= g R','negative_sublevel_ball':'Phi_g(theta)<=0 => ||theta7|| <= 2 g R'}
(ROOT/'results/R7P-015_compact_domain_constants.json').write_text(json.dumps(domain,indent=2)+'\n')

# R7P-016: numerical exposed-support census for positive C4 rays.
# J and s coordinates differ by positive diagonal scaling, so positive normal cones have same support combinatorics.
points=C.copy(); n=12; d=4
supports=[]
for mask in range(1,1<<n):
    S=[i for i in range(n) if mask>>i&1]
    # Variables v[0:4], m, delta. Max delta.
    c=np.zeros(6);c[5]=-1
    Aeq=[];beq=[]
    row=np.zeros(6);row[:4]=1;Aeq.append(row);beq.append(1)
    for i in S:
        row=np.zeros(6);row[:4]=points[i];row[4]=-1;Aeq.append(row);beq.append(0)
    Aub=[];bub=[]
    for j in range(n):
        if j in S:continue
        row=np.zeros(6);row[:4]=points[j];row[4]=-1;row[5]=1
        Aub.append(row);bub.append(0)
    res=linprog(c,A_ub=np.array(Aub) if Aub else None,b_ub=np.array(bub) if Aub else None,
                A_eq=np.array(Aeq),b_eq=np.array(beq),bounds=[(0,None)]*4+[(None,None),(0,10)],method='highs')
    if res.success and res.x[5]>1e-9:
        supports.append({'labels':S,'size':len(S),'margin':float(res.x[5]),'positive_normal':[float(x) for x in res.x[:4]]})
# sort concise
supports.sort(key=lambda z:(z['size'],z['labels']))
summary={'status':'NUMERICAL_CANDIDATE_CENSUS','support_count':len(supports),'size_histogram':{},'supports':supports,
 'chart_rule':'Use shared-field relative gaps z_j=exp(-(max h-h_j)); do not make irrational mode-5 exponentials independent.',
 'scope':'positive four-amplitude rays only; not arbitrary signed theta7 limits.'}
for s in supports:summary['size_histogram'][str(s['size'])]=summary['size_histogram'].get(str(s['size']),0)+1
(ROOT/'results/R7P-016_positive_ray_support_census.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({'R7P-009':coordinate,'R7P-010':d12,'R7P-015':domain,'R7P-016_support_count':len(supports),'hist':summary['size_histogram']},indent=2))
