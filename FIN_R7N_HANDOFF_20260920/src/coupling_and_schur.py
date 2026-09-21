from pathlib import Path
import sys,json,math
import numpy as np
ROOT=Path(__file__).resolve().parents[1];H=ROOT/'inputs/FR223_20260916';sys.path.insert(0,str(H/'src'))
import off_face
from compact_model import p_formula
rng=np.random.default_rng(11011)
# R7N-011 direct physical coupling diagnostics
errs=[]
for _ in range(128):
 r=math.exp(rng.uniform(math.log(1/900),0));s=math.exp(rng.uniform(math.log(1/128),0));t=math.exp(rng.uniform(math.log(1/9),0));y=math.exp(rng.uniform(math.log(1e-6),0))
 A=1+2*s*t**3+r*t**4+2*r*s*t
 B=2*math.sqrt(r)*(s*t**(2+math.sqrt(3))+t**2+s*t**(2-math.sqrt(3)))
 q=A/(A+y*B);e=y*B/(A+y*B);p=p_formula(r,s,t,y)
 errs.append(max(abs(q-p[::2].sum()),abs(e-p[1::2].sum()),abs(q+e-1)))
coupling={'task':'R7N-011','samples':128,'max_abs_direct_probability_discrepancy':max(errs),
 'exact_structure':['A_even>0 on the closed hull because anchor term 1 is present','B_odd>=0','q=A/(A+yB)','e=yB/(A+yB)=1-q','for A,B>0, de/dy=A*B/(A+yB)^2>0, so physical e is monotone in y at fixed r,s,t'],
 'warning':'Freeing q/e independently is only an outer relaxation; violations there are not automatically physical counterexamples.',
 'scientific_status':'CONDITIONAL_LEMMA_PLUS_NUMERICAL_CROSSCHECK'}
json.dump(coupling,open(ROOT/'results/R7N-011_physical_coupling.json','w'),indent=2)
# R7N-012 numerical Schur identity cross-check at two thresholds
rows=[]; thresholds=[('sigma',off_face.constants()[1]),('tau0',67/250)]
for _ in range(64):
 x=np.exp(np.log(np.array([1/900,1/128,1/9,1e-6]))+rng.random(4)*(0-np.log(np.array([1/900,1/128,1/9,1e-6]))))
 M,_=off_face.covariance_direct_compact(x)
 rec={'x':x.tolist(),'thresholds':{}}
 for name,a in thresholds:
  K=a*np.eye(4)-M;d=K[3,3];S=K[:3,:3]-np.outer(K[:3,3],K[3,:3])/d
  ik=int(np.sum(np.linalg.eigvalsh(K)<-1e-10));is_=int(np.sum(np.linalg.eigvalsh(S)<-1e-10))
  rec['thresholds'][name]={'pivot':float(d),'neg4':ik,'neg3':is_,'inertia_match':ik==is_}
 rows.append(rec)
assert all(z['thresholds'][n]['pivot']>0 and z['thresholds'][n]['inertia_match'] for z in rows for n,_ in thresholds)
# Demonstrate that transformed matrices differ when threshold changes; maximize over sampled points.
diff=0.0; diff_x=None
for rr in rows:
 M,_=off_face.covariance_direct_compact(rr['x']); mats=[]
 for name,a in thresholds:
  K=a*np.eye(4)-M;d=K[3,3];S=K[:3,:3]-np.outer(K[:3,3],K[3,:3])/d; mats.append((name,a*np.eye(3)-S))
 dd=float(np.max(np.abs(mats[0][1]-mats[1][1])))
 if dd>diff: diff=dd;diff_x=rr['x']
schur={'task':'R7N-012','samples':64,'thresholds':{k:float(v) for k,v in thresholds},'all_sampled_pivots_positive':True,'all_sampled_inertia_matches':True,
       'max_abs_difference_Mtilde_sigma_vs_tau0_over_samples':diff,'difference_witness_x':diff_x,
       'implementation':'src/generic_threshold_shifted.py rebuilds c=lambda6/(3*a) and K=aI-Mtilde for every explicit a; a=None reproduces historical sigma checker exactly on FR26.',
       'scientific_status':'NUMERICAL_CROSSCHECK_PLUS_INTERVAL_IMPLEMENTATION','note':'Exact Schur congruence is algebraic; sampled inertia checks are regression diagnostics, not the proof of the identity.'}
json.dump(schur,open(ROOT/'results/R7N-012_schur_identity_diagnostics.json','w'),indent=2)
print(json.dumps({'R7N-011':coupling,'R7N-012':schur},indent=2))
