"""Reproducibility diagnostics, not a verifier of the analytical proof."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math
P=Path(__file__).resolve().parent;ROOT=P.parent
checks=[]
def check(name,ok):
 checks.append({'name':name,'pass':bool(ok)})
 if not ok: raise AssertionError(name)
# Exact exponent bookkeeping: T^2 sqrt(T) = R^5.
for name,source_power,bound in [('continuum',F(-1,2),4),('direct',F(-1),3)]:
 check(name+' squared trace exponent',5+2*source_power==bound)
 check(name+' amplitude exponent',6-bound=={'continuum':2,'direct':3}[name])
check('shape exponent',-3+4==1 and 6-1==5)
check('shape squared carrier gap at minimum mass',4*(4+F(5,4))-3==18)
check('highest nonlinearity first contributes at order five',2*2+1==5)
D=json.loads((P/'feedback_preview.json').read_text())['results']
for n,d in enumerate(D):
 G=d['q_gram'];tag=f'run {n}'
 check(tag+' time window',d['T']==d['R']**2)
 check(tag+' Gram symmetry',max(abs(G[i][j]-G[j][i]) for i in range(3) for j in range(3))<1e-12)
 check(tag+' positive Gram diagonal',all(G[i][i]>=0 for i in range(3)))
 check(tag+' Gram two by two minors',all(G[i][i]*G[j][j]-G[i][j]**2>=-1e-10 for i in range(3) for j in range(i)))
 det=G[0][0]*(G[1][1]*G[2][2]-G[1][2]**2)-G[0][1]*(G[0][1]*G[2][2]-G[1][2]*G[0][2])+G[0][2]*(G[0][1]*G[1][2]-G[1][1]*G[0][2])
 check(tag+' Gram determinant',det>=-1e-9)
 check(tag+' recombination',math.isclose(sum(map(sum,G)),d['q_total_g4_one_squared_integral'],rel_tol=1e-12))
 check(tag+' work balance',max(d['relative_energy_work_defects'])<.0003)
 check(tag+' even projection',max(d['max_even_bound_projections'])<1e-10)
 check(tag+' odd projection',d['max_odd_shape_projection']<1e-10)
 check(tag+' outer band',max(d['final_outer_band_phase_squared'])<1e-15)
for a,b in [(0,1),(4,5)]:
 check(f'time refinement {a}',abs(D[a]['q_total_g4_one_squared_integral']/D[b]['q_total_g4_one_squared_integral']-1)<.0004)
 check(f'work refinement {a}',all(3.5<x/y<4.5 for x,y in zip(D[a]['relative_energy_work_defects'],D[b]['relative_energy_work_defects'])))
check('grid quadrature refinement',abs(D[2]['q_total_g4_one_squared_integral']/D[1]['q_total_g4_one_squared_integral']-1)<.0003)
check('domain refinement',abs(D[3]['q_total_g4_one_squared_integral']/D[1]['q_total_g4_one_squared_integral']-1)<1e-10)
for f in json.loads((P/'INPUTS.json').read_text())['files']:
 b=(ROOT/f['path']).read_bytes()
 check('input '+f['path'],hashlib.sha256(b).hexdigest()==f['sha256'] and len(b)==f['size_bytes'])
result={'scope':'Exact scaling identities, finite numerical diagnostics and input provenance; not formal proof verification','passed':len(checks),'failed':0,'checks':checks}
(P/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'passed':len(checks),'failed':0}))
