"""Exact multilinear bookkeeping, numerical diagnostics and source hashes."""
from pathlib import Path
from fractions import Fraction as F
import itertools,json,hashlib,math
P=Path(__file__).resolve().parent;checks=[]
def check(name,value):
 checks.append({'name':name,'pass':bool(value)})
 if not value:raise AssertionError(name)
# An arbitrary ordered trilinear tensor tests the telescoping identity without
# assuming symmetry; rational coefficients retain every ordered cross term.
a=[F(3,5),F(2,5)];b=[F(1,2),F(-1,7)]
def tri(x,y,z):return sum(F(1+4*i+2*j+k)*x[i]*y[j]*z[k] for i,j,k in itertools.product(range(2),repeat=3))
e=[x-y for x,y in zip(a,b)]
check('ordered trilinear telescoping',tri(a,a,a)-tri(b,b,b)==tri(e,a,a)+tri(b,e,a)+tri(b,b,e))
check('absolute triple sum mass',sum(abs(a[i]*a[j]*a[k]) for i,j,k in itertools.product(range(2),repeat=3))==sum(map(abs,a))**3)
check('coherent mixed terms are retained',tri(a,a,a)!=sum(a[i]**3*tri([int(j==i) for j in range(2)],[int(j==i) for j in range(2)],[int(j==i) for j in range(2)]) for i in range(2)))
check('critical squared scaling',6-2*2==2)
check('shape carrier squared gap',4*(4+F(5,4))-3==18)
D=json.loads((P/'superposition_preview.json').read_text())['results']
for j,d in enumerate(D):
 check(f'nonzero mixed response {j}',d['mixed_cross_trace_norm']>.5*d['mixture_trace_norm'])
 check(f'cubic scaling {j}',d['cubic_homogeneity_relative_error']<1e-11)
 for i,run in enumerate(d['runs']):
  check(f'work balance {j}/{i}',max(run['relative_energy_work_defects'])<.0006)
  check(f'bound projections {j}/{i}',max(run['max_even_bound_projections'])<1e-10 and run['max_odd_shape_projection']<1e-10)
  check(f'boundary tails {j}/{i}',max(run['final_outer_band_phase_squared'])<1e-15)
for key in ['mixture_trace_norm','mixed_cross_trace_norm','diagonal_only_trace_norm']:
 check('time refinement '+key,abs(D[0][key]/D[1][key]-1)<.0002)
for i in range(4):
 check(f'work convergence {i}',all(3.5<a/b<4.5 for a,b in zip(D[0]['runs'][i]['relative_energy_work_defects'],D[1]['runs'][i]['relative_energy_work_defects'])))
for f in json.loads((P/'INPUTS.json').read_text())['files']:
 raw=(P.parent/f['path']).read_bytes()
 sha=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
 check('source '+f['path'],hashlib.sha256(raw).hexdigest()==f['sha256'] and len(raw)==f['size_bytes'] and sha==f['git_blob_sha'])
(P/'verification.json').write_text(json.dumps({'scope':'Finite diagnostics and provenance, not formal proof verification','passed':len(checks),'failed':0,'checks':checks},indent=2)+'\n')
print(json.dumps({'passed':len(checks),'failed':0}))
