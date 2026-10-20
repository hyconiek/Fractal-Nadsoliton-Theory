#!/usr/bin/env python3
"""Exact finite and interval controls for the paid spectator-limit proof."""
import argparse,ast,hashlib,json,re,subprocess
from fractions import Fraction as F
from itertools import product
from math import factorial
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
INPUTS=[
 'AGENTS.md','research_tasks/00_EXECUTION_RULES.md',
 'fin_spectator_composition_research_20261020/NEXT_LIMIT_PROOF_ORDER.md',
 'fin_spectator_composition_research_20261020/PROOF_OR_OPEN_LEMMA.md',
 'fin_hyp_history_entropy_20261020/MODEL_CONTRACT.json',
 'fin_hyp_history_entropy_20261020/AUTHORIZATION_AND_MODEL.md',
 'fin_hyp_history_entropy_review_20261020/PROOF_AUDIT_AND_DUAL.md',
 'fin_reference_transport_research_20261020/PROOF_OR_OPEN_LEMMA.md',
 'fin_kinetic_invariants_campaign_20261014/KI-002/REVERSIBILITY_AND_MOBILITY.md',
 'fin_rank7_mathphysics_next/proofs/MP7-031_finite_N_variational_model.md',
 'FIN_ST2532_ST2621_Configuration_Generator_and_Dynamic_RG_Report_EN.md',
 'fin_reversible_source_transfer_research_20261020/SEARCH_RESULTS.json']
PROTECTED=['fin_spectator_composition_research_20261020','fin_reference_transport_research_20261020',
 'fin_hyp_history_entropy_20261020','fin_hyp_history_entropy_review_20261020']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def manifest(folder):
 seen=set()
 for line in (folder/'MANIFEST.sha256').read_text().splitlines():
  h,n=line.split('  ',1);p=Path(n)
  assert n not in seen and not p.is_absolute() and '..' not in p.parts
  assert sha(folder/n)==h,(folder,n);seen.add(n)
 actual={p.relative_to(folder).as_posix() for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
 assert actual==seen|{'MANIFEST.sha256'},folder
 return len(seen)
def norm(a):return max(sum(abs(x) for x in row) for row in a)
def exp_minus_box(x):
 assert 0<=x<=1
 lo=sum(((-x)**k/F(factorial(k)) for k in range(26)),F(0))
 return lo,lo+x**26/F(factorial(26))
def edge_flow_root(c,pi):
 a=c*(c+2);lo=F(0);hi=min(pi);f=lambda t:a*t*t+t-pi[0]*pi[1]
 assert f(lo)<0<f(hi)
 for _ in range(64):
  mid=(lo+hi)/2
  if f(mid)<0:lo=mid
  elif f(mid)>0:hi=mid
  else:return mid,mid
 assert f(lo)<0<f(hi)
 return lo,hi
def state_at(t,jumps,values):return values[sum(x<=t for x in jumps)]
def warp(t,knots):
 for (a,b),(c,d) in zip(knots,knots[1:]):
  if a<=t<=c:return b+(t-a)*(d-b)/(c-a)
 raise ValueError(t)
def block_law(mu,p,n):
 out={}
 for w in product(range(len(mu)),repeat=n+1):
  x=mu[w[0]]
  for i,j in zip(w,w[1:]):x*=p[i][j]
  out[w]=x
 return out
def check_markdown():
 for path in HERE.glob('*.md'):
  lines=path.read_text().splitlines();assert all(x==x.rstrip() for x in lines),path
  prose=[];fence=False
  for line in lines:
   if line.lstrip().startswith('```'):fence=not fence
   elif not fence and not line.startswith(('    ','\t')):prose.append(line)
  text=re.sub(r'`[^`\n]*`','', '\n'.join(prose))
  for link in re.findall(r'\[[^\]\n]+\]\(([^)\n]+)\)',text):
   link=link.strip('<>').split('#')[0]
   if link and not link.startswith(('http:','https:')):assert (path.parent/link).exists(),(path,link)
def duplicate_search():
 source=json.loads((ROOT/'fin_reversible_source_transfer_research_20261020/SEARCH_RESULTS.json').read_text())
 names=[]
 for name in source['searched_files']:
  p=Path(name)
  if len(p.parts)==1 and re.search(r'(Process|Generator|Markov|Kinetic|Entropy|Reference|Transfer|Bridge|Clock|Source)',p.name,re.I):names.append(name)
 for d in ('fin_hyp_history_entropy_20261020','fin_hyp_history_entropy_review_20261020','fin_reference_transport_research_20261020','fin_spectator_composition_research_20261020'):
  names.extend(str(p.relative_to(ROOT)) for p in (ROOT/d).glob('*.md'))
 names=sorted(set(names));hits=[]
 pat=r'square.root.{0,55}(activity|rate|limit)|uniformization|Skorokhod|jump.times.{0,30}grid|diagonal.{0,45}(limit|completion)|spectator.{0,50}limit'
 for i in range(0,len(names),90):
  r=subprocess.run(['rg','-n','-i','--max-columns','330','--',pat,*names[i:i+90]],cwd=ROOT,text=True,capture_output=True,timeout=60)
  assert r.returncode in (0,1);hits.extend(r.stdout.splitlines())
 return {'selected_files':len(names),'pattern':pat,'matching_lines':len(hits),'lines':hits,
  'scope':'bounded duplicate search, not global novelty or source absence',
  'reused':'square-root-rate shape and finite Markov methods; new paid spectator completion, bounds and topology'}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--record',action='store_true');ap.add_argument('--check',action='store_true');ap.add_argument('--output',type=Path)
 args=ap.parse_args();checks=[]
 def ck(name,value):
  if not value:raise AssertionError(name)
  assert name not in checks;checks.append(name)
 inputs={n:sha(ROOT/n) for n in INPUTS}
 old={n:{'entries':manifest(ROOT/n),'manifest_sha256':sha(ROOT/n/'MANIFEST.sha256')} for n in PROTECTED}
 for c in (F(2),F(6),F(10)):
  q=3;a=F(q)
  p=[[((c+1) if i==j else F(1))/(c+q) for j in range(q)] for i in range(q)]
  qc=[[c*(p[i][j]-F(i==j)) for j in range(q)] for i in range(q)]
  limit=[[F(-2) if i==j else F(1) for j in range(q)] for i in range(q)]
  ck(f'uniform_source_rows_{c}',all(sum(row)==1 and min(row)>=0 for row in p))
  ck(f'uniform_generator_diagonal_{c}',all(sum(row)==0 for row in qc))
  ck(f'uniform_generator_error_{c}',norm([[qc[i][j]-limit[i][j] for j in range(q)] for i in range(q)])<=4*a*a/c)
  ck(f'uniform_scaling_squared_bound_{c}',c/(c+q)<=1 and (1-a/c<=0 or c/(c+q)>=(1-a/c)**2))
 # Fixed nonuniform pi; enclose a real nonlinear source root at each named c.
 pi=[F(1,5),F(4,5)];a=F(3);certificates=[]
 for c in (F(4),F(16),F(64)):
  lo,hi=edge_flow_root(c,pi)
  ck(f'nonuniform_root_{c}',0<lo<=hi<min(pi) and hi-lo<F(1,2**60))
  for x in range(2):
   rlo=c*(pi[x]-hi)/((c+1)*pi[x]);rhi=c*(pi[x]-lo)/((c+1)*pi[x])
   ck(f'nonuniform_scaling_bound_{c}_{x}',rhi<=1 and rlo>=(1-a/c)**2)
  target=[F(2),F(1,2)]
  ck(f'nonuniform_generator_bound_{c}',max(2*(target[x]-c*lo/pi[x]) for x in range(2))<=4*a*a/c)
  elo,ehi=exp_minus_box(F(5,2)/c);errors=[]
  for x in range(2):
   y=1-x;source=(lo/pi[x],hi/pi[x]);ref=(pi[y]*(1-ehi),pi[y]*(1-elo))
   errors.append(max(abs(u-v) for u in source for v in ref))
  ck(f'nonuniform_one_step_grid_TV_bound_{c}',max(errors)<=3*a*a/(c*c))
  certificates.append({'c':str(c),'flow_lower':str(lo),'flow_upper':str(hi),'pi':['1/5','4/5']})
 ck('limit_reversibility_nonuniform',pi[0]*2==pi[1]*F(1,2))
 ck('nonuniform_limit_activity_fixed',F(2)*F(1,2)==1)
 grid=[F(1,3),F(2,3)];actual=[F(1,4),F(3,5)];knots=[(F(0),F(0)),*zip(grid,actual),(F(1),F(1))]
 ck('warp_order_and_endpoint',all(knots[i][0]<knots[i+1][0] and knots[i][1]<knots[i+1][1] for i in range(3)))
 ck('warp_displacement',max(abs(x-y) for x,y in knots)<=F(1,3))
 times=sorted(set([x for x,y in knots]+[(knots[i][0]+knots[i+1][0])/2 for i in range(3)]))
 ck('warp_state_alignment',all(state_at(t,grid,[0,1,2])==state_at(warp(t,knots),actual,[0,1,2]) for t in times))
 e1lo,e1hi=exp_minus_box(F(1))
 for c in (2,6,10):ck(f'grid_nojump_common_mass_{c}',F(c+1,c+2)**c>e1hi)
 ck('full_path_TV_above_3_5',1-e1hi>F(3,5))
 shared=[[F(2,3),F(1,3)],[F(1,3),F(2,3)]]
 p0=block_law([F(1),F(0)],shared,3);q0=block_law([F(1,2)]*2,shared,3)
 ck('initial_TV_paid_once',sum(abs(p0[w]-q0[w]) for w in p0)/2==F(1,2))
 b=[1,2,3];kp=[[F(x*y) for y in b] for x in b];dp=[F(x*x) for x in b]
 ck('completion_normalized_reference_invariant',all(kp[i][j]/F(b[i]*b[j])==1 for i,j in product(range(3),repeat=2)))
 for c in (2,6):
  pc=[[F(c+1,c+3) if i==j else F(1,c+3) for j in range(3)] for i in range(3)]
  transformed=[[F(3)*(kp[i][j]+(c*dp[i] if i==j else 0))*F(1,3*(c+3)*b[i]*b[j]) for j in range(3)] for i in range(3)]
  ck(f'joint_gauge_finite_source_{c}',transformed==pc)
 ck('changing_only_K_changes_relative_activity',kp[0][1]/kp[0][2]==F(2,3)!=1)
 weights=[1,2,3];pi3=[F(x*x,14) for x in weights]
 q3=[[F(weights[j],weights[i]) if i!=j else F(0) for j in range(3)] for i in range(3)]
 for i in range(3):q3[i][i]=-sum(q3[i])
 ck('derived_three_state_conservative',all(sum(row)==0 for row in q3))
 ck('derived_three_state_balance',all(pi3[i]*q3[i][j]==pi3[j]*q3[j][i] for i,j in product(range(3),repeat=2)))
 jl=q3[2][0]/(-q3[2][2]);jr=pi3[0]/(1-pi3[2])
 ck('retained_clock_free_reset_mismatch',jl==F(1,3) and jr==F(1,5) and jl-jr==F(2,15))
 # Labelled q3/M2 -> counts with departure and reverse multiplicities.
 states=list(product(range(3),repeat=2));s=[1,2,2];w=[s[(y-x)%3] for x,y in states];pl=[F(x*x,27) for x in w]
 ql=[[F(w[j],w[i]) if i!=j and sum(a!=b for a,b in zip(states[i],states[j]))==1 else F(0) for j in range(9)] for i in range(9)]
 for i in range(9):ql[i][i]=-sum(ql[i])
 count=lambda x:tuple(x.count(a) for a in range(3));types=sorted({count(x) for x in states});groups=[[i for i,x in enumerate(states) if count(x)==n] for n in types]
 pm=[sum(pl[i] for i in g) for g in groups];qm=[]
 for g in groups:
  i=g[0];row=[sum(ql[i][j] for j in h) for h in groups];qm.append(row)
  ck(f'actual_count_lumping_{i}',all([sum(ql[x][j] for j in h) for h in groups]==row for x in g))
 valid=True
 for i,n in enumerate(types):
  for j,np in enumerate(types):
   if i!=j and qm[i][j]:
    oldlabel=[a for a in range(3) if np[a]-n[a]==-1][0];newlabel=[a for a in range(3) if np[a]-n[a]==1][0]
    valid=valid and qm[i][j]==F(n[oldlabel]*w[groups[j][0]],w[groups[i][0]])
    valid=valid and qm[i][j]*qm[j][i]==n[oldlabel]*(n[newlabel]+1)
 ck('count_activity_and_cavity_ratio',valid)
 ck('count_exact_invariant_measure',sum(pm)==1 and all(sum(pm[i]*qm[i][j] for i in range(6))==0 for j in range(6)))
 ck('old_sources_unchanged',inputs=={n:sha(ROOT/n) for n in INPUTS})
 ck('old_packages_unchanged',old=={n:{'entries':manifest(ROOT/n),'manifest_sha256':sha(ROOT/n/'MANIFEST.sha256')} for n in PROTECTED})
 result={'date':'2026-10-20','source_model_id':'FIN-HYP-HIST-ME-001','derived_limit_id':'FIN-HYP-HIST-ME-001-SL-SQRT',
  'status':'PASS_BOUNDED_SCALING_GENERATOR_GRID_AND_ACTIVITY_CONTROLS','bounded_checks':len(checks),'checks':checks,
  'nonuniform_root_certificates':certificates,'arithmetic':'exact Fraction with outward rational root/exp bounds; no binary64 theorem certificate',
  'general_path_theorem':'analytical J1 time-warp coupling, not a simulation proof','full_path_TV_countercontrol':'two-state T=1 exact 1-exp(-1), bounded below by 3/5',
  'grid_row_TV_constant':'3 A^2/c^2','J1_bad_probability_constant':'((7/2) A^2 T+2 A)/c plus initial TV',
  'new_physical_postulate_or_clock':False,'primary_completion':'D=I unchanged',
  'conditional_activity_provider':'DERIVED_FROM_SUPPLIED_HYPOTHESIS_AND_COMPLETION','strict_F01':'OPEN_SOURCE',
  'native_HB_adopted_or_repaired':False,'source_inputs':inputs,'protected_packages':old,'independent_peer_review':False,
  'goal_objective':'badaj teorię','goal_finished':False}
 if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
 if args.record:
  assert args.output and args.output.resolve()==HERE/'RESULTS.json'
  (HERE/'DUPLICATE_SEARCH.json').write_text(json.dumps(duplicate_search(),indent=2)+'\n')
  (HERE/'INPUTS.sha256').write_text(''.join(f'{h}  {n}\n' for n,h in inputs.items()))
  (HERE/'MANIFEST.sha256').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='MANIFEST.sha256'))
 if args.check:
  manifest(HERE)
  assert json.loads((HERE/'RESULTS.json').read_text())==result
  assert (HERE/'INPUTS.sha256').read_text()==''.join(f'{h}  {n}\n' for n,h in inputs.items())
  for p in HERE.glob('*.json'):json.loads(p.read_text())
  for p in HERE.glob('*.py'):ast.parse(p.read_text())
  check_markdown()
 print(json.dumps({'status':result['status'],'bounded_checks':len(checks),'conditional_limit':'DERIVED_IN_PROOF','strict_F01':'OPEN_SOURCE'}))
if __name__=='__main__':main()
