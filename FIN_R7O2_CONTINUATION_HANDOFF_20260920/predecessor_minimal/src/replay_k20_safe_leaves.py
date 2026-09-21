from __future__ import annotations
from pathlib import Path
import sys,json,time
ROOT=Path(__file__).resolve().parents[1];IR=ROOT/'inputs/intake_review_20260919';sys.path.insert(0,str(IR));import scientific_rechecks as sr
sr.iv.dps=30
sys.path.insert(0,str(ROOT/'src'));import k20_surrogate_eval as k20
D=json.load(open(ROOT/'checkpoints/R7N-044_K20_residual.json'));A=json.load(open(ROOT/'checkpoints/R7N-044_K20_adaptive.json'))
CELLS=[('direct',i,r) for i,r in enumerate(D['safe_leaves'])]+[('adaptive',i,r) for i,r in enumerate(A['safe_leaves'])]
OUT=ROOT/'checkpoints/R7N-046_K20_formula_replay.json'
def run(n=5000):
 d=json.load(open(OUT)) if OUT.exists() else {'task':'R7N-046-K20-formula-replay','total':len(CELLS),'processed':0,'passed':0,'failed':[]}
 st=time.time();start=d['processed'];end=min(len(CELLS),start+n)
 for j in range(start,end):
  layer,i,r=CELLS[j];ok,k,iv=k20.classify(r['box'])
  if ok:d['passed']+=1
  else:d['failed'].append({'global_index':j,'layer':layer,'index':i,'box':r['box']})
 d['processed']=end;d['complete']=end==len(CELLS);d['elapsed_last']=time.time()-st
 OUT.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({k:v for k,v in d.items() if k!='failed'},indent=2));print('failed_count',len(d['failed']))
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 5000)
