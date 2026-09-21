from __future__ import annotations
from pathlib import Path
import sys,json,time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import k16_surrogate_eval as k16
D=json.load(open(ROOT/'checkpoints/R7N-044_full_cover_k16.json'))
CELLS=D['safe_leaves']
OUT=ROOT/'checkpoints/R7N-046_K16_formula_replay.json'
def run(n=3000):
 d=json.load(open(OUT)) if OUT.exists() else {'task':'R7N-046-K16-formula-replay','total':len(CELLS),'processed':0,'passed':0,'failed':[]}
 st=time.time();start=d['processed'];end=min(len(CELLS),start+n)
 for j in range(start,end):
  r=CELLS[j];ok,k,iv=k16.classify(r['box'])
  if ok:d['passed']+=1
  else:d['failed'].append({'index':j,'box':r['box'],'stored_component':r.get('gradient_component'),'stored_interval':r.get('surrogate_interval')})
 d['processed']=end;d['complete']=end==len(CELLS);d['elapsed_last']=time.time()-st
 OUT.write_text(json.dumps(d,indent=2)+'\n')
 print(json.dumps({k:v for k,v in d.items() if k!='failed'},indent=2)); print('failed_count',len(d['failed']))
if __name__=='__main__': run(int(sys.argv[1]) if len(sys.argv)>1 else 3000)
