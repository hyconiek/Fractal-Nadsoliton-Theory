from __future__ import annotations
from pathlib import Path
import sys,json,time
ROOT=Path('/mnt/data/fin_rank7_next_campaign');sys.path.insert(0,str(ROOT/'src'));import k20_surrogate_eval as k20
SRC=ROOT/'checkpoints/R7N-044_full_cover_k16.json'; OUT=ROOT/'checkpoints/R7N-044_K20_residual.json'

def load():
    src=json.load(open(SRC)); cells=src['unresolved_leaves']
    if OUT.exists():
        d=json.load(open(OUT)); return cells,d
    d={'task':'R7N-044-K20-residual','source_checkpoint':str(SRC),'source_processed_total':src['processed_total'],
       'source_unresolved_count':len(cells),'epsilon_k20':str(k20.EPS),'processed_count':0,'safe_count':0,'failed_count':0,
       'safe_leaves':[],'failed_leaves':[],'complete':False}
    return cells,d

def run(n=1000):
    cells,d=load(); start=d['processed_count']; end=min(len(cells),start+n); st=time.time(); safe=d['safe_leaves']; failed=d['failed_leaves']
    for i in range(start,end):
        rec=cells[i]; ok,k,iv=k20.classify(rec['box'])
        row={'source_index':i,'path':rec['path'],'box':rec['box']}
        if ok:
            row.update({'reason':'FULL_SAFE_BY_K20','gradient_component':k,'surrogate_interval':iv}); safe.append(row)
        else:
            failed.append(row)
    d.update({'processed_count':end,'safe_count':len(safe),'failed_count':len(failed),'safe_leaves':safe,'failed_leaves':failed,
              'complete':end==len(cells),'elapsed_seconds_last_chunk':time.time()-st})
    OUT.write_text(json.dumps(d,indent=2)+'\n')
    print(json.dumps({k:v for k,v in d.items() if k not in ('safe_leaves','failed_leaves')},indent=2))
if __name__=='__main__': run(int(sys.argv[1]) if len(sys.argv)>1 else 1000)
