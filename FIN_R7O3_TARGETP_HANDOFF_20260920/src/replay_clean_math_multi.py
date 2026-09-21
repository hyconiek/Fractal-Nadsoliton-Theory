from __future__ import annotations
from pathlib import Path
import json,sys,time
from concurrent.futures import ProcessPoolExecutor
ROOT=Path(__file__).resolve().parents[1]; CERT=ROOT/'certificates/active_leaf_certificates.jsonl'
sys.path.insert(0,str(ROOT/'src'))
from replay_clean_math import one

def main():
    start=int(sys.argv[1]); stop=int(sys.argv[2]); step=int(sys.argv[3]); workers=int(sys.argv[4])
    lines=CERT.read_text().splitlines()
    outdir=ROOT/'checkpoints/clean_math_shards';outdir.mkdir(parents=True,exist_ok=True)
    with ProcessPoolExecutor(max_workers=workers) as ex:
      for a in range(start,stop,step):
        b=min(a+step,stop); p=outdir/f'shard_{a:05d}_{b-1:05d}.json'
        if p.exists():
          print('SKIP',a,b,flush=True);continue
        t=time.time(); rows=list(ex.map(one,lines[a:b],chunksize=3))
        out={'range':[a,b],'count':len(rows),'pass':sum(r['ok'] for r in rows),'fail':sum(not r['ok'] for r in rows),'seconds':time.time()-t,'failures':[r for r in rows if not r['ok']]}
        p.write_text(json.dumps(out,indent=2)+'\n'); print('DONE',a,b,out['pass'],out['fail'],round(out['seconds'],2),flush=True)
if __name__=='__main__':main()
