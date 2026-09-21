from __future__ import annotations
from pathlib import Path
import json,sys,time,hashlib
from concurrent.futures import ProcessPoolExecutor
ROOT=Path(__file__).resolve().parents[1]
CERT=ROOT/'certificates/active_leaf_certificates.jsonl'
sys.path.insert(0,str(ROOT/'src'))

def one(line):
    import json,sys
    from pathlib import Path
    root=Path(__file__).resolve().parents[1]
    sys.path.insert(0,str(root/'src'))
    from fixed_witness_checker import certify_fixed, exact_match
    c=json.loads(line);o=certify_fixed(c,9);m=exact_match(c,o)
    ok=o['ok'] and all(m.values()) and o['reason']==c['pd_method'].replace('SYLVESTER','PHYSICAL_CENTERED_SYLVESTER_PD').replace('GERSHGORIN','PHYSICAL_CENTERED_GERSHGORIN_PD')
    # producer pd_method names may simply be SYLVESTER/GERSHGORIN
    if c['pd_method']=='SYLVESTER': reason_ok=o['reason']=='PHYSICAL_CENTERED_SYLVESTER_PD'
    elif c['pd_method']=='GERSHGORIN': reason_ok=o['reason']=='PHYSICAL_CENTERED_GERSHGORIN_PD'
    else: reason_ok=False
    ok=o['ok'] and all(m.values()) and reason_ok
    return {'leaf_id':c['leaf_id'],'original_index':c['original_index'],'ok':ok,'reason':o['reason'],'exact_match':m}

def main():
    start=int(sys.argv[1]);end=int(sys.argv[2]);workers=int(sys.argv[3]) if len(sys.argv)>3 else 8
    lines=[]
    with CERT.open() as f:
        for i,line in enumerate(f):
            if i>=end:break
            if i>=start:lines.append(line)
    t=time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex: rows=list(ex.map(one,lines,chunksize=4))
    out={'range':[start,end],'count':len(rows),'pass':sum(r['ok'] for r in rows),'fail':sum(not r['ok'] for r in rows),'seconds':time.time()-t,'failures':[r for r in rows if not r['ok']]}
    outdir=ROOT/'checkpoints/clean_math_shards';outdir.mkdir(parents=True,exist_ok=True)
    p=outdir/f'shard_{start:05d}_{end-1:05d}.json';p.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({**out,'failures':out['failures'][:3]},indent=2))
if __name__=='__main__':main()
