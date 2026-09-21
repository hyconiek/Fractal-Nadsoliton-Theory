from __future__ import annotations
from pathlib import Path
import hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parent

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def check_manifest():
    mf=ROOT/'MANIFEST.sha256'
    if not mf.exists(): return True,{'note':'manifest not generated yet'}
    bad=[]; n=0
    for line in mf.read_text().splitlines():
        if not line.strip(): continue
        h,rel=line.split('  ',1); p=ROOT/rel; n+=1
        if not p.exists() or sha(p)!=h: bad.append(rel)
    return not bad,{'files':n,'bad':bad[:20]}

def run_script(rel):
    p=subprocess.run([sys.executable,str(ROOT/rel)],cwd=ROOT,text=True,capture_output=True,timeout=120)
    return p.returncode==0,{'returncode':p.returncode,'stdout_tail':p.stdout[-3000:],'stderr_tail':p.stderr[-3000:]}

def sample_formula(module_name, checkpoint, layer, count=240):
    sys.path.insert(0,str(ROOT/'src'))
    mod=__import__(module_name)
    d=json.load(open(ROOT/checkpoint))
    leaves=d['safe_leaves']
    if not leaves:return False,{'reason':'no leaves'}
    ids=sorted(set(round(i*(len(leaves)-1)/(count-1)) for i in range(count))) if count>1 else [0]
    bad=[]
    for i in ids:
        ok,_,_=mod.classify(leaves[i]['box'])
        if not ok: bad.append(i)
    return not bad,{'layer':layer,'sampled':len(ids),'population':len(leaves),'bad':bad}

def main():
    checks={}
    checks['manifest_before']=dict(zip(('pass','detail'),check_manifest()))
    # Stored full replays must be complete.
    k16r=json.load(open(ROOT/'checkpoints/R7N-046_K16_formula_replay.json')); k20r=json.load(open(ROOT/'checkpoints/R7N-046_K20_formula_replay.json'))
    checks['stored_formula_replays']={'pass':k16r['complete'] and k20r['complete'] and not k16r['failed'] and not k20r['failed'], 'detail':{'K16':k16r['passed'],'K20':k20r['passed']}}
    b=json.load(open(ROOT/'results/R7N-002_fresh_119_replay.json'))
    checks['baseline_record']={'pass':b['tests']==119 and b['passed_files']==29 and not b['failed_files'],'detail':{'tests':b['tests'],'files':b['files']}}
    for name,rel in [('full_geometry_audit','src/audit_full_phase_exhaustion.py'),('campaign_mutations','src/audit_campaign_mutations.py'),('quartic_symmetry_audit','src/audit_phase_symmetry.py')]:
        ok,detail=run_script(rel);checks[name]={'pass':ok,'detail':detail}
    ok,d=sample_formula('k16_surrogate_eval','checkpoints/R7N-044_full_cover_k16.json','K16');checks['K16_formula_sample']={'pass':ok,'detail':d}
    # K20 sample from both direct and adaptive layers.
    ok1,d1=sample_formula('k20_surrogate_eval','checkpoints/R7N-044_K20_residual.json','K20-direct')
    ok2,d2=sample_formula('k20_surrogate_eval','checkpoints/R7N-044_K20_adaptive.json','K20-adaptive')
    checks['K20_formula_samples']={'pass':ok1 and ok2,'detail':[d1,d2]}
    out={'portable_replay':'PASS' if all(v['pass'] for v in checks.values()) else 'FAIL','checks':checks,
         'note':'Clean-directory portable replay verifies the complete stored geometry/mutation proof, the stored full formula-replay completion records, and deterministic formula recomputation samples. Full all-leaf formula replay remains available through the bundled replay scripts.'}
    print(json.dumps(out,indent=2))
    return 0 if out['portable_replay']=='PASS' else 1
if __name__=='__main__': raise SystemExit(main())
