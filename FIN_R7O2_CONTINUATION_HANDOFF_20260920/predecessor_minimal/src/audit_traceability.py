from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
SCIENTIFIC={'EXACT_PROVED','INTERVAL_CERTIFIED','CONDITIONAL_LEMMA','NUMERICAL_REPRODUCED','NUMERICAL_NEW','COUNTEREXAMPLE_CERTIFIED','COUNTEREXAMPLE_NUMERICAL','UNRESOLVED'}
def exists(p): return (ROOT/p).exists()
def main(require_all_terminal=False):
    tasks=json.load(open(ROOT/'TASKS.json'))['tasks']; claims=json.load(open(ROOT/'CLAIMS.json')); reg=json.load(open(ROOT/'THEOREM_REGISTER.json'))
    missing=[]; orphan=[]; pending=[]
    terminal={'DONE','RESOURCE_STOP','BLOCKED_DEPENDENCY','SUPERSEDED_BY_COUNTEREXAMPLE','OUT_OF_SCOPE'}
    for t in tasks:
        if t['execution_status'] not in terminal: pending.append(t['ID'])
        if t['scientific_status'] in SCIENTIFIC and t['scientific_status']!='NOT_A_SCIENTIFIC_CLAIM':
            if not t.get('evidence') and t['execution_status']=='DONE': orphan.append((t['ID'],'no evidence'))
        for p in t.get('evidence',[]):
            if not exists(p): missing.append((t['ID'],p))
    for name,c in claims.items():
        p=c.get('evidence')
        if p and not exists(p): missing.append((name,p))
        if c.get('status','').startswith('INTERVAL_CERTIFIED') and not p: orphan.append((name,'accepted claim without evidence'))
    for name,r in reg.items():
        for p in r.get('proof_files',[])+([r['certificate_file']] if r.get('certificate_file') else []):
            if not exists(p): missing.append((name,p))
    ok=not missing and not orphan and (not require_all_terminal or not pending)
    out={'task':'R7N-057','all_paths_resolve':not missing,'zero_orphan_accepted_claims':not orphan,'all_tasks_terminal':not pending,'pending_tasks':pending,'missing_paths':missing,'orphan_claims':orphan,'theorem_register_entries':len(reg),'pass':ok}
    (ROOT/'results/R7N-057_traceability.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
    if not ok and require_all_terminal: raise SystemExit(1)
if __name__=='__main__': main(False)
