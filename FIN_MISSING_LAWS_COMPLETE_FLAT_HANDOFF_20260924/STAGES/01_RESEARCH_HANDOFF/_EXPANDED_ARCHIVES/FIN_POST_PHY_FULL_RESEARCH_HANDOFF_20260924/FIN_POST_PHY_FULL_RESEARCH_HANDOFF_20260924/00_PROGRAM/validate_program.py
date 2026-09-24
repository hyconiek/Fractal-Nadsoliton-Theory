"""Structural/provenance validation only; not a theorem verifier."""
from pathlib import Path
import hashlib,json,re,subprocess,sys
ROOT=Path(__file__).resolve().parents[1];HERE=Path(__file__).resolve().parent
DOCS=['FIN_POST_PHY_RESEARCH_PROGRAM.md','FIN_POST_PHY_TASK_CATALOG.md',
      'FIN_POST_PHY_DEDUP_MATRIX.md','FIN_POST_PHY_KILL_SOURCE_GATES.md']
JSONS=['FIN_POST_PHY_TASKS.json','FIN_POST_PHY_DEPENDENCY_DAG.json']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    data=json.loads((ROOT/JSONS[0]).read_text());dag=json.loads((ROOT/JSONS[1]).read_text())
    tasks=data['tasks'];assert len(tasks)==42 and len({t['id'] for t in tasks})==42
    assert sha(ROOT/'AGENTS.md')==data['agents_sha256']
    required=['title','central_question','importance','input_status','repo_dedup','claim_type',
      'mathematical_object','method','acceptance','refutation','held_out','stop_rule',
      'computational_class','artifacts','allowed_conclusion','forbidden_conclusion',
      'dependencies','score_dimensions','pretarget_data','activation']
    for t in tasks:
        assert all(k in t and (t[k] or k=='dependencies') for k in required),t['id']
        assert t['status']=='PLANNED_NOT_EXECUTED'
        assert t['priority'] in ['P0','P1','P2','P3']
    positions={x:i for i,x in enumerate(dag['topological_order'])}
    assert len(positions)==42
    assert all(positions[e['from']]<positions[e['to']] for e in dag['edges'])
    catalogue=(ROOT/DOCS[1]).read_text()
    assert re.findall(r'^## ([A-Z]+-\d{3}) —',catalogue,re.M)==[t['id'] for t in tasks]
    report=(ROOT/DOCS[0]).read_text()
    assert re.findall(r'^## ([A-K])\. ',report,re.M)==list('ABCDEFGHIJK')
    assert all(i in report for i in data['first_batch']) and 5<=len(data['first_batch'])<=10
    links=0
    for name in DOCS:
        txt=(ROOT/name).read_text()
        assert len(re.findall(r'^```',txt,re.M))%2==0,name
        assert not any(line.rstrip()!=line for line in txt.splitlines()),name
        for target in re.findall(r'\]\(([^)]+)\)',txt):
            if target.startswith(('http:','https:','#')):continue
            path=target.split('#')[0]
            assert (ROOT/path).exists(),(name,path)
            links+=1
    for ref,row in json.loads((HERE/'selected_sources.json').read_text()).items():
        assert sha(ROOT/row['path'])==row['sha256'],ref
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    changed=subprocess.check_output(['git','diff','--name-only',data['reviewed_head']+'..HEAD'],cwd=ROOT,text=True).splitlines()
    # Intervening commits may preserve partial planning work; scientific inputs must stay distinct.
    nonplan=[p for p in changed if not p.startswith(('fin_post_phy_review/','FIN_POST_PHY_'))]
    assert not nonplan,('Scientific baseline changed; re-audit',nonplan)
    out={'status':'PASS_DESIGN_STRUCTURE_AND_SOURCE_IDENTITY','date':'2026-09-24',
         'reviewed_scientific_head':data['reviewed_head'],'head_at_validation':head,
         'intervening_commit_files':changed,'scientific_input_changes_in_intervening_commits':nonplan,
         'task_count':42,'substantive_atoms':37,'audit_synthesis_nodes':3,'dormant_gates':2,
         'first_batch':data['first_batch'],'acyclic':True,'local_links_checked':links,
         'agents_unchanged':True,'agents_sha256':sha(ROOT/'AGENTS.md'),
         'research_executed_by_this_design':False,
         'scope':'Checks document fields, DAG, references and frozen inputs, not the future theorem targets.',
         'output_sha256':{n:sha(ROOT/n) for n in DOCS+JSONS},
         'design_code_sha256':{p.name:sha(p) for p in HERE.glob('*.py')}}
    (HERE/'validation.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['output_sha256','design_code_sha256']},indent=2))
if __name__=='__main__':main()
