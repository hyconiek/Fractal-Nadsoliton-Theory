"""Validate the design documents and generate the machine-readable task plan.

This generates design artifacts only; no research task is executed.
"""
from pathlib import Path
from collections import Counter, deque
import hashlib
import json
import re
import sys

ROOT=Path(__file__).resolve().parents[1]
EXPECTED_AGENTS='3317ae9c465a249532ec1fe1f5d79c443a17be8ba9089a55dc1c0948f6c1bf6c'
DOCS=['FIN_PHYSICS_GAP_ANALYSIS.md','FIN_PHYSICS_NEXT_CAMPAIGN_PLAN.md',
      'FIN_PHYSICS_KILL_TESTS.md','FIN_MATH_TO_PHYSICS_BRIDGE_MAP.md',
      'FIN_PHYSICS_CLAIM_LADDER.md','AGENTS_PHYSICS_PROPOSED_PATCH.md']
FIELDS={'Title':'title','Scientific question':'scientific_question','Why it matters':'why_it_matters',
        'Exact inputs':'exact_inputs','Dependencies':'dependencies_text','Method':'method',
        'Deliverables':'deliverables','Acceptance criterion':'acceptance_criterion',
        'Failure / refutation criterion':'failure_refutation_criterion',
        'What becomes possible if successful':'enables','What must NOT be concluded':'nonconclusions',
        'Estimated computational class':'computational_class','Priority':'priority'}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    assert sha(ROOT/'AGENTS.md')==EXPECTED_AGENTS,'AGENTS.md changed: review before generating'
    plan=(ROOT/DOCS[1]).read_text()
    heads=list(re.finditer(r'^### (PHY-\d{3})\s*$',plan,re.M));tasks=[]
    for n,h in enumerate(heads):
        body=plan[h.end():heads[n+1].start() if n+1<len(heads) else len(plan)]
        fields=dict(re.findall(r'^- \*\*([^*]+):\*\* (.+)$',body,re.M))
        assert set(fields)==set(FIELDS),(h[1],set(fields)^set(FIELDS))
        row={'id':h[1],**{FIELDS[k]:v for k,v in fields.items()}}
        row['computational_class']=row['computational_class'].rstrip('.')
        row['priority']=row['priority'].rstrip('.')
        assert row['computational_class'] in ['S','M','L'],row['id']
        assert row['priority'] in ['P0','P1','P2','exploratory'],row['id']
        packages=re.findall(r'^## Package ([A-Z]) — (.+)$',plan[:h.start()],re.M)
        row['package_id'],row['package_name']=packages[-1]
        row['dependencies']=list(dict.fromkeys(re.findall(r'PHY-\d{3}',row['dependencies_text'])))
        row['dependency_policy']='Outcome review; rejected/skipped branches allowed' if row['id'] in ['PHY-011','PHY-027','PHY-030','PHY-031'] else 'Use proved premises only; negative outcomes may license a documented alternative, never a false theorem'
        row['execution_status']='NOT_STARTED'
        row['claim_status']='RESEARCH_QUESTION_NOT_RESULT'
        tasks.append(row)
    assert [r['id'] for r in tasks]==[f'PHY-{n:03d}' for n in range(1,32)]
    ids={r['id'] for r in tasks};deg={r['id']:len(r['dependencies']) for r in tasks}
    kids={i:[] for i in ids}
    for r in tasks:
        assert set(r['dependencies'])<=ids and r['id'] not in r['dependencies']
        for d in r['dependencies']:kids[d].append(r['id'])
    q=deque(sorted(i for i,d in deg.items() if not d));order=[]
    while q:
        i=q.popleft();order.append(i)
        for j in kids[i]:
            deg[j]-=1
            if not deg[j]:q.append(j)
    assert len(order)==len(tasks),'Dependency cycle'
    links=0
    for name in DOCS:
        s=(ROOT/name).read_text()
        assert len(re.findall(r'^```',s,re.M))%2==0,name
        for target in re.findall(r'\]\(([^)]+)\)',s):
            if target.startswith(('http://','https://','#')):continue
            path=target.split('#')[0]
            assert (ROOT/path).exists(),(name,target)
            links+=1
    kills=re.findall(r'^## (KT-\d{2}) —', (ROOT/DOCS[2]).read_text(),re.M)
    assert kills==[f'KT-{i:02d}' for i in range(1,11)]
    audit=json.loads((ROOT/'fin_physics_review/evidence_audit.json').read_text())
    assert audit['agents_sha256']==EXPECTED_AGENTS
    inputs=dict(re.findall(r'^\| (E\d) \| (.+) \|$',plan,re.M))
    assert set(inputs)=={f'E{i}' for i in range(9)}
    result={
      'date':'2026-09-23','namespace':'PHY','status':'DESIGN_ONLY_NOT_EXECUTED',
      'verdict':'MATHEMATICAL CORE CLOSED EXCEPT FOR NON-BLOCKING ITEMS — MOVE TO PHYSICS',
      'verdict_scope':'Sufficient finite reference interface for conditional bridge research, not a physical-theory claim',
      'agents_unchanged_sha256':EXPECTED_AGENTS,
      'recommended_next_campaign':{'package_id':'R','title':'Relational composition, informational continuity and causal structure',
        'preflight':['PHY-001'],'first_tasks':['PHY-002','PHY-003','PHY-004','PHY-005','PHY-006'],
        'first_decision_gate':'After PHY-006; stop or narrow failed lanes before geometry/field promotion',
        'next_if_pass':['PHY-007','PHY-008'],
        'deferred':['PHY-010','PHY-026','PHY-027']},
      'computational_classes':{'S':{'minutes_per_run':5,'memory_gib':2},
                              'M':{'minutes_per_run':30,'memory_gib':4},
                              'L':{'minutes_total_per_task':120,'minutes_per_subjob':30,'memory_gib':6}},
      'authorization':'Plan only; no campaign launch, AGENTS edits, uploads, installations, purchases or empirical collection',
      'input_catalog':inputs,'tasks':tasks,'topological_order':order,'kill_tests':kills,
      'validation':{'task_count':len(tasks),'required_fields_per_task':13,'acyclic':True,
                    'local_links_checked':links,'agents_unchanged':True},
      'source_documents_sha256':{name:sha(ROOT/name) for name in DOCS},
      'evidence_audit_sha256':sha(ROOT/'fin_physics_review/evidence_audit.json')}
    path=ROOT/'FIN_PHYSICS_PRIORITY_TASKS.json'
    if '--write' in sys.argv:path.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    else:assert json.loads(path.read_text())==result,'Priority JSON is stale'
    print(json.dumps(result['validation'],indent=2))

if __name__=='__main__':main()
