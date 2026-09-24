"""Render and validate the post-PHY DESIGN. Never executes research tasks."""
from pathlib import Path
from collections import Counter,deque
import hashlib,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from fin_post_phy_review.program_spec import TASKS,HEAD,AGENTS_HASH,DATE
ROOT=Path(__file__).resolve().parents[1]; HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
TYPES={'THEOREM_TARGET','NUMERICAL_CERTIFICATE_TARGET','CONDITIONAL_MODEL','NO_GO_TARGET',
       'IDENTIFIABILITY_TARGET','OPERATIONAL_PROTOCOL','EMPIRICAL_GATE','SPECULATIVE_ONLY'}
DIMENSIONS=['information_gain','branch_falsification','foundational_importance','assumption_independence',
            'computational_cost_penalty','duplication_risk_penalty','heldout_transfer','premise_reduction']
FIRST=['AUD-001','SRC-001','REF-002','REF-004','MEM-001','MEM-003','OP-001','OP-002']

def main():
    assert sha(ROOT/'AGENTS.md')==AGENTS_HASH,'Guardrails changed; inspect before rebuilding'
    scan=json.loads((HERE/'search_index.json').read_text());assert scan['head']==HEAD
    logical=scan['logical_sources']
    campaigns=sorted({k.split('::')[0] for k in logical})
    aliases={'C'+c[:2]:c for c in campaigns}
    def reference(s):
        if len(s)>4 and s[:1]=='C' and s[1:3].isdigit() and s[3]==':':
            key=aliases[s[:3]]+'::'+s[4:]
            row=logical[key];p=ROOT/row['path'];assert sha(p)==row['sha256'],key
            return f'[{s}]({row["path"]})'
        assert (ROOT/s).exists(),s
        return f'[{s}]({s})'
    tasks=json.loads(json.dumps(TASKS));ids={t['id'] for t in tasks};assert len(ids)==len(tasks)==42
    for t in tasks:
        if t['id']=='SYN-001':
            t['dependencies']=[r['id'] for r in tasks if not r['id'].startswith('SYN-')]
        assert t['claim_type'] in TYPES and t['computational_class'] in ['S','M','L','X']
        assert len(t['scores'])==8 and all(0<=x<=5 for x in t['scores'])
        a,b,c,d,e,f,g,h=t['scores'];t['priority_score']=2*a+2*b+2*c+d+g+h-e-f
        score=t['priority_score'];t['priority']='P0' if score>=40 else 'P1' if score>=33 else 'P2' if score>=26 else 'P3'
        t['score_dimensions']=dict(zip(DIMENSIONS,t['scores']))
        assert set(t['dependencies'])<=ids and t['id'] not in t['dependencies']
        for s in t['repo_dedup']['prior_work']:reference(s)
        t['dependency_mode']='outcome_records_including_skips' if t['id'].startswith('SYN-') else 'proved_premises_or_explicitly_scoped_negative_outcomes'
        t['activation']='PROPOSED_ONLY; authorize execution separately'
        if t['id'].startswith('GATE-'):t['activation']='DORMANT; require a genuinely new typed object; otherwise STOP—DO NOT PROMOTE'
        if t['id'].startswith('EMP-'):t['activation']='External data/coordination not authorized; readiness may be assessed symbolically'
        if t['computational_class']=='X':t['activation']+='; separate bounded campaign/budget approval required'
    by={t['id']:t for t in tasks};deg={t['id']:len(t['dependencies']) for t in tasks};children={i:[] for i in ids}
    for t in tasks:
        for d in t['dependencies']:
            assert by[d]['wave']<=t['wave'],(d,t['id'])
            children[d].append(t['id'])
    q=deque(sorted(i for i,d in deg.items() if not d));order=[]
    while q:
        i=q.popleft();order.append(i)
        for j in children[i]:
            deg[j]-=1
            if not deg[j]:q.append(j)
    assert len(order)==42,'Dependency cycle'
    seen=set()
    for i in FIRST:
        assert set(by[i]['dependencies'])<=seen,(i,by[i]['dependencies'])
        seen.add(i)
    ranked=sorted([t for t in tasks if not t['id'].startswith(('AUD-','SYN-','GATE-','EMP-'))],
                  key=lambda t:(-t['priority_score'],-t['scores'][0],t['id']))
    result={'date':DATE,'reviewed_head':HEAD,'head_scope':'Local repository snapshot; historical handoff audit HEAD was 9ad5cd5, not the current HEAD.',
            'status':'RESEARCH_DESIGN_ONLY','agents_sha256':AGENTS_HASH,'task_count':42,
            'task_count_scope':'37 substantive research/protocol atoms, 3 audit/synthesis nodes, 2 dormant promotion gates',
            'priority_rule':'2*information+2*falsification+2*foundation+independence+transfer+premise_reduction-cost-duplication',
            'score_status':'Explicit ordinal scientific judgment, not measured bits of expected information or success probabilities.',
            'priority_thresholds':{'P0':40,'P1':33,'P2':26,'P3':'below 26'},
            'first_batch':FIRST,'ranked_decisive_top10':[t['id'] for t in ranked[:10]],
            'topological_order':order,'tasks':tasks}
    (ROOT/'FIN_POST_PHY_TASKS.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    dag={'reviewed_head':HEAD,'acyclic':True,'nodes':[{'id':t['id'],'wave':t['wave'],'priority':t['priority'],'activation':t['activation']} for t in tasks],
         'edges':[{'from':d,'to':t['id'],'mode':t['dependency_mode']} for t in tasks for d in t['dependencies']],
         'topological_order':order,'first_batch':FIRST}
    (ROOT/'FIN_POST_PHY_DEPENDENCY_DAG.json').write_text(json.dumps(dag,indent=2,ensure_ascii=False)+'\n')
    intro='''# FIN post-PHY task catalogue

Date: %s. Reviewed local HEAD: `%s`.

**Design only.** Every task is PLANNED_NOT_EXECUTED. See
[programme diagnosis](FIN_POST_PHY_RESEARCH_PROGRAM.md),
[de-dup matrix](FIN_POST_PHY_DEDUP_MATRIX.md),
[kill/source gates](FIN_POST_PHY_KILL_SOURCE_GATES.md) and
[machine-readable tasks](FIN_POST_PHY_TASKS.json).

The catalogue has 37 substantive research/protocol atoms, three audit/synthesis
nodes and two dormant promotion gates. It does not manufacture 42 new theorems.
New means a proposed extension relative to the searched FIN snapshot, not a
world-priority claim. Literature methods are benchmarks with hypotheses to pay.

## Common execution contract

Write outputs under `fin_post_phy_campaign/<ID>/`; never edit predecessor
payloads. Before research record exact source hashes, mathematical domains,
PM assumptions, observation/reset laws and numerical acceptance inequalities.
Keep PROVEN, STRONG_NUMERICAL_EVIDENCE, CONDITIONAL, HYPOTHESIS and SPECULATION
separate from execution status. A negative theorem is a successful outcome.

S: normally <=5 min/run and 2 GiB. M: <=30 min/run and 4 GiB.
L: <=2 h/task in <=30-minute checkpointed jobs, <=6 GiB.
X: separate explicit campaign/budget approval; do not launch automatically.
Resource ceilings are not proof-time promises. Use one heavy worker by default.
Stop after two failed attempts with the same cause without a new mathematical
input. A timeout or failed sufficient bound is not a refutation.

Source tasks must state the pretarget information and freeze the law/grammar
before inspecting held-out targets. Do not pick a law because it returns rank
seven, dimension three or a desired constant. A missing independently motivated
source object is a valid STOP outcome, not permission for unlimited inverse design.

Each task returns REPORT.md, a theorem/counterexample note, results.json, replay.py
and a manifest. An analytic-only replay may check exact symbolic identities;
do not invent a numerical certificate for a conceptual proposition. Protocol
tasks additionally return preregistration and raw-record schemas; they collect
no empirical data without authorization. Numeric targets require outward bounds,
not tolerances alone. Never use a local likelihood maximum as an upper bound.

Priority scores are ordinal judgments on eight explicit dimensions. They rank
scientific value, not immediate executability. Dependencies and dormant/external
gates override scores. SYN-001 consumes outcomes, including honest skips, not
an assumption that every upstream hypothesis was proved.

'''%(DATE,HEAD)
    sections=[intro]
    labels=[('central_question','Central question'),('importance','Why it matters'),('input_status','Input status'),
            ('claim_type','Claim type'),('mathematical_object','Exact mathematical object'),('method','Exact method'),
            ('acceptance','Acceptance criterion'),('refutation','Failure/refutation criterion'),
            ('pretarget_data','Pretarget information'),('held_out','Held-out test'),('stop_rule','Stop rule'),
            ('computational_class','Computational class'),('allowed_conclusion','What may be concluded'),
            ('forbidden_conclusion','What must NOT be concluded')]
    for t in tasks:
        sections.append(f'## {t["id"]} — {t["title"]}\n\n')
        sections.append(f'Priority: **{t["priority"]}**, score {t["priority_score"]}; wave {t["wave"]}.\n\n')
        sections.append('- **Dependencies:** '+(', '.join(t['dependencies']) or 'None')+'.\n')
        for key,label in labels:sections.append(f'- **{label}:** {t[key]}\n')
        sections.append('- **Repo de-dup / prior:** '+'; '.join(reference(p) for p in t['repo_dedup']['prior_work'])+'.\n')
        sections.append('- **Genuinely new atom:** '+t['repo_dedup']['new_atom']+'\n')
        sections.append('- **Artifacts:** '+', '.join('`fin_post_phy_campaign/'+p+'`' for p in t['artifacts'])+'.\n')
        sections.append('- **Activation:** '+t['activation']+'.\n')
        sections.append('- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** '+str(t['scores'])+'.\n\n')
    (ROOT/'FIN_POST_PHY_TASK_CATALOG.md').write_text(''.join(sections))
    dedup=['# FIN post-PHY de-dup and source index\n\n',
           f'Date: {DATE}. Local HEAD `{HEAD}`; source/guardrail bytes are preserved.\n\n',
           'A repo-wide **content** search was performed with rg, not only filename matching. ',
           'Search coverage and full file lists are in [search_index.json](fin_post_phy_review/search_index.json). ',
           'Counts overlap and are not counts of independently inspected proofs. Ignored/binary files were excluded; ',
           'the 565 logical text sources in the content-addressed handoff were indexed separately.\n\n',
           '| Search class | Matching files |\n|---|---:|\n']
    for k,v in scan['groups'].items():dedup.append(f'| {k} | {v["file_count"]} |\n')
    dedup.append('\n## Canonical campaign precedence\n\nC13 supersedes overlapping C12 numerical designs. HANKEL-03B is already executed; its nine-checkpoint guarantee is not an all-reset-time power theorem. The first PHY-001–031 queue is not reopened.\n\n| Alias | Master campaign |\n|---|---|\n')
    for a,c in aliases.items():dedup.append(f'| {a} | [{c}](FIN_PHY_MASTER_HANDOFF_20260923/campaigns/{c}/INDEX.md) |\n')
    dedup.append('\n## Task-by-task non-duplication\n\n| New task | Prior exact programme/file | Increment beyond prior work |\n|---|---|---|\n')
    for t in tasks:
        dedup.append(f'| {t["id"]} | '+ '; '.join(reference(p) for p in t['repo_dedup']['prior_work'])+' | '+t['repo_dedup']['new_atom']+' |\n')
    dedup.append('\n## Novelty limitation\n\nThe review checks the closest located FIN results and explicit master precedence; it is not a proof that no differently worded lemma exists anywhere. Before costly execution, inspect the cited originals and the targeted search class. If the exact target is found, mark DUPLICATE and replace it only with a genuinely stricter domain, better certificate, new observation or sharper impossibility statement. Do not count rerunning a producer as a new theorem.\n')
    (ROOT/'FIN_POST_PHY_DEDUP_MATRIX.md').write_text(''.join(dedup))
    selected={}
    for t in tasks:
        for p in t['repo_dedup']['prior_work']:
            if p.startswith('C') and ':' in p:
                key=aliases[p[:3]]+'::'+p[4:];selected[p]=logical[key]
    (HERE/'selected_sources.json').write_text(json.dumps(selected,indent=2)+'\n')
    main_report=ROOT/'FIN_POST_PHY_RESEARCH_PROGRAM.md'
    if main_report.exists():
        text=main_report.read_text()
        start=text.index('## F. Kill-test matrix')
        end=text.index('## J. If I had 100')
        (ROOT/'FIN_POST_PHY_KILL_SOURCE_GATES.md').write_text(
            '# FIN post-PHY kill, source and physical gates\n\n'
            'Generated from sections F–I of [the main programme](FIN_POST_PHY_RESEARCH_PROGRAM.md). '
            'Research design only; no new physical claim is promoted.\n\n'+text[start:end])
    print(json.dumps({'tasks':len(tasks),'acyclic':True,'first_batch':FIRST,
                      'top10':[(t['id'],t['priority_score']) for t in ranked[:10]],
                      'priorities':dict(Counter(t['priority'] for t in tasks))},indent=2))

if __name__=='__main__':main()
