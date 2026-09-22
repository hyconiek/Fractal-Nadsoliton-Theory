"""Bounded read-only evidence audit for the post-MP7 research design.

No old producer script is executed and no mathematical claim is auto-promoted.
The only output is this new campaign's evidence ledger.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter, deque
import hashlib, json, math

ROOT=Path(__file__).resolve().parents[1]
PACK=ROOT/'fin_rank7_mathphysics_next'
HERE=Path(__file__).resolve().parent
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p): return json.loads(p.read_text())

def main():
    out={'scope':'Bounded integrity, DAG, selected certificate and exact arithmetic checks; no large replay; not repository-wide acceptance',
         'agents_sha256':sha(ROOT/'AGENTS.md')}
    present=0;missing=[];changed=[]
    for line in (PACK/'MANIFEST.sha256').read_text().splitlines():
        if not line.strip(): continue
        want,rel=line.split(maxsplit=1);p=PACK/rel.lstrip('*')
        if not p.exists(): missing.append(rel)
        else:
            present+=1
            if sha(p)!=want: changed.append(rel)
    out['package_manifest']={'present':present,'missing':missing,'changed':changed}
    claims=load(PACK/'CLAIM_REGISTER.json')['claims'];bad=[];checked=0
    for c in claims:
        for rel,want in c['artifact_hashes'].items():
            checked+=1
            if not (PACK/rel).exists() or sha(PACK/rel)!=want:bad.append([c['id'],rel])
    out['claim_hashes']={'entries':checked,'mismatches':bad}
    routes={'inputs/r7n/':ROOT/'FIN_R7N_HANDOFF_20260920',
            'inputs/r7o3/':ROOT/'FIN_R7O3_TARGETP_HANDOFF_20260920',
            'inputs/r7p/':ROOT/'fin_rank7_followup','inputs/intake/':ROOT/'fin_rank7_intake_review',
            'work/':PACK}
    mappings=[]
    for line in (PACK/'replay/MP7-045_INPUT_MANIFEST.sha256').read_text().splitlines():
        if not line.strip():continue
        want,rel=line.split(maxsplit=1)
        prefix=next((k for k in routes if rel.startswith(k)),None)
        p=routes[prefix]/rel[len(prefix):] if prefix else None
        if rel.startswith('inputs/r7n/inputs/intake_review_20260919/'):
            p=ROOT/'fin_rank7_intake_review'/Path(rel).name
        mappings.append({'path':rel,'exists':p is not None and p.exists(),
                         'match':p is not None and p.exists() and sha(p)==want,
                         'expected':want,'actual':sha(p) if p is not None and p.exists() else None})
    out['clean_inputs']={'count':len(mappings),'predecessors':sum(not r['path'].startswith('work/') for r in mappings),
                         'issues':[r for r in mappings if not r['match']]}
    graph=load(PACK/'audit/MP7-043_dependency_graph.json');ids={n['id'] for n in graph['nodes']}
    indeg=Counter({k:0 for k in ids});children={k:[] for k in ids}
    for a,b in graph['edges']:
        assert a in ids and b in ids
        indeg[b]+=1;children[a].append(b)
    queue=deque(sorted(k for k in ids if not indeg[k]));order=[]
    while queue:
        a=queue.popleft();order.append(a)
        for b in children[a]:
            indeg[b]-=1
            if not indeg[b]:queue.append(b)
    out['dag']={'nodes':len(ids),'edges':len(graph['edges']),'acyclic':len(order)==len(ids),
                'semantic_correction':'Split MP7-034 local orbit-stabilizer lemma from global phase-weight theorem. The former, independent of GLOBAL_TRANSITION, supplies its multiplicity.',
                'stale_open_atoms':graph.get('open_quantitative_atoms')}
    a=load(PACK/'results/MP7-016_g37_global_exhaustion.json')
    b=load(PACK/'results/MP7-017_first_global_transition.json')
    tree=a['tree_records'];index={r['path']:r for r in tree}
    assert len(index)==len(tree)==253 and 'R' in index
    for p,r in index.items():
        if r['status']=='SPLIT':assert p+'0' in index and p+'1' in index
        else:assert p+'0' not in index and p+'1' not in index
    assert sum(r['status']=='SPLIT' for r in tree)==126
    assert a['uniform_contraction_q_upper']<1 and b['uniform_banach_q_upper']<1
    assert all(v[0]>0 for v in a['energy_nonzero_root_intervals'].values())
    assert b['saddle_energy_lower_over_event_box']>0
    assert all(r['strict_inclusion'] for r in b['stationary_branch_tubes'])
    out['global_certificates']={'g37_tree_records':len(tree),'event_records_reported':b['global_exhaustion_records'],
        'event_full_tree_serialized':'tree_records' in b,
        'scope':'Checks recorded structure/signs and inspected code; accepted prior local replay is reused, not repeated'}
    coex=load(ROOT/'fin_rank7_followup/certificates/R7P-026_equal_energy_event.json')
    checks=[]
    for mode,n in [(3,6),(4,6),(5,6),(6,12)]:
        lo,hi=map(F,coex['spectral_intervals'][str(mode)])
        dl=math.nextafter(float(lo)/n,-math.inf);dh=math.nextafter(float(hi)/n,math.inf)
        checks.append({'mode':mode,'dlo_encloses':F(dl)<=lo/n,'dhi_encloses':F(dh)>=hi/n,
                       'g37_product_lo':F(3.7*dl)<=F(37,10)*lo/n,
                       'g37_product_hi':F(3.7*dh)>=F(37,10)*hi/n})
    out['selected_exact_float_parameter_checks']=checks
    margin=F(6500000000000000000000000000,121203210760018863485060548125326941)
    threshold=F(67,250)-margin
    out['rounding_correction']={'exact_ceiling':str(threshold),
          'safe_decimal':'0.267999946371058',
          'claim_C06_decimal_0.267999946371_is_below_proved_ceiling':F('0.267999946371')<threshold}
    out['selected_input_hashes']={str(p.relative_to(ROOT)):sha(p) for p in [
        ROOT/'AGENTS.md',PACK/'HANDOFF.md',PACK/'REPORT.md',PACK/'CLAIM_REGISTER.json',
        PACK/'NONCONCLUSIONS.md',PACK/'NEXT_ATOMS.md',PACK/'REPLAY.md',PACK/'AGENTS_PROPOSED_PATCH.md',
        ROOT/'fin_rank7_intake_review/MP7_ANALYTIC_INTAKE.md',ROOT/'fin_r7o3_review/verification.json']}
    assert not changed and not bad and out['dag']['acyclic']
    (HERE/'evidence_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='selected_input_hashes'},indent=2))

if __name__=='__main__':main()
