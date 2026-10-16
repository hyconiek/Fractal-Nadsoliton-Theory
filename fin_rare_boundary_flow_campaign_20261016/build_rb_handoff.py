#!/usr/bin/env python3
"""Generate only RB campaign bookkeeping. Does not change any predecessor/AGENTS."""
import pathlib,json,hashlib,zipfile,datetime
P=pathlib.Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
write=lambda path,content:(path.parent.mkdir(parents=True,exist_ok=True),path.write_text(content,encoding='utf8'))
jwrite=lambda path,data:write(path,json.dumps(data,ensure_ascii=False,indent=2,sort_keys=True)+'\n')
inputs={
'archive AGENTS snapshot (not current repository)': '/mnt/data/files/AGENTS.md',
'FP frozen partition': '/mnt/data/fin_count_partition_feasibility_campaign_20261016/FP-001/FROZEN_PARTITION.json',
'FP exact source/domain': '/mnt/data/fin_count_partition_feasibility_campaign_20261016/FP-001/EXACT_DOMAIN_AND_SOURCE.md',
'FP primal/dual': '/mnt/data/fin_count_partition_feasibility_campaign_20261016/FP-002/CONDUCTANCE_CONSTRUCTION_OR_DUAL.md',
'FP necessary cap': '/mnt/data/fin_count_partition_feasibility_campaign_20261016/FP-002/POSITIVITY_AND_CAP_PROOF.md',
'BA parity flux': '/mnt/data/fin_reversible_bracket_averaging_campaign_20261016/BA-001/COMPONENT_LEDGER.md',
'BA flux contract': '/mnt/data/fin_reversible_bracket_averaging_campaign_20261016/BA-001/CLASS_AND_FLUX_CONTRACT.json',
'KR theorem': '/mnt/data/fin_stationary_kinetic_rigidity_campaign_20261016/KR-002/STATIONARY_CLASS_PROCESS_THEOREM_OR_OPEN_LEMMA.md',
'MT original available archive': '/mnt/data/files/FIN_MT_001_MT_005_FULL_HANDOFF_20261014.zip',
'UF original available archive': '/mnt/data/files/FIN_UF_001_UF_005_FULL_HANDOFF_20261014.zip',
'FP source archive': '/mnt/data/FIN_FP_001_FP_004_ONLY_HANDOFF_20261016.zip',
'BA source archive': '/mnt/data/FIN_BA_001_BA_004_ONLY_HANDOFF_20261016.zip',
'KR source archive': '/mnt/data/FIN_KR_001_KR_004_ONLY_HANDOFF_20261016.zip',
}
missing=[
 'fin_count_partition_feasibility_review_20261016/INTAKE_20261016.md',
 'fin_count_partition_feasibility_review_20261016/PROOF_AUDIT_AND_CAP_ENCLOSURE.md',
 'fin_count_partition_feasibility_review_20261016/BOUNDARY_LEAKAGE_AND_IDEAL_CUT_SEED.md',
 'fin_count_partition_feasibility_review_20261016/BOUNDARY_DECISION.md',
 'fin_reversible_bracket_averaging_review_20261016/COUNT_PARTITION_FEASIBILITY_SEED.md',
 'fin_reversible_bracket_averaging_review_20261016/PROOF_AUDIT_AND_SCOPE.md',
 'live checked-out AGENTS.md (historical Library snapshot available only)',
 'current checked-out project source tree and original user pasted plan bytes',
]
existing=[{'alias':k,'path':v,'sha256':sha(v),'access':'actually locally available bytes'} for k,v in inputs.items() if pathlib.Path(v).exists()]
# Read original archived MT statements, hash member BYTE contents, do not redistribute predecessor files.
arch_members={
 'MT original available archive': ['fin_critical_metastability_campaign_20261014/MT-001_COERCIVITY_CONNECTIVITY.md','fin_critical_metastability_campaign_20261014/MT-003_CONDITIONAL_GIBBS.md'],
 'UF original available archive':['fin_uniform_critical_campaign_20261014/UF-002_PARITY_COEFFICIENTS_AND_CONTROLS.md']
}
members=[]
for archive, names in arch_members.items():
 with zipfile.ZipFile(inputs[archive]) as z:
  for n in names:
   members.append({'archive':archive,'member':n,'sha256':hashlib.sha256(z.read(n)).hexdigest()})
jwrite(P/'SOURCE_PROVENANCE.json',{'available':existing,'missing_exact_files':missing,'archive_members_inspected':members,'plan_scope':'User-pasted RB and FP intake supplied in conversation; no file-byte custody assertion','research_authorization':'RB intrinsic mathematics only; no new diagnostic law or physical source approval'})
write(P/'INPUTS.sha256',''.join(f"{e['sha256']}  {e['path']}\n" for e in existing))

stage={
'RB-001':{
 'execution':'EXECUTED','science':'EXACT_FINITE_GRAPH_CONTRACT_AND_CUT_DECOMPOSITION','summary':'Bound H_N to actual MT departure E_N and derived exact native-flow defect for every lower ideal, plus automatic global-level cut positivity.','claims':['Exact source s=(N/2)W d0','H_N support and D12 symmetry','Native defect identity for any lower ideal'],'not_proved':['True all-N H_N ideal sign','Actual seed N>=10^8 directional lemma independently checked','Any diagnostic conductance assignment']},
'RB-002':{
 'execution':'EXECUTED','science':'OPEN_ACTUAL_Q12_ALL_IDEAL_SIGN','summary':'Finite max-flow/min-cut equivalence and exact 12-vertex negative control, but neither actual all-N sign theorem nor native cut witness.','claims':['Finite H downset maxflow characterization','Exact native deleted-edge defect criterion','Synthetic rational 12-vertex irregular cut -9/10'],'not_proved':['All actual q12 H_N downsets pass','One negative ideal in actual q12 H_N','Any uniform-cap generator']},
'RB-003':{
 'execution':'NOT_ACTIVATED','science':'DEPENDENT_ON_RB_002','summary':'No positive capacity-free route exists yet, so no cap or IBA theorem is claimed.','claims':['Only inherited if-feasible rare-leak implication'],'not_proved':['Any sequence-uniform C','IBA obstruction','Actual q12 diagnostic family']},
'RB-004':{
 'execution':'EXECUTED','science':'OPEN_IDEAL_OR_UNIFORM_CAP_LEMMA','summary':'Exact task status and dependency ledger preserved.','claims':['No physical status change','KR under NC/IBA remains intact'],'not_proved':['Universal IBA','Native IBA counterfamily','Activity source']}
}
for name,data in stage.items():
 d=P/name
 jwrite(d/'RESULTS.json',{'id':name,'status':data['science'],'execution':data['execution'],'model':'q12/A7 exact g_c, native MT D_N, diagnostic H_N only','inputs':{'supplied':'User-provided frozen RB plan, upstream model and exposed local files','derived':data['claims'],'new_modelling_premise':False,'fitted':False,'empirical':False},'scope':{'initial':'stationary pi_D,N for hypothetical diagnostic','clock':'native model clock; fixed sqrt(N) comparison only if a family exists','support':'H_N only; source on true D_N','rate_cap':'not established'},'claims':data['claims'],'open_or_disallowed':data['not_proved'],'next':'RB-002 actual all-ideal sign proof or witness, then uniform-cap gate'})
 write(d/'REPORT.md',f"# {name} — {data['science']}\n\n{data['summary']}\n\n"+'\n'.join('- '+x for x in data['claims'])+'\n')
 write(d/'PROOF_OR_OPEN_LEMMA.md',f"# {name} proof/open lemma\n\n"+('See `EXACT_SOURCE_AND_CUTS.md` for the full exact native cut identity.' if name=='RB-001' else 'See `ORDER_IDEAL_PROOF_OR_WITNESS.md` for the proof and precise open native lemma.' if name=='RB-002' else 'NOT_ACTIVATED; see `NOT_ACTIVATED.md`.' if name=='RB-003' else 'See root `NEXT_OPEN_ATOM.md` and `DEPENDENCY_DAG.md`; no all-N theorem or countermodel was claimed.')+'\n')
 write(d/'NONCONCLUSIONS.md', '# '+name+' nonconclusions\n\n'+'\n'.join('- NOT PROVED: '+x for x in data['not_proved'])+'\n\nNo physical model validation or primary generator replacement.\n')
 write(d/'REPLAY.md', '# '+name+' replay\n\nRun `PYTHONDONTWRITEBYTECODE=1 python3 verify_rb.py` from the root of this extracted package and check `EXACT_DIAGNOSTICS.json`. This is a deliberately finite toy/elementary identity checker, not an all-N q12 solver. Review the analytic proof directly. Verify archived inputs separately using the paths and hashes in root `SOURCE_PROVENANCE.json`; originals need not accompany this incremental handoff.\n')
 write(d/'INPUTS.sha256',''.join(f"{e['sha256']}  {e['path']}\n" for e in existing))

jwrite(P/'LEAK_GRAPH_CONTRACT.json',{
 'campaign':'RB_20261016','original_parent':'q12 A7 native LOO HB; exact g_c=12/Lambda6; N multiple of 12','domain':'actual MT connected local sublevel D_N','departure':'E_N={n in D_N: exists native legal n->n-prime with n-prime outside D_N}','fixed_partition':{'epsilon':'1/8','Z':'N^(1/4)*m','R':'2(C^2+S^2)/N','C':'sum_j n_j cos(pi*j/2)','S':'sum_j n_j sin(pi*j/2)','p':'1/2 + exp(-Z^2)/8','T':'-2 log(1-p)/(1-Lambda3/Lambda6)','A':'R<=T (ties retained)','B':'R>T'},'H_N':'legal one-copy inside D_N edges within a phase, plus cross-phase if an endpoint is in E_N','drift_source':'s_n=(N/2)*W_n*L0_N m(n)','orientation':'from lower parity to higher parity','allowed_new_rates':'nonnegative reversible count rates q=c/W, only H_N, with total rate<=CN for one sequence-uniform finite C','status':'allowed graph specified; actual rates NOT CONSTRUCTED'})
jwrite(P/'NONNEGATIVE_FLOW_STATUS.json',{'RB-002':'OPEN_ACTUAL_Q12_ALL_IDEAL_SIGN','all_N_q12_primal':False,'all_N_q12_negative_ideal':False,'toy_negative_ideal':{'graph':'synthetic 12-vertex exact rational, NOT actual q12','source_sum':'-9/10','dual':{'phi':'1_{a0,a1}','beta':'0'}},'criteria':{'capacity_free':'total s=0 and s(L)>=0 for every H-lower ideal','rate_cap':'separate CNW vertex constraints; not paid'},'FP_strict_status':'OPEN_UNIFORM_FLOW_LEMMA'})

jwrite(P/'CLAIM_REGISTER.json',{
 'RB001-A':{'status':'EXACT_FINITE_THEOREM','claim':'For every H_N-lower ideal, signed native source equals positive H-up flow + forbidden-up flow - forbidden-down flow.'},
 'RB001-B':{'status':'EXACT_FINITE_COROLLARY','claim':'Every global parity prefix has nonnegative native source on any restricted allowed support.'},
 'RB002-A':{'status':'EXACT_FINITE_THEOREM','claim':'Nonnegative capacity-free H-flow iff every H-lower ideal has nonnegative source; toy exact negative dual verifies necessity of irregular ideals.'},
 'RB002-Q12':{'status':'OPEN','claim':'No verified all-N q12/A7 ideal signs or native cut countercertificate.'},
 'RB003':{'status':'NOT_ACTIVATED','claim':'No admissible uniformly capped H_N count generator, no IBA counterexample.'},
 'PHYSICAL':{'status':'UNCHANGED_OPEN','claim':'Physical source/units/energy/locality and QM/GR/ToE unchanged.'}
})
write(P/'REPORT.md','''# FIN rare boundary flow RB-001–RB-004 — report

**Verdict: OPEN_IDEAL_OR_UNIFORM_CAP_LEMMA** (first unpaid gate: exact sign for every native H_N lower ideal).

- RB-001 derives a new exact microscopic deleted-edge native source identity, demonstrates parity-prefix sign automation, and binds the MT departure layer and same FP phase rule.
- RB-002 proves the finite max-flow/ideal characterization and builds an exact 12-vertex **synthetic** negative-control certificate where all ordinary phase-by-level cuts are positive, yet an irregular ideal has sum -9/10. There is no actual native q12/A7 all-N flow construction or negative ideal.
- RB-003 **NOT_ACTIVATED**. There is no finite uniform rate cap or temporal IBA counterfamily.
- RB-004 records the open gate without perturbing the original strict FP status, KR sufficient theorem or physical-source gates.

Finite source checks in verify_rb.py have scope confined to integer Fourier/parity conventions and the synthetic rational cut certificate; none proves the actual all-N positivity statement. Independently audit the exact mathematical proofs and unavailable newer review files before admission.
''')
write(P/'HANDOFF.md','''# Incremental handoff — FIN RB-001–RB-004

Date: 2026-10-16. Destination: `fin_rare_boundary_flow_campaign_20261016/`.

**Terminal architectural disposition: `OPEN_IDEAL_OR_UNIFORM_CAP_LEMMA`.** RB-001/RB-002/RB-004 executed; RB-003 NOT_ACTIVATED. No changes to primary FIN, historical review results, AGENTS or physical status.

## New results and exact first blocker

1. Exact native deleted-edge identity and automatic global-parity-prefix sign theorem (`RB-001/EXACT_SOURCE_AND_CUTS.md`).
2. Native ideal-sign feasibility gate plus exact rational twelve-vertex failure of coarse cuts (`RB-002/ORDER_IDEAL_PROOF_OR_WITNESS.md`). The synthetic failure is **not** actual FIN.
3. The actual all-N q12 `sum_{L} s>=0` question remains open. No uniform cap attempt is licensed until capacity-free feasibility is paid. Thus no leaky IBA counterfamily is submitted.

The inherited rare-leak static/temporal implication is only conditional on a still-missing positive-rate construction. Rejected/failing examples and missing source originals are documented explicitly.

## Dependencies and provenance

`SOURCE_PROVENANCE.json` and `INPUTS.sha256` bind actually accessible original source files, archived MT/UF members and the previous incremental packages. New FP-review `BOUNDARY_LEAKAGE_AND_IDEAL_CUT_SEED.md`, new BA/FP review completions and current checked-out `AGENTS.md` were **not** accessible as verifiable bytes; the Library AGENTS is a historical snapshot. The current user's pasted plan was read as the authorization, not assigned an invented file identity. Scientific claims in the actual unavailable review seed were reused only with explicit inherited status.

See `REPLAY.md`, `CLAIM_REGISTER.json`, `NEXT_OPEN_ATOM.md`, `NONCONCLUSIONS.md`, `PROPOSED_AGENTS_PATCH.md`. This ZIP contains **only new RB** products; it does not bundle old campaigns or predecessor research.
''')
write(P/'REPLAY.md','''# Independent replay and integrity

Run from extracted directory (no network necessary):

```bash
PYTHONDONTWRITEBYTECODE=1 python3 verify_rb.py
PYTHONDONTWRITEBYTECODE=1 python3 verify_package.py
```

`verify_rb.py` checks integer Fourier/parity conventions for the q12 label map and an exact synthetic 12-vertex rational native-flow/dual fixture. Its result is bounded finite diagnostics and does not establish any true q12 large-N feasibility theorem.

`verify_package.py` verifies manifest SHA-256 coverage of all RB artifacts other than self-referential metadata. Source scientific documents are *not* silently included; to audit them, retrieve the exact upstream bytes and compare with the digests in `SOURCE_PROVENANCE.json`/`INPUTS.sha256`. Missing newer review sources cannot be proved identical by this package.

The only accepted model-specific new formula is an exact finite identity conditional on the already specified native generator and support. The general directed graph maxflow theorem is proved in RB-002. The unresolved all-N native ideal-sign property remains `OPEN`.
''')
# the rebuilt verification is deterministic, no dependency on source archives
r=__import__('subprocess').run(['python3',str(P/'verify_rb.py')],check=True,capture_output=True,text=True)
write(P/'VERIFY_LOG.txt',r.stdout)

# Actual manifest includes all deliverables, but excludes itself and machine-generated package-check metadata
excluded={'MANIFEST.sha256','PACKAGE_VERIFICATION.json'}
paths=sorted([x for x in P.rglob('*') if x.is_file() and x.name not in excluded and '__pycache__' not in x.parts])
write(P/'MANIFEST.sha256',''.join(f'{sha(x)}  {x.relative_to(P)}\n' for x in paths))
print('WROTE',len(paths),'manifest entries')
