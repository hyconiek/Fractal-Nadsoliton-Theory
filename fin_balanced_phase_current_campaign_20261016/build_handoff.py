from pathlib import Path
import json, hashlib, zipfile, textwrap, subprocess, sys, os
D=Path(__file__).parent

def put(name,text):
 p=D/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(textwrap.dedent(text).strip()+'\n',encoding='utf-8')

def putj(name,obj):
 p=D/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(obj,indent=2,sort_keys=True,ensure_ascii=False)+'\n')

# Record exact available source hashes (no retroactive claims for unavailable predecessors).
local={
 'ARCHIVED_AGENTS_NOT_CURRENT':'/mnt/data/files/AGENTS.md',
 'SC001_LOCAL_THEOREM':'/mnt/data/fin_signed_interface_orbit_cut_campaign_20261016/SC-001/INCOMING_COMPENSATION_PROOF_OR_OPEN_LEMMA.md',
 'SC_HANDOFF':'/mnt/data/fin_signed_interface_orbit_cut_campaign_20261016/HANDOFF.md',
 'RB_WEIGHT_RATIO':'/mnt/data/fin_rb_orbit_weighted_degree_campaign_20261016/UNIFORM_NATIVE_EDGE_ELLIPTICITY.md',
 'RB_ORBIT_LIFT':'/mnt/data/fin_rb_orbit_flow_reduction_campaign_20261016/ORBIT_QUOTIENT_AND_LIFT_THEOREM.md',
 'FP_CONDUCTANCE':'/mnt/data/fin_count_partition_feasibility_campaign_20261016/FP-002/CONDUCTANCE_CONSTRUCTION_OR_DUAL.md',
}
inputs={}
for k,v in local.items():
 p=Path(v)
 inputs[k]={'path':v,'available':p.is_file(),'sha256':hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None}
mtz=Path('/mnt/data/files/FIN_MT_001_MT_005_FULL_HANDOFF_20261014.zip')
if mtz.exists():
 with zipfile.ZipFile(mtz) as z:
  for n in ['fin_critical_metastability_campaign_20261014/MT-001_COERCIVITY_CONNECTIVITY.md','fin_critical_metastability_campaign_20261014/MT-003_CONDITIONAL_GIBBS.md']:
   raw=z.read(n);inputs[n]={'archive':str(mtz),'archive_member':n,'sha256':hashlib.sha256(raw).hexdigest(),'available':True}
for code,name in [
 ('PC01','fin_signed_interface_orbit_cut_review_20261016/INTAKE_20261016.md and PROOF_AUDIT.md'),
 ('PC02','fin_balanced_phase_cut_research_20261016/ACTUAL_COUNTERCUT_TO_STRONG_SHIELDING.md'),
 ('PC03','fin_balanced_phase_cut_research_20261016/BALANCED_PHASE_DEFECT_THEOREM.md and COMPACT_PHASE_PREFIX_LOCAL_LIMIT.md'),
 ('PC04','fin_balanced_phase_cut_research_20261016/EXACT_COMPACT_CURRENT_MARGINALS_AND_GLUING.md and CURRENT_CARRY_CONTRACT.json'),
 ('PC05','fin_balanced_phase_cut_research_20261016/RESULTS.json and MANIFEST.sha256'),
 ('PC06','fin_rbcut_weighted_degree_review_20261016/PROOF_AUDIT_AND_CORRECTIONS.md'),
 ('PC07_review','fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md and fin_stationary_kinetic_rigidity_review_20261016/PROOF_COMPLETIONS_AND_SCOPE.md'),
 ('PC08','fin_count_partition_feasibility_review_20261016/BOUNDARY_LEAKAGE_AND_IDEAL_CUT_SEED.md'),
]:
 inputs[code]={'reference':name,'available':False,'sha256':None,'note':'User supplied plan and summary but not matching source bytes; not retrospectively inspected.'}
putj('SOURCE_PROVENANCE.json',{'actual_hashed_inputs':inputs,'current_AGENTS_available':False,'PC03_PC04_source_bytes_available':False,'model_data_from':'Current user message plus accessible predecessor snapshots. Source-dependent statements explicitly conditional.'})
with (D/'INPUTS.sha256').open('w') as f:
 for key,row in sorted(inputs.items()):
  if row.get('sha256'): f.write(f"{row['sha256']}  {key}\n")
putj('MARGINAL_AND_SUPPORT_CONTRACT.json',{
 'model':'q12 A7 native LOO reflected MT; diagnostic H; not primary generator',
 'lambda':0,'g':'12/Lambda6 exact','N':'multiples of 12','epsilon':'1/8','z_band':'adjacent z_l,z_l+1 in [-1,1]',
 'phase_A':'R=2(C^2+S^2)/N <= T(z)=-(2/gamma3)*log(1-[1/2+exp(-z^2)/8])',
 'phase_B':'strict complement inside actual MT D_N',
 'source':'s_x=(N/2)W_x L0^D m(x)',
 'fixed_incoming':'I_P,l(x)=J_P,l-1 * W_x/W_P,l',
 'fixed_correction':'a_x=[-s_x-I_P,l(x)]_+ on E; I_E\'=I+a; I_bulk\'=I-(sum_E a)*W_x/W_bulk; O\'=I\'+s',
 'edge_support':'native odd->even single-copy edges, within phase; phase crossings only incident to actual E_N',
 'target':'O_l\'(T) <= I_l+1\'(N_H(T)) for all T subset of entire X_l',
 'result':'only singleton and bounded-overlap nondeparture subclasses proved conditionally',
 'dependency':'uniform compact phase-prefix current ratio inherited from unavailable PC03/04 originals; uniform E_N fibre suppression inherited MT/PC03'
})
put('ELEMENTARY_TRANSPORT_TESTS.md', '''
# PC-001 — elementary actual-fibre checks

The frozen row/column identities are kept **exact**: every phase P at level l has incoming total `J_{P,l-1}`, outgoing total `J_{P,l}`, and their difference is the exact sum of microscopic sources `s_x`. Corrections `a_x` merely redistribute incoming marginals within an individual phase fibre, preserving their sum. Edge support is full H, including E_N-incident phase crossings. Zeros and `R=T(z)` phase-A ties are retained. This accounting is algebraic and does not by itself supply c_l.

**Actual geometric check.** Every upward odd→even count replacement changes the integer Fourier pair by exactly one of `(±1,±1)`; 9 distinct odd/even label pairs realize each. At x∉E_N, the SC shell argument gives at least 9 retained same-phase H successors. Exact Gibbs ratios on these edges exceed .9 for all sufficiently large N. See `PARTIAL_HALL_PROOF.md`.

**Marginal window check (dependent on PC03/PC04).** The supplied compact phase-prefix asymptotic and an analytically bounded quartic integral imply the incoming bulk interval [.18,.35]×NW and all-column lower .18×NW. E corrections cannot spoil the column lower, and their compact-fibre subtraction is negligible by the actual MT departure exponential estimate. No phase is discarded for small mass.

**Negative check.** Local degree, singleton Hall or the positivity of the marginals does not establish full weighted Hall. `verify_pc.py` supplies a small exact abstract example with individually feasible singletons and a failing doubleton; it is *not* a FIN generator.

The user-requested original PC01–PC05 file bytes were not mounted in this execution, so their hashes, original numerical margins and independently checked PC04 source proof are not asserted. The new result is conditional on those precise inherited local-limit inputs.
''')
put('HALL_PROOF_OR_ACTUAL_COUNTERFAMILY.md', '''
# PC-002 — partial actual-H Hall proof, full question OPEN

**Verdict:** `OPEN_ADJACENT_FIBRE_COUPLING` for the required quantifier **every** T⊆X_l. No genuine q12 counterfamily has been constructed.

**New nontrivial positive theorem in the actual model, conditional on the inherited PC03/PC04 compact fibre local limit:** Weighted Hall holds for each singleton T={x} with x outside E_N on the frozen z band. It also holds for any set T⊆X_l\E_N for which one can select nine phase-safe real H-up successors per x with **at most four selected incoming incidences at each right vertex**. This uses the true full-H neighbor set, exact source/corrections and MT weights, not a proxy cylinder or within-phase-only replacement.

**Proof and strict rational margins:** `PARTIAL_HALL_PROOF.md`, Sections 1–4; `verification.json`. For these T, LHS≤357/1000 NW(T) and RHS≥729/2000 NW(T), leaving exact positive margin 3/400 NW(T). The ratio 0.9 is an intentionally weak lower bound for the actual Gibbs ratio on every edge.

**Remaining proof gate:** Show this overlap condition (or a different signed source-sensitive transport inequality) sufficiently controls *every* dense subset T with its actual prescribed marginals, or produce an explicit true q12/A7 sequence of violating adjacent-fibre sets on an unbounded N subsequence. Nothing from the limited proof covers arbitrary dense T.

**Failure controls:** Singleton degree alone does not prove Hall; E_N source corrections and outward domain rejection cannot be deleted; a failure for a stricter within-phase-only graph would not refute full H. No failure of the *frozen ansatz* was found, much less an obstruction to arbitrary H-supported rates.
''')
put('COUNTERTESTS.md', '''
# Negative controls and tested hypotheses

1. **Overlapping successors:** The exact abstract 2-left/2-right fixture has left masses 1,1, right masses 3/2,1/2 and both left vertices connected only to the first right vertex. Each singleton passes weighted Hall, but the doubleton fails: 2>3/2. This is a logic check, **not a frozen-q12 counterfamily**.
2. **No source relabelling:** Using native full-fibre upward conductance row sums instead of `I',O'` would solve a different problem. Not done.
3. **No omission of E:** For an E vertex, correction `a=[-s−I]_+` may be strictly positive and its upward degree may be much smaller than 9 because reflection rejects exits. Our nine-edge theorem expressly excludes E tails. E heads are retained with incoming correction.
4. **No arbitrary normalization:** The p_A(z)=1/2+exp(−z²)/8, p_B=1−p_A, J_P source and exact finite label clock remain frozen.
5. **No synthetic law:** The g=0 untruncated parity-stratum multinomial coupling is an exact baseline for the unrestricted graph. It is *not* a coupling for the g_c/MT phase-constrained H graph.
6. **No broad finite-N claim:** N0 in the compact-band local limit is unknown, no practical N or probability rate is asserted.
''')
put('PROOF_OR_OPEN_LEMMA.md', '''
# Exact unpaid PC-002 lemma

On the actual q12/A7 MT graph and the fixed alpha=1/8 phase partition, for every sufficiently large N multiple of 12 and every adjacent z-level pair in [-1,1], show

`Σ_{x∈T} O'_l(x) ≤ Σ_{y∈N_H(T)} I'_{l+1}(y)` for all T⊆X_l,

with actual `W`, `s`, `E`, full-H support (including E-crossings), and the *unmodified* PC04 corrected carries. Or exhibit an unbounded admissible N subsequence and actual subsets T_N with the opposite strict inequality. The proved subclass has a nine-safe-successor selection with at most 4-fold overlap and no left E vertices. The remaining sets can be dense enough for many chosen safe paths to converge at the same right states. No valid universal congestion estimate follows from nine-way local availability or a mean-field limit alone.

If full Hall is paid, a separate compact c_l construction and still separate full-D_N gluing are needed. No finite-N output flow is claimed by the present partial result.
''')
put('COMPACT_BAND_FLOW.md', '''
# PC-003 — NOT_ACTIVATED

The required full weighted-Hall inequalities on all subsets of every adjacent fibre were not proved. Therefore no nonnegative edge coupling c_l with exact prescribed `I',O'` has been constructed. No 5NW vertex cap, exact full generator, far-soft-level gluing, or IBA process obstruction is claimed. The positive singleton and sparse-subset inequalities are insufficient to activate this task.
''')
put('FULL_DOMAIN_EXTENSION_OR_OPEN_ATOM.md', '''
# Unpaid full-domain work after a hypothetical compact-Hall proof

1. Prove Hall for arbitrary dense T in each adjacent compact-band pair using the *original* PC04 corrected carries (or refute the frozen carry ansatz with actual N-indexed sets).
2. Only on a proved compact coupling, establish a consistent gluing to far soft levels and all actual E_N departure vertices. Preserve exact statewise drift and one N-independent incident rate cap in D_N.
3. Only after a genuine full-D_N generator family exists, invoke the inherited boundary-leak temporal IBA argument and pay its preparation/jump-bracket conditions. There is no current authorization to promote any of these implications.
''')
putj('SOURCE_AND_CAP_PAYMENT.json',{'pc003_status':'NOT_ACTIVATED','compact_full_hall':False,'actual_c_l_constructed':False,'full_domain_gluing':False,'full_domain_uniform_cap':False,'IBA_counterfamily':False,'conditional_subclass_only':True})
put('IBA_NOT_ACTIVATED.md', '''
The fixed carry has no proven full adjacent-fibre coupling, much less a full MT-domain generator; the rare-leak stationary Jensen argument is therefore **NOT_ACTIVATED**. Stationary phase predictions or limits alone cannot replace exact statewise rates. No iid preparation, natural activity, physical clock or detector evidence was derived.
''')
putj('RESULTS.json',{'verdict':'OPEN_ADJACENT_FIBRE_COUPLING','tasks':{'PC-001':'BOUND_WITH_UNAVAILABLE_SOURCE_PROVENANCE','PC-002':'PROVED_CONDITIONAL_SINGLETON_AND_OVERLAP_LE4; FULL_HALL_OPEN','PC-003':'NOT_ACTIVATED','PC-004':'SCOPED_HANDOFF'},'new_quantitative_claims':{'q12_native_weight_ratio_lower':'9/10 for retained edges; N>=10000 for ratio bound','phase_incoming_lower':'9/50 NW conditional on current local limit for sufficiently large N','phase_bulk_incoming_upper':'7/20 NW conditional','bulk_outgoing_upper':'357/1000 NW for sufficiently large N','restricted_Hall_right_lower':'729/2000 NW(T) for overlap<=4','positive_hall_margin':'3/400 NW(T)'},'full_hall':None,'actual_N_counterfamily':None,'full_generator':False,'original_PC03_PC04_available':False})
putj('CLAIM_REGISTER.json',[
 {'id':'PC-C1','claim':'Exact single-edge Gibbs weight ratio >9/10 on genuine MT edges for N>=10000','status':'DERIVED_FINITE_ALGEBRA','type':'M'},
 {'id':'PC-C2','claim':'Quartic phase current ratio bounds A:[1/5,5/16], B:[3/16,1/3]','status':'PROVED_ANALYTIC_FROM_FROZEN_PHASE_FORM','type':'M'},
 {'id':'PC-C3','claim':'Uniform compact-band incoming current window [.18,.35] NW','status':'CONDITIONAL_ON_UNAVAILABLE_PC03_PC04_LOCAL_LIMIT_AND_MT_FIBRE_TAIL','type':'M'},
 {'id':'PC-C4','claim':'Hall for singletons outside E and overlap<=4 selected-successor sets','status':'PROVED_CONDITIONAL_ON_PC-C3','type':'M'},
 {'id':'PC-O1','claim':'Full Hall for every T in whole fibre','status':'OPEN','type':'M'},
 {'id':'PC-O2','claim':'Actual compact coupling, full domain generator and IBA counterfamily','status':'NOT_ACTIVATED','type':'M'}
])
put('NONCONCLUSIONS.md', '''
# Nonconclusions

No weighted-Hall theorem for arbitrary subsets T of any whole actual q12 phase fibre; no growing-N actual violating T; no refutation of the frozen carry ansatz; no compact exact coupling; no full-domain flows, 5NW cap or kinetic counterexample to universal IBA. No late iid theorem, new FIN kinetic law, empirical records, physical source of A7/activity/clock/energy, spatial incidence, conserved physical energy, apparatus/calibration, QM/GR/SM/L_total/ToE, or ontology/selector/legacy-role changes. The source current-prefix input PC03/PC04 was not byte-verifiable in this environment; the partial theorem is conditional on the stated content and should be independently audited in the actual repo.
''')
put('NEXT_OPEN_ATOM.md', '''
# Next strictly frozen open atom

**Actual dense-subset weighted Hall for PC04 marginals.** After paying source custody, obtain either (i) a *uniform actual-q12* transport/capacity argument for `O'_l(T)≤I'_{l+1}(N_H(T))` for all whole-fibre sets T in the fixed adjacent critical band, including arbitrary overlap and E-correction, or (ii) an exact sequence of actual q12 count-level subsets violating that inequality for unbounded N. The singleton and at-most-four-overlap subclasses are closed conditionally. Do not substitute a further degree-only test or restart threshold tuning. If the carry is refuted, do not infer infeasibility of general H-supported flows. Full-D_N gluing and the IBA consequence are separate gates.
''')
put('PROPOSED_AGENTS_PATCH.md', '''
# Proposed guardrail — not applied

PC-001/002 research has a conditional actual-H partial Hall lemma: on each fixed compact critical band and subject to original PC03/PC04 uniform phase-current fibre asymptotics, singleton sets outside MT departure E and nondeparture sets admitting nine safe successors with at most fourfold right overlap satisfy weighted Hall. This does not pay full PC-002, PC-003 or IBA; no full-domain flow was constructed. Preserve the frozen PC04 marginals and all E_N crossings. Source PC03/PC04 bytes were unavailable to the executor, so independent audit must check the inherited limit before accepting even this conditional partial lemma.
''')
put('DEPENDENCY_DAG.md', '''
# Dependency graph

`frozen q12/A7, MT cylinder, exact Gibbs W` -> new edge ratio >.9 (direct algebra)
`frozen q12 Fourier R/T and SC nine-way shell geometry` -> 9 safe H successors outside E
`PC03 compact phase-prefix current limit [source not byte-accessible]` + `quartic phase bounds` -> incoming current window
`MT exit-layer mass [archived]` + `PC04 corrected I,O [user-supplied]` -> all-column lower and nondeparture row upper
`weight ratio + 9 edges + current window + <=4 overlap` -> scoped weighted-Hall lemma
`arbitrary dense overlap / E departure left tails` -> **OPEN PC-002** -> no PC-003 -> no IBA
''')
put('HANDOFF.md', '''
# FIN PC-001–PC-004 — handoff

**Date:** 2026-10-17. **Scientific disposition:** `OPEN_ADJACENT_FIBRE_COUPLING` with a **new quantitative conditional partial-Hall theorem**. Strictly no full compact transport or native-IBA counterfamily.

New exact results: native MT opposite-parity Gibbs edge-weight ratio >.9 (N≥10000); the fixed phase functions yield analytic compact current coefficient enclosures A:[1/5,5/16] and B:[3/16,1/3] once the inherited phase-current local-limit normalization is paid. On the exact frozen PC04 marginals and actual full-H support, these yield strict Hall for every nondeparture singleton and every nondeparture set whose nine selected legal upward successors overlap at most 4 times, with margin ≥3/400 N W(T). The proof keeps the original E_N corrections and phase ties.

Remaining: **full weighted Hall for *every* actual adjacent-fibre subset** or an actual growing-N adverse family. The overlap-four property is not known for dense subsets. PC-003 NOT_ACTIVATED, no full-D_N gluing and no IBA claim.

Source note: archived MT/SC/RBCUT files were hashed, but the actual new PC01–PC05, several current reviews and current AGENTS were not available; user supplied their formulas and summary in the current message. `PARTIAL_HALL_PROOF.md` marks the corresponding conditional dependence. Exact arithmetic checker is supplementary only.

Primary FIN, older sources, reviews, and AGENTS unchanged. This package contains only new PC results and a proposed patch.
''')
put('REPORT.md', '''
# PC campaign report

Tasks PC-001/002 executed analytically within source constraints; PC-002 closed **only** for restricted actual-q12 subset classes. The full target remains open. PC-003 NOT_ACTIVATED; PC-004 completed via this handoff and an explicit next atom.

**Positive result:** On compact consecutive soft levels, exact local Gibbs ratios and the inherited phasewise current window give genuine Hall slack for each nondeparture singleton and for any nondeparture T with an at-most-four incidence selection of nine legal H successors per x. No marginals, phases or domain were changed.

**Negative control:** Individual Hall tests can pass while a union fails; an abstract 2-left weighted graph certifies this logic and is not presented as FIN evidence.

**Provenance constraint:** The user's supplied plan/summary gives the PC04 formulas, but the actual review and new research input files are absent in the mounted environment; no fictitious hashes or source replay. Refer to SOURCE_PROVENANCE.json.
''')
put('RESEARCH_SUMMARY_PL.md', '''
# FIN — PC-001–PC-004, podsumowanie po polsku

Uzyskano nowe, częściowe twierdzenie dotyczące rzeczywistego grafu q12/A7 i zamrożonych marginesów prądu z instrukcji. Na wszystkich dozwolonych krawędziach lokalnej domeny MT iloraz wag stanów sąsiednich jest większy od 0,9 dla N≥10000. Cztery kombinacje przyrostów (C,S) każda posiadają dziewięć jednokopiowych realizacji. Po wykorzystaniu przyjętej, ale w tym środowisku niezweryfikowanej z oryginalnych plików granicy prądów fazowych można udowodnić ważoną nierówność Halla dla pojedynczych stanów spoza warstwy wyjścia E_N oraz dla zbiorów, w których wybrane dziewięć następców nakłada się nie więcej niż czterokrotnie. Zapewnia to dodatni margines co najmniej 3/400 razy N W(T).

WYNIK NIE WYSTARCZA dla dowolnego zbioru całego włókna, gdzie nakładanie może być znacznie większe. Nie ma pełnego sprzężenia prądów, globalnego generatora ani kontrprzykładu dla IBA. Status PC-002: OPEN_ADJACENT_FIBRE_COUPLING, PC-003: NOT_ACTIVATED. Brak oryginalnych aktualnych plików PC03/PC04 i aktualnego AGENTS do audytu, dlatego wynik zależny od ich granicy lokalnej jest warunkowy. Nie zmieniono pierwotnego FIN ani wcześniejszych pakietów.
''')
put('REPLAY.md', '''
# Replay

Run `python verify_pc.py` (Python standard library, Fraction exact arithmetic); it writes `verification.json`. The tests check rational weight/current/Hall margins, 36 Fourier directions, a strict negative Hall-control, and one exact iid q12 multinomial balance fixture. They are not a proof for all N nor a replay of unavailable source PC03/PC04.

Then use `python check_package.py <zip_path>` after package creation to validate hashes and ZIP CRC. `MANIFEST.sha256` covers emitted new campaign outputs, excluding itself and the archival ZIP wrapper, without altering original source files.
''')
put('SCOPED_VERDICT.json', '''{"decision":"OPEN_ADJACENT_FIBRE_COUPLING","partial_hall":true,"partial_hall_scope":"nondeparture singleton and nine-safe-edge subsets with overlap<=4, conditional on inherited PC03/PC04 compact current limit","full_weighted_hall":false,"PC003":"NOT_ACTIVATED","source_custody":"archived MT/SC sources checked; newer PC03/PC04 originals absent"}''')
# rerun tests
subprocess.run([sys.executable,str(D/'verify_pc.py')],check=True,stdout=subprocess.DEVNULL)
# final manifest - build script excluded to avoid self-modifying content concerns? Include it.
files=sorted(p for p in D.rglob('*') if p.is_file() and p.name not in {'MANIFEST.sha256','check_package.py'} and not p.name.endswith('.zip'))
manifest=''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(D))+'\n' for p in files)
(D/'MANIFEST.sha256').write_text(manifest)
put('check_package.py', '''
#!/usr/bin/env python3
import hashlib, pathlib, sys, zipfile
zpath=pathlib.Path(sys.argv[1]); assert zpath.is_file()
with zipfile.ZipFile(zpath) as z:
 assert z.testzip() is None
 names=set(z.namelist()); root='fin_balanced_phase_current_campaign_20261016/'
 lines=z.read(root+'MANIFEST.sha256').decode().splitlines()
 for row in lines:
  digest,name=row.split('  ',1)
  assert root+name in names, name
  assert hashlib.sha256(z.read(root+name)).hexdigest()==digest,name
 assert len(lines)>=20
 print('ZIP PASS:',len(lines),'scientific output manifest hashes; entries',len(names))
''')
# add checker separately (excluded from original manifest, hashed separately in PACKAGE_VERIFICATION.json)
putj('PACKAGE_VERIFICATION.json',{'manifest_entries':len(files),'manifest_scope':'All outputs except package checker, PACKAGE_VERIFICATION, MANIFEST; checker separately bound below','check_package_sha256':hashlib.sha256((D/'check_package.py').read_bytes()).hexdigest(),'verify_pc_run':'PASS'})
zipout=Path('/mnt/data/FIN_PC_001_PC_004_ONLY_HANDOFF_20261017.zip')
with zipfile.ZipFile(zipout,'w',zipfile.ZIP_DEFLATED,compresslevel=9) as z:
 for p in sorted(D.rglob('*')):
  if p.is_file(): z.write(p, arcname='fin_balanced_phase_current_campaign_20261016/'+str(p.relative_to(D)))
subprocess.run([sys.executable,str(D/'check_package.py'),str(zipout)],check=True)
print('zip',zipout,zipout.stat().st_size,'total files',sum(p.is_file() for p in D.rglob('*')))
