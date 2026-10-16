#!/usr/bin/env python3
"""Generate HM-only handoff from already-authored HM master proof and finite checks.
Does not modify prior evidence. Uses only files present in this directory / mounted inputs.
"""
from pathlib import Path
import json, hashlib, zipfile, os
root=Path(__file__).resolve().parent
root.mkdir(exist_ok=True)
source_paths={
 'HM user-supplied plan':Path('/mnt/data/Wklejony Markdown(20261015-233307).md'),
 'available archived AGENTS snapshot (not current commit)':Path('/mnt/data/files/AGENTS.md'),
 'MT archive':Path('/mnt/data/files/FIN_MT_001_MT_005_FULL_HANDOFF_20261014.zip'),
 'UF archive':Path('/mnt/data/files/FIN_UF_001_UF_005_FULL_HANDOFF_20261014.zip'),
 'OBS archive':Path('/mnt/data/files/FIN_OBS_C01_C05_FULL_HANDOFF_20261015.zip'),
 'NT package':Path('/mnt/data/FIN_NT_001_NT_006_ONLY_HANDOFF_20261015.zip'),
}
def sha(b):return hashlib.sha256(b).hexdigest()
provenance=[]
for desc,p in source_paths.items():
    provenance.append({'description':desc,'path':str(p),'available':p.exists(),'sha256':sha(p.read_bytes()) if p.exists() else None})
mt=source_paths['MT archive']
with zipfile.ZipFile(mt) as z:
    member='fin_critical_metastability_campaign_20261014/MT-003_CONDITIONAL_GIBBS.md'
    provenance.append({'description':'MT-003 local Gibbs source member actually inspected','archive':str(mt),'member':member,'sha256_member_bytes':sha(z.read(member))})
with zipfile.ZipFile(source_paths['NT package']) as z:
    mem='fin_native_transfer_campaign_20261015/verify_nt.py'
    provenance.append({'description':'NT model builder reused for finite checks','archive':str(source_paths['NT package']),'member':mem,'sha256_member_bytes':sha(z.read(mem))})
missing=['actual repository HEAD 3be481a7 / latest AGENTS.md','fin_native_transfer_review_20261016/INTAKE_20261016.md','fin_native_transfer_review_20261016/MEDIATOR_BRIDGE_SEED.md','fin_rank7_mathphysics_next/proofs/MP7-031_finite_N_variational_model.md','fin_rank7_mathphysics_next/proofs/MP7-032_gaussian_auxiliary_representation.md','fin_kinetic_invariants_review_20261014/INTAKE_20261014.md','fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md']
(root/'INPUT_PROVENANCE.json').write_text(json.dumps({'available':provenance,'not_verified_directly':missing,'claim':'No reconstructed historical checkout or source-review access. MT-003 source is read from the preserved archive, independent of current review. All new HM calculations specify premises.'},ensure_ascii=False,indent=2)+'\n')

stages={
'HM-001':dict(name='Exact cavity update',status='EXACT_FINITE_N_IDENTITY',brief='Completing the square gives the exact native LOO softmax after eliminating a freshly sampled cavity field. θ|n and ν_r are different; no (N−1) denominator or allowed-label renormalization.',files={
'JOINT_AND_CAVITY_CONTRACT.json':{'law':'Pi(x,θ)∝exp[−N|θ|²/(2g)+θ·X^Tn]','posterior':'N((g/N)X^Tn,(g/N)I_7)','cavity':'ν_r(θ)∝exp[−N|θ|²/(2g)+θ·X^Tr] Σ_j exp(X_j·θ)','primary':'p(j|r)=softmax_j[(g/N)Ar]','rule':'update-only θ; original denominator N; reject outgoing reflections'},
'EXACT_UPDATE_BRIDGE.md':'# HM-001 exact update bridge\n\nThe normalizing Gaussian integral for candidate j yields `exp[g/(2N)||X^T(r+e_j)||²]`. The j-dependent part is `(g/N)(Ar)_j+(g/(2N))A_jj`; constant diagonal cancels. Consequently an ephemeral ν_r draw followed by categorical conditional output gives exactly the original LOO refresh law, including null refreshes. Any attempt to integrate using ρ_n instead changes the kernel. See HM_THEOREM_PROOFS.md §HM-001 for the full derivation.\n',
'SMALL_FINITE_CHECKS.json':{'source':'verification.json','fixture':'20 fixed 11-slot cavity counts, N=12, A7 floating reconstruction','statement':'Integrated completed-square categorical probabilities equal LOO softmax to binary64 tolerance, NOT Gaussian Monte Carlo'},
}),
'HM-002':dict(name='Persistent feedback and memory',status='EXACT_NONMARKOV_PROJECTION_NO_GO',brief='Stationary joint Pi is reversible for the persistent field sampler but P exp(tL) P fails the count semigroup identity due to a nonzero second-derivative defect.',files={
'DIAGNOSTIC_MODEL_CONTRACT.json':{'ID':'AUX-HB-OU-kappa','primary_FIN_adopted':False,'kappa':[1,4],'label_attempt_rate':1,'label_proposal':'softmax(Xθ)','field_SDE':'dθ=κ((g/N)X^T n−θ)dt+sqrt(2κg/N)dB','stationary':'Pi(dn,dθ)=pi(dn)rho_n(dθ)'},
'REVERSIBILITY_AND_DOMAINS.md':'# HM-002 reversibility and domains\n\nAt fixed θ the labelled categorical generator is reversible under the product conditional Π(x|θ). At fixed labels the OU part is reversible for ρ_n. Their sum is reversible for Π; reflecting count jumps by rejection preserves Π_D. For count parity m, Lm=φ(θ)−m is bounded and L²m has at most linear θ-growth, so Gaussian L² domains justify K′(0), K′′(0) on this test.\n',
'COMPRESSION_RIGIDITY_OR_COUNTEREXAMPLE.md':'# Strict stationary count compression obstruction\n\n`||QLPm||²=E_pi Var_{rho_n}[φ(θ)]>0`; θ has Gaussian full support, and on the parity ray φ=tanh(bθ6), so φ is not constant. Since L is reversible, `<m,(K′′(0)−K′(0)^2)m>=||QLPm||²>0`, where K(t)=P e^{tL}P. Therefore K(t) is not a semigroup. It cannot equal the native count semigroup. This is a statement about the specified persistent sampler only; no universal mediator no-go follows.\n',
'MEMORY_SOURCE_WITNESS.json':{'quantity':'<m,(K_doubleprime(0)-K_prime(0)^2)m>','exact_expression':'E_pi Var_{rho_n}[Σ_j σ_j softmax_j(Xθ)]','strictly_positive':True,'scope':'global full Pi; local Pi_D when a positive-mass interior count type is included'},
}),
'HM-003':dict(name='Collective generator coefficients',status='EXACT_CHARACTERISTICS_FORMAL_CRITICAL_COEFFICIENTS',brief='18-dimensional Jacobian has one soft mode; projection gives rκ=κ/(1+κ), drift rκ(aZ−Z³/3), quadratic variation 2rκ. At this stage this is a coefficient target, not a path theorem.',files={
'JOINT_CHARACTERISTICS.md':'# HM-003 joint characteristic identities\n\nWith q=softmax(Xθ), `Lp_j=q_j−p_j`, `Lθ=κ(gX^Tp−θ)`. Jump bracket: `d<p_j,p_k>/dt=N^−1[diag(p)+diag(q)−pq^T−qp^T]_{jk}`; θ bracket `2κg I_7/N`; mixed count/θ brackets zero. For m=σ·p and h=√(Λ6/12)θ6, drift: `Lm=φ−m`, `Lh=κ[(g/gc)m−h]`. Brackets: `(2/N)(1−mφ)`, `2κ(g/gc)/N`, zero mixed. See HM_THEOREM_PROOFS.md.\n',
'SOFT_FAST_CHART.json':{'J':'[[-1,1],[kappa,-kappa]]','right_soft':[1,1],'left_soft':'[kappa/(kappa+1),1/(kappa+1)]','other_eigenvalue':'-(1+kappa)','full_tangent_dimension':18,'stable_tangent_dimensions':17,'active_transverse_characteristic_polynomial':'s²+(1+kappa)s+kappa(1−Lambda_k/Lambda6)'},
'DERIVED_MOBILITY_NOISE.json':{'r_kappa':'kappa/(1+kappa)','kappa_1':{'r':0.5,'bracket':1.0},'kappa_4':{'r':0.8,'bracket':1.6},'scaled_candidate_drift':'r_kappa*(chi*lambda*z-z³/3)','scaled_candidate_diffusion':'sqrt(2*r_kappa)','scope':'generator coefficients; HM-004 separately pays process convergence'},
'COEFFICIENT_NEGATIVE_CONTROLS.json':{'reference_native_critical_diffusion':1.0,'kappa1_rate':0.5,'kappa4_rate':0.8,'same_stationary_quartic':True,'same_native_clock':False,'g0_rule':'define θ=0 categorical refresh, not singular g=0 Gaussian posterior'},
}),
'HM-004':dict(name='Stationary critical process',status='CONDITIONAL_STATIONARY_REFLECTED_PROCESS_THEOREM_SUBMITTED_FOR_INTAKE',brief='A stationary L1 martingale argument gives the quartic slow diffusion with mobility rκ for each fixed κ=1,4, assuming MT-003 polynomial bounds and explicitly stated boundary-layer probability decay. New-model iid starts are NOT handled.',files={
'STATIONARY_JOINT_LIMIT_OR_OPEN_LEMMA.md':'# HM-004 stationary reflected path-limit argument\n\nStart Π_D,N=π_D,Nρ_n and κ fixed. Let s=(κm+h)/(κ+1) and Z_N(u)=N^1/4 s(√N u). For the exact generator the drift is rκ[φ(θ)−h+(g_N/gc−1)m]. On the parity line φ=tanh h. Transverse θ has mean-square O(N^-1) under the stationary joint law, and symmetry kills the first transverse derivative at every h. Therefore φ=tanh h+O_L1(N^-1). The h−s gap is O_L2(N^-1/2). Taylor remainder after critical rescaling vanishes by the MT polynomial moments. The scaled bracket converges in stationary L1 to 2rκ; maximal jumps vanish. The MT local-Gibbs initial law converges to ρ_lambda. Standard martingale tightness and uniqueness then identify dZ=rκ(chi*lambda*Z−Z³/3)du+sqrt(2rκ)dB. Exact quantitative ledger and caveat in HM_THEOREM_PROOFS.md §HM-004. CONDITIONAL on N^3/4 π_D(boundary layer)→0 and MT uniform moments; the inherited review/current checkout could not be independently consulted in this run.\n',
'ERROR_AND_SCOPE.json':{'scaled_drift_L1_components':{'transverse_field':'O(N^(-1/4))','fifth_parity_Taylor':'O(N^(-1/2))','h_cube_to_s_cube':'O(N^(-1/4))','gain_times_m_to_s':'O(N^(-1/4))','reflected_boundary':'O(N^(3/4)*pi_D(B_D))'},'scaled_quadratic_variation_limit':'2*r_kappa','jump_size':'O(N^(-3/4))','initial':'stationary reflected Pi_D,N only','kappa':[1,4],'N':'multiples of 12','lambda':'fixed bounded','if_missing_boundary_payment':'OPEN_BOUNDARY_OR_MOMENT_PREMISE'},
'BOUNDARY_AND_FIELD_TAIL_PAYMENT.md':'# HM-004 boundary/field tail contract\n\nThe field θ remains unbounded: its full conditional Gaussian law is included exactly; stationary conditional moments control |θ_perp|² and |h−m|^q. The reflected jump generator equals the unrestricted one off the one-step departure layer. At the layer, a count-parity drift magnitude ≤2 bounds the correction, hence in the scaled drift the contribution ≤2N^(3/4)π_D(B_D) on fixed times. The MT fixed local free-energy barrier should pay this; the relevant current review file was unavailable and the boundary condition is kept explicit in the theorem. This proof does NOT transfer to diagnostic iid starts or unreflected new-model paths without new trapping/coupling work.\n',
}),
'HM-005':dict(name='Mediator nonidentifiability',status='EXACT_HIDDEN_CLOCK_NONIDENTIFIABILITY',brief='There are infinitely many reversible exact-native count-path lifts with the same Pi and distinct hidden refresh rate η. The field-clock witness separates them, while every count-path law is identical.',files={
'EXACT_LIFT_COMPARISON.md':'# HM-005 comparison of distinct auxiliaries\n\n(1) Ephemeral cavity ν_r+softmax: exactly reproduces single native LOO update; does not persist. (2) Reset lift: for each native count transition n→n′ redraw θ′∼ρ_n′, and optionally refresh θ∼ρ_n at arbitrary hidden rate η≥0. Joint Π_D reversible by native detailed balance; count process exactly native independent of η. (3) Persistent AUX-HB-OU-κ: same Π_D, but count projection has memory and different mobility. The alternatives show both *exact path nonidentifiability of hidden clock* (reset lifts) and *static nonidentifiability of kinetics* (persistent sampler). None selects a physical mediator.\n',
'HIDDEN_CLOCK_WITNESS.json':{'centered_field':'v_l=theta_l-(g/N)*(X^Tn)_l','event':'no count change during [0,t]','correlation_formula':'(g/N)*exp(-eta*t)*E_{pi_D}[exp(-d(n)*t)]','d(n)':'sum_{n_prime != n} q_D(n,n_prime)','count_path_law':'exact native for every eta','distinguishable_eta':True,'physical_mediator_selected':False},
}),
}

for stage,entry in stages.items():
    d=root/stage;d.mkdir(exist_ok=True)
    (d/'RESULTS.json').write_text(json.dumps({'task':stage,'title':entry['name'],'execution':'EXECUTED','science_status':entry['status'],'novel_HM_result':entry['brief'],'laws':{'original':'q12 labelled LOO HB FIN, unchanged','comparison':'AUX-HB-OU-kappa or separate exact-analysis lift','comparison_adopted_as_FIN':False},'scope':{'kappa':[1,4],'time':'model time','input_predecessor_status':'conditional','physical_evidence':'none'},'proof':'../HM_THEOREM_PROOFS.md','checks':'../verification.json','input_provenance':'../INPUT_PROVENANCE.json'},indent=2,ensure_ascii=False)+'\n')
    (d/'REPORT.md').write_text(f'# {stage}: {entry["name"]}\n\n**Scientific disposition:** `{entry["status"]}`.\n\n{entry["brief"]}\n\nSee the dedicated deliverables in this task, the complete derivations in [HM_THEOREM_PROOFS.md](../HM_THEOREM_PROOFS.md), and the source-access disclosures in [INPUT_PROVENANCE.json](../INPUT_PROVENANCE.json). The controls in `verify_hm.py` are bounded convention checks rather than a mechanical proof of the growing-N limits.\n')
    (d/'PROOF_OR_OPEN_LEMMA.md').write_text(f'# {stage}: proof, claims and exceptions\n\n{entry["brief"]}\n\n**Status:** `{entry["status"]}`.\n\nThe supporting equations, exact steps, and every scoped dependency are developed in [HM_THEOREM_PROOFS.md](../HM_THEOREM_PROOFS.md) under heading {stage}.\n\n'+('**Open independent-intake requirement:** confirm the exact MT-003 moment and boundary-layer hypotheses in the current review checkout; the diagnostic stationary critical law is conditional until that audit.\n' if stage=='HM-004' else '**Main constraint:** this does not grant a new fundamental FIN law or a physical mediator.\n'))
    (d/'NONCONCLUSIONS.md').write_text(f'# {stage} nonconclusions\n\n- No emergence of physical mediators, source of A7, real-time clock, conserved energy, or spatial locality.\n- No QM/GR/Standard Model/ToE inference.\n- A static Gaussian identity is not a persistent dynamic law.\n- Diagnostic κ and hidden η are external mathematical model parameters, not derived FIN constants.\n- No unreflected iid-start theorem for the diagnostic sampler.\n- No additional numerical precision claim beyond bounded finite computations.\n')
    (d/'REPLAY.md').write_text('Run `OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 ../verify_hm.py` from this task folder. Read `../HM_THEOREM_PROOFS.md` for analytical hypotheses; the script does not verify martingale convergence or MT premises.\n')
    for name,content in entry['files'].items():
        p=d/name
        if isinstance(content,dict):p.write_text(json.dumps(content,indent=2,ensure_ascii=False)+'\n')
        else:p.write_text(content)
    (d/'INPUTS.sha256').write_text('\n'.join(f"{q['sha256']}  {q['path']}" for q in provenance if 'sha256' in q and q['sha256'])+'\n')

(root/'CLAIM_REGISTER.json').write_text(json.dumps({k:{'status':e['status'],'scope':'finite-N exact or conditional stationary local limit; see task report','no_physical_promotion':True} for k,e in stages.items()},indent=2,ensure_ascii=False)+'\n')
(root/'HANDOFF.md').write_text('''# FIN HM-001–HM-005 — HM-only handoff

Date: 2026-10-16. This package contains only **new HM analysis**, never copies of historical prior packages, and makes **no edit** to existing AGENTS or primary FIN law. No Git checkout was claimed.

## Verdict

**ACCEPT_EXACT_CAVITY_BRIDGE_AND_PERSISTENT_PROJECTION_NO_GO; CONDITIONAL_STATIONARY_DIAGNOSTIC_LIMIT; AUXILIARY_KINETIC_SOURCE_UNDERDETERMINED.**

HM-001: a fresh cavity Gaussian latent draw, with the sum-of-exponentials tilted law, integrates to the original full-denominator LOO categorical distribution. This is an update identity, not a persistent medium.

HM-002: AUX-HB-OU-κ has the same joint static Π but the count compression is not a semigroup: `||QLPm||²=E Var(φ(θ)|n)>0`. The result is a rigorous non-Markov obstruction for **this** diagnostic, not every mediator.

HM-003: the joint generator has one soft and 17 stable coordinates for κ=1,4. Projection produces `r_κ=κ/(1+κ)`, specifically 1/2 and 4/5. This is genuine kinetic sensitivity with unchanged statics.

HM-004: the stationary reflected diagnostic admits a conditional critical process theorem: `dZ = rκ(χλ Z−Z³/3)ds+√(2rκ)dB` under the MT-003 polynomial moment theorem and explicitly displayed boundary rate. This is a new **model-theoretic** statement, not primary FIN dynamics. The current NT/MT review checkout was unavailable in this execution; independent review must confirm that the actual local-domain premises pay the bound. Neither iid initialization nor unreflected comparator dynamics is covered.

HM-005: an exact count-path reset lift exists for every independent hidden refresh rate η≥0; its count law is unchanged, whereas a hidden Gaussian residual correlation varies as `exp(−ηt)`. Hence the native path and auxiliary static measure do **not** source one unique hidden field clock.

## Novelty and evidence

The exact cavity integral, projection second-derivative defect, 18-coordinate coefficient derivation, conditional stationary path argument, and exact reset-lift witness are **HM-only** work. Inherited A7/MT critical premises are not discoveries here. The finite checks in verification.json test formulas and small examples; they are not certificate of MT hypotheses or path convergence. See INPUT_PROVENANCE.json for actually accessed predecessor sources and missing current review/MP7 originals.

## Outcome and next atom

**AUXILIARY_KINETIC_SOURCE_UNDERDETERMINED.** A target-independent, physically licensed source of a mediator clock, dynamical coupling, physical time and field observability is still absent. Do not postulate one under the primary FIN model ID; new PM requires separate human architectural approval. Stop here; do not initiate a larger kappa sweep, mediation hardware experiment, or QM/GR branch.
''')
(root/'REPORT.md').write_text('# HM campaign execution report\n\nFive stages executed as mathematics/proof gates. Exact finite identities in HM-001/002/003/005; conditional stationary local theorem submitted for independent MT audit in HM-004. No predecessor review/source files were modified. Outputs use English per plan. Run verify_hm.py for bounded tests; inspect INPUT_PROVENANCE.json for incomplete current-repository custody.\n')
(root/'NONCONCLUSIONS.md').write_text('# HM conclusions that are NOT licensed\n\n- Neither the static Hubbard–Stratonovich identity nor the diagnostic OU law derives a physical mediator.\n- The persistent AUX process is not primary FIN, even though it shares Π.\n- The diagnostic stationary theorem is not an iid-start or unreflected comparator theorem.\n- The exact reset lifts show hidden kinetic nonidentifiability, not the nonexistence of physical mediators.\n- No new source for A7, κ, g, physical clock, physical units, conservation, composition, spatial incidence, observable instruments or EV.\n- No quantum mechanics, gravity, SM, L_total or ToE closure; all legacy/selector/ontology/FIN-M0/B4/task335 boundaries remain unchanged.\n')
(root/'NEXT_OPEN_ATOM.md').write_text('# Next open atom after HM\n\n**HIDDEN_MEDIATOR_KINETIC_SOURCE_AND_PHYSICAL_INTERFACE.** Only if independently supplied, audit one typed object predicting the auxiliary field clock/coupling and a distinguishable observational consequence, without reading back the desired native law. No new object was supplied by HM. A narrow internal mathematical alternative is an independent check of HM-004 stationary boundary/MT assumptions; that review is **not** a new physical bridge. No automatic next campaign.\n')
(root/'PROPOSED_AGENTS_PATCH.md').write_text('''# Proposed ledger patch — NOT APPLIED

- 2026-10-16 HM: Exact cavity Gaussian refresh integrates to native LOO update, same denominator N and null attempts included; update-only auxiliary.
- Persistent AUX-HB-OU-κ (κ=1,4) has existing Π static invariant but fails exact autonomous count compression: second derivative obstruction `||QLPm||²>0`.
- Joint comparison generator supplies rκ=κ/(1+κ) and critical stationary comparison **conditional** on MT moment and boundary-layer hypotheses; current review checkout not independently audited.
- Exact reset lifts preserve native count paths for arbitrary hidden refresh η; static/path count data do not source hidden clock.
- Verdict AUXILIARY_KINETIC_SOURCE_UNDERDETERMINED; primary FIN unchanged; no physics promotion.
''')
(root/'REPLAY.md').write_text('# Reproduction\n\nFrom this directory: `OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 verify_hm.py`. For scientific verification, audit the exact proof sections in HM_THEOREM_PROOFS.md and compare the boundary/moment ledger with the actual current MT/NT reviews. Manifest: `sha256sum -c MANIFEST.sha256` from this directory. `package_check.py` also validates the expected manifest members. No past package is required to run small numerical checks, but independent source custody is explicitly limited.\n')
(root/'INPUTS.sha256').write_text('\n'.join(f"{q['sha256']}  {q['path']}" for q in provenance if 'sha256' in q and q['sha256'])+'\n')
(root/'package_check.py').write_text('''#!/usr/bin/env python3
from pathlib import Path
import hashlib
root=Path(__file__).parent
rows=(root/'MANIFEST.sha256').read_text().splitlines()
for row in rows:
    digest,name=row.split('  ',1)
    p=root/name
    assert p.is_file(), f'missing: {name}'
    assert hashlib.sha256(p.read_bytes()).hexdigest()==digest, f'hash mismatch {name}'
actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p.name not in ('MANIFEST.sha256','PACKAGE_VERIFICATION.json','__pycache__') and '__pycache__' not in p.parts}
expected={row.split('  ',1)[1] for row in rows}
assert actual==expected, ('uncovered files',actual-expected,expected-actual)
print(f'HM manifest: {len(rows)}/{len(rows)} valid with coverage')
''')
# no __pycache__; stale files ignored only if metadata, not in manifest.
files=sorted(p for p in root.rglob('*') if p.is_file() and p.name not in ('MANIFEST.sha256','PACKAGE_VERIFICATION.json') and '__pycache__' not in p.parts)
(root/'MANIFEST.sha256').write_text(''.join(f'{sha(p.read_bytes())}  {p.relative_to(root).as_posix()}\n' for p in files))
(root/'PACKAGE_VERIFICATION.json').write_text(json.dumps({'manifest_count':len(files),'excluded_from_manifest':['MANIFEST.sha256','PACKAGE_VERIFICATION.json'],'science_scope':'Exact HM algebra and conditional MT-based comparator theorem, no physical claims','generated_by':'finalize_hm.py'},indent=2)+'\n')
print('Manifest entries:',len(files))
