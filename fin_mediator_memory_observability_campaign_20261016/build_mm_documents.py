#!/usr/bin/env python3
"""Create MM-only documentation and hashes; never edits prior packages or AGENTS."""
from pathlib import Path
from hashlib import sha256
import json, zipfile
ROOT=Path(__file__).resolve().parent
DATA=Path('/mnt/data')
def save(rel, text):
 p=ROOT/rel;p.parent.mkdir(parents=True,exist_ok=True)
 p.write_text(text.strip()+'\n',encoding='utf-8')
def js(rel,obj):save(rel,json.dumps(obj,indent=2,ensure_ascii=False))
inputs={
 'MM supplied plan (attached pasted text)':DATA/'Wklejony Markdown(20261016-002418).md',
 'HM predecessor user handoff archive':DATA/'FIN_HM_001_HM_005_ONLY_HANDOFF_20261016.zip',
 'DKO predecessor user handoff archive':DATA/'FIN_DKO_001_DKO_005_ONLY_HANDOFF_20261015.zip',
 'KO earlier source handoff archive':DATA/'FIN_KO_001_KO_005_ONLY_HANDOFF_20261015.zip',
 'MT archived handoff (not full present review)':DATA/'files/FIN_MT_001_MT_005_FULL_HANDOFF_20261014.zip',
 'AGENTS library snapshot (NOT current Git checkout)':DATA/'files/AGENTS.md',
}
provenance={k:{'path':str(v),'sha256':sha256(v.read_bytes()).hexdigest(),'size':v.stat().st_size} for k,v in inputs.items() if v.exists()}
js('INPUT_PROVENANCE.json',{'actual_source_file_bytes':provenance,'unavailable_to_executor':['current checkout at 3be481a7 with active HM review folder','fin_auxiliary_mediator_review_20261016/PROOF_COMPLETIONS_AND_CORRECTIONS.md, POSTERIOR_RETENTION_AND_MM_SEED.md, verification.json as original bytes','fin_native_transfer_review/AGENTS live ledger','historical original uploaded wrapper custody'],'scope':'Do not claim proof/hashes for unavailable local review files. HM intake and MM plan were supplied as text by user.'})
lines=[f'{v["sha256"]}  {v["path"]}' for k,v in provenance.items()]
save('INPUTS.sha256','\n'.join(lines))

save('MM-001/OBSERVATION_CONTRACT.json', json.dumps({
 'model_ids':['FIN_NATIVE_A7_LOO_HB','AUX_HB_OU_KAPPA_1','AUX_HB_OU_KAPPA_4','EXACT_NATIVE_ETA_RESET_LIFTS'],
 'preparation':{'diagnostic':'stationary reflected Pi_D,N=pi_D,N rho_n','native_stationary':'reflected pi_D,N','native_iid':'late iid-to-local coupling only if inherited MT admits it'},
 'readout':{'fixed_block':'B_N chosen before observing labels','alpha':'k_N/N -> 1/2 (full alpha=1 separate negative control)','features':['f1(j)=sqrt(2)cos(2*pi*j/12)','f3(j)=sqrt(2)cos(6*pi*j/12)'],'times':[0,1,2,3,4],'observed':'V_r,N=sqrt(N)/k_N sum_{a in B_N} f_r(x_a(t))','binary_projection':'sign in {-1,0,+1} including zero'},
 'nuisance':{'allowed':'unknown fixed nonzero real gain separately for f1 and f3; one shared unit rescaling of clock','excluded':['offsets','time-varying sign flips','erasure','additive/coloured detector noise','unknown labels','resampled block','arbitrary calibration']},
 'reference_rule':{'rho_rj':'sin(pi*mean(sign(V_r(0))*sign(V_r(j)))/2)','q':'rho_null,2/rho_null,1 with nonzero guard and rho_null,3=rho_null,1*q^2 check','if_failed':'UNCLASSIFIED / REFERENCE_UNSUPPORTED'},
 'horizon':'fixed fast model time [0,4]; no growing M or shrinking alpha','status':'frozen before new MM certificate; existing supplied plan'
 },indent=2,ensure_ascii=False))

save('MM-001/MODEL_AND_MEMORY_CONTROL_MAP.md',r'''
# MM-001 — Actual experiment and control models

**Frozen readout:** five simultaneous observations of the f1 (null) and f3 (active) q12 Fourier harmonics, at times 0,1,2,3,4, from *one fixed predetermined* labelled block B_N, k_N/N→1/2. The observer is allowed only the ten signs sign(V_r,N(t)), taking sign(0)=0. Neither theta, population F, private U, complete label lists, nor a hidden terminal value is part of that observation. The analysis may use the full two-block labelled generator. For native FIN the original LOO target uses n-e_i, denominator N, unit attempts. For the diagnostic the same selected slots refresh from softmax(X theta), while theta obeys its separately declared OU feedback SDE with kappa=1 or 4.

**Preparations:** both processes are compared in their respective stationary *reflected local sectors*: native pi_D,N (labelled lift), diagnostic Pi_D,N=pi_D,N rho_n. Only the native late-iid process may inherit the already-admitted MT coupling; an analogous diagnostic iid theorem has not been supplied. The critical soft coordinate remains a random constant over each fixed fast-time interval.

**False-positive control:** the native partially observed f3 readout already satisfies V=F+U. Its limiting C(t)=gamma^-1 e^-gamma|t|+c_alpha e^-|t| with c_alpha=(1-alpha)/alpha, producing TWO distinct time scales whenever alpha<1 and gamma!=1. Its ordinary two-pole Hankel curvature C1*C3-C2^2 is strictly positive. Thus generic curvature/two-pole plots cannot be attributed to feedback. Native full readout alpha=1 is an additional control. In exact-native HM-005 eta reset lifts, every count-path observable (including these signed blocks) has the same native law at every eta. A different hidden clock without changed data is therefore possible.

**Boundary of nuisance group:** fixed nonzero real gains on either observed channel preserve pairwise products of its signs (even negative gain changes both signs); a common clock-unit change is calibrated internally by null-mode decay. Separate time-varying sign changes, offsets, unknown channel labels, erasures and coloured readout noise are not covered. The diagnostic kappa is a supplied comparison constant, not natural time.

**No physical promotion:** a negative filtered determinant distinguishes only the named ideal mathematical process classes, not a unique physical medium or an apparatus. Older generic multi-time identities and the DKO native memory theorem are reused, not rediscoveries.
''')

save('MM-001/PRIOR_WORK_AND_NUISANCE_SCOPE.md',r'''
# Reused rather than new

HM-001–005 supply the exact static MP7 joint density, persistent feedback generator and stationary critical limit; MM does not resubmit them. DKO-001/002 supply native fixed-alpha joint F/U fast-path theorem. KO supplies half-block initial hypergeometric sampling and MT the local stationary Gibbs/moment/boundary control. Old report 298 concerns generic multi-time semigroup fingerprints; neither generic Hankel algebra nor the existence of observation-induced memory is new MM work.

New MM scope: fixed-horizon *diagnostic* three-component coupled (F,H,U) process, a *reference-pole-filtered sign-law functional*, its native full/partial/eta controls and bounded-observation transfer. No fresh physical units, apparatus, independent mediator source or finite-data calibration is admitted.
''')

save('MM-002/ACTUAL_FAST_CHARACTERISTICS.md',r'''
# MM-002 — Generator-derived fast characteristics

Write alpha=k/N and f_r(j)=sqrt(2)cos(2*pi*r*j/12), r=1,3. For a predetermined block B define a_a=N/k-1 if a in B, and -1 otherwise. Then F_r=N^-1/2 sum_a f_r(x_a), U_r=N^-1/2 sum_a a_a f_r(x_a), V_r=F_r+U_r. For the AUX-HB-OU labelled generator, q_j(theta)=softmax_j(X theta), whose proposal does *not* depend on the refreshed slot's old label. On any unrestricted count configuration (or reflected interior),

    L F_r = sqrt(N) sum_j f_r(j) q_j(theta) - F_r,
    L U_r = -U_r,                           (EXACT)

because sum_a a_a=0. Reflection rejects whole proposals whose total n leaves D_N; it does NOT renormalize q. All boundary terms must be retained as O(sqrt(N)*Pi_D(boundary)) in stationary L1 (fixed alpha), assuming the inherited exponential local-Gibbs departure bound.

Take an orthonormal Fourier X representation of A7: in the active cosine r=3 direction X_{j,3}=sqrt(Lambda3/12) f3(j). Let c=g_c*Lambda3/12=Lambda3/Lambda6 in (0,1), gamma=1-c, and H_3=sqrt(N)*sqrt(Lambda3/12)*theta_3. The exact field drift is

    L H_3 = kappa (c_N F_3-H_3),    c_N=g_N Lambda3/12,
    d<M_H3>/dt = 2 kappa c_N.

At critical lambda=0, c_N=c. Cross jump/continuous brackets between H and F,U vanish identically. In the null r=1 direction A7 has Lambda1=0, so no H1 drives F1.

The exact predictable count-jump brackets are, with w_a(F)=1 and w_a(U)=a_a,

    d<M_{X_r},M_{Y_s}>/dt
      =N^-1 sum_a w_a(X)w_a(Y) sum_j q_j(theta)
        [f_r(j)-f_r(x_a)][f_s(j)-f_s(x_a)].

Here X,Y are F or U. For r=s the inner square is E_q f_r(new)^2 + f_r(old)^2 - 2 f_r(old) E_q f_r(new), NOT identically 2 statewise. The N=12 independent jump fixture with block indices all even, uniform counts and theta=0 gives FF=2, UU=2 and FU=1: mixed brackets vanish **only in the stationary/conditional-block limit**, not as exact finite-state identities.

On the admitted MT local positive cylinder, the parity coordinate has theta6=O_Lp(N^-1/4) while all transverse field coordinates have O_Lp(N^-1/2). Parity-axis symmetry implies the f1/f3 mean proposals vanish at every pure soft configuration. A Taylor expansion around that axis, bounded softmax derivatives, and stationary Hölder estimates give for fixed T

    E_Pi integral_0^T |sqrt(N) f3·q(theta)-H3| dt
       <= C_T [N^-1/4 + N^-1/2 + |c_N-c|],

and the corresponding f1 term tends to 0. Conditional hypergeometric block moments and full stationary local moments yield, in integrated L1,

    <M_F3> -> 2t,    <M_U3> -> 2 c_alpha t,
    <M_F3,M_U3> -> 0,    <M_H3> -> 2 kappa c t,

with zero cross-mode limits. The same local estimates pay bracket errors O(N^-1/4)+O(N^-1/2) (not a claimed universal finite-N bound). Gaussian theta tails are controlled by its known conditional normal rho_n and stationarity, not dropped. For each fixed alpha away from zero, jumps of F/U are O(N^-1/2); the OU field has continuous Brownian martingale part. These are generator calculations plus stated inherited stationary moment hypotheses, not a guessed two-OU covariance.
''')

save('MM-002/STATIONARY_INITIAL_LAW.md',r'''
# MM-002 — Stationary joint initialization

The exact diagnostic stationary local joint measure is Pi_D(n,dtheta)=pi_D(n)rho_n(dtheta), with rho_n=N((g_N/N)X^T n,(g_N/N)I_7). Conditional on n, the chosen block is a hypergeometric sample (labelled permutation symmetry) and theta is independent of its selection. Under the admitted MT global-fast Gaussian limit, F3 ⇒ N(0,1/gamma), independently of the quartic soft limit, while

    H3 | F3 ~ c F3 + sqrt(c) G,  G~N(0,1),

at criticality, with G independent of F3. Hence (F3,H3) is jointly Gaussian with covariance gamma^-1 [[1,c],[c,c]]. Fixed-fraction KO hypergeometric CLT gives U3 ⇒ N(0,c_alpha), independently of the joint (F3,H3) limit, through conditional characteristic functions. The null F1 and U1 have variance 1 and c_alpha and are independent of the active pair in the limit. The soft variable is random but asymptotically frozen over model fast time O(1).

Finite-N independence does NOT hold and zero covariance alone does not prove the Gaussian factorization; both depend on the inherited joint Gibbs and conditional hypergeometric limit. The diagnostic is initialized exactly in Pi_D, not from an unproved new-model iid mixing window.
''')

save('MM-002/JOINT_FAST_LIMIT_OR_OPEN_LEMMA.md',r'''
# MM-002 — Conditional diagnostic fast-path theorem submitted for audit

**Claim (analytical, dependent on inherited MT/KO bounds):** For each kappa∈{1,4}, fixed alpha∈(0,1), fixed T<infinity, N divisible by 12, g_N=g_c at lambda=0, and stationary reflected Pi_D,N preparation, the joint soft/global/field/private observation path converges on D([0,T]) to a random frozen quartic soft coordinate plus the stationary continuous Gaussian Markov system

    dF3=(-F3+H3)dt+sqrt(2)dB,
    dH3=kappa(c F3-H3)dt+sqrt(2 kappa c)dW,
    dU3=-U3 dt+sqrt(2 c_alpha)dV,

with (F3,H3)(0) covariance gamma^-1[[1,c],[c,c]], U3(0)~N(0,c_alpha), independent initial blocks/drivers, and V3=F3+U3. Null f1: independent rate-one stationary F1 and U1 with variances 1 and c_alpha, so V1 has covariance (1+c_alpha)e^-|t|. Brownian symbols represent emergent limiting martingales; they are not added to native FIN.

**Payment route:** (1) exact characteristics and stationary initial laws are in adjacent files; (2) stationary moment+Gaussian conditional tails imply integrated drift error ->0, including the parity/fast cross terms and unbounded theta; (3) conditional hypergeometric concentration and MT moments make integrated mixed bracket errors vanish; (4) reflected boundary drift/brackets are bounded by a polynomial in N times the inherited exponentially small boundary layer; (5) jump sizes vanish O(N^-1/2), and Doob/BDG on stopped semimartingales with local high moments yields maximum tightness and square-martingale UI; (6) limiting continuous martingale problem has unique law (linear Lipschitz drift, positive constant covariance), so cross-zero *joint martingale* brackets and Gaussian initial law yield independent limiting drivers. Remove localization by stationary moment/tail bounds.

**Unpaid source verification:** This execution does not possess current HM/MT review completions as original files, only the user-supplied HM intake assertions and archived original HM/MT campaign outputs. Thus inherited boundary rate, conditional joint Gaussian MT theorem, and the complete tightness hypotheses were not independently re-audited against the active Git checkout. The preceding steps form a conditional analytical proof submission, not a mechanically verified theorem. If any inherited high-moment/local-domain hypothesis fails for this stationary comparison, downgrade this process limit to OPEN_STATIONARY_MOMENT_OR_BOUNDARY_GATE and do not promote MM-004.

**Prohibited upgrades:** no iid-start diagnostic transfer, unreflected comparator limit, shrinking alpha_N, increasing observation horizon, finite-N rate/N0, physical field/clock or source claim.
''')

save('MM-002/BOUNDARY_UI_PAYMENT.md',r'''
# MM-002 — Boundary, tails, and UI ledger

The diagnostic's invariant Pi_D has marginal pi_D in counts; its field conditional is an unbounded Gaussian rho_n. All moments of theta are controlled by E_piD ||(g/N)X^T n||^p plus Gaussian variance powers. In the local MT sector, the soft theta6 is O_Lp(N^-1/4) and transverse theta is O_Lp(N^-1/2). Reflection removes only count proposals leaving D_N. For fixed alpha the absolute drift of F,U due to suppressed proposals is ≤C sqrt(N) 1_{n∈B_D}, and bracket differences are ≤C 1_{n∈B_D}. The inherited exponential departure estimate pays their integrated stationary expectations on fixed T (and stronger polynomial horizons where independently proved). Field OU continues unbounded with no reflection. Maximal semimartingale tightness follows from stationary p-th moments, BDG, bounded F/U jumps, OU fourth moments and localization; square-martingale UI requires fourth-moment bounds, not merely convergence of first brackets. Native actual late iid-path transfer may use its own MT coupling, but NO diagnostic iid analogue is available. This ledger is conditional on exactly the inherited MT high moments and positive local cylinder.
''')

save('MM-003/FILTERED_CONSTRAINT_PROOF.md',r'''
# MM-003 — Reference-pole-filtered constraint (exact Gaussian algebra)

Let c=Lambda3/Lambda6, gamma=1-c, 0<c<1. The diagnostic global fast drift matrix is J=[[-1,1],[kappa*c,-kappa]], the diffusion bracket diag(2,2*kappa*c). The Lyapunov equation J Sigma+Sigma J^T+Q=0 has stationary Sigma=gamma^-1[[1,c],[c,c]], as independently checked by exact symbolic simplification. The eigenvalues of -J are

    beta_s,f = [(1+kappa) ∓ sqrt((1-kappa)^2+4*kappa*c)]/2.

Evaluating the characteristic polynomial at gamma and 1 shows beta_s<gamma<1<beta_f for positive c and kappa. Stationary F covariance is gamma^-1[w_s exp(-beta_s|t|)+w_f exp(-beta_f|t|)], with w_s=(beta_f-gamma)/(beta_f-beta_s), w_f=1-w_s, both >0. The independently sampled private U contributes c_alpha exp(-|t|), c_alpha=(1-alpha)/alpha. Thus the partial diagnostic readout has THREE component rates (two feedback poles and the private pole), whereas the native readout has the native global rate gamma plus the same private rate1. Full alpha=1 omits the private term.

A common positive model-time rescaling t->nu t changes x=e^-beta_s nu h, y=e^-beta_f nu h, q=e^-nu h without changing the ordering x>q>y. For active unnormalized C_j=A x^j+B y^j+D q^j, define d_j=C_{j+1}-q C_j for j=1,2,3. Algebra gives

    H_ref=d1*d3-d2^2=A*B*x*y*(x-q)*(y-q)*(x-y)^2 < 0.

The coefficient D cancels exactly, regardless of observation fraction. If C_j=A x^j+D q^j (native partial) or C_j=A x^j (native full), H_ref is identically zero. In contrast ordinary C1*C3-C2^2 for native partial equals A*D*x*q*(x-q)^2>0 and is a false-positive memory indicator.

**Bounded-sign law:** At each mode r and lag j=1..4, let s_rj=E[sign V_r(0) sign V_r(j)], with sign(0)=0. For a nonsingular Gaussian pair s=(2/pi)arcsin(rho), so R_rj=sin(pi s_rj/2)=Corr(V_r(0),V_r(j)). Set q=R_null,2/R_null,1 only when R_null,1≠0 and the null r=1 one-rate checks are supported, including R_null,3=R_null,1*q^2. Replace C_j by the normalized R_active,j and evaluate H_ref with this q. Because normalization divides all C_j by the same C_0>0, the sign of H_ref is unchanged and the expression is a functional of bounded ternary observations *only*. It is not an estimate from one five-time observed tuple.

**Gauge restrictions:** Each observed channel may have its own unknown FIXED nonzero real multiplier: sign products cancel its sign. A single common change of model-time units changes q and active rates consistently and is removed by null calibration. No offset, time-dependent sign-changing gain, detector noise, label reassignment or independent channel clock is included. Under an unsupported null reference, the rule returns REFERENCE_UNSUPPORTED, never refits q.

**Margin and scope:** At alpha=1/2,h=1,kappa1/4, nominal H_ref is strictly negative. The adjacent interval file uses mpmath.iv with an expressly ASSUMED input gamma interval [0.16257283698729,0.16257283698731] inferred from the prior nominal spectrum; source spectral enclosure must be independently bound for a certified physical-parameter result. No uniform margin for alpha->0 or q->0,1 has been proved.
''')

ver=json.loads((ROOT/'verification.json').read_text())
js('MM-003/BOUNDED_LAW_FUNCTIONAL.json',{'observed':'10 ternary signs at five times; seven pair mean functionals', 'null_mode_pair_lags':[1,2,3], 'active_pair_lags':[1,2,3,4], 'from_ensemble_signs':'R_rj=sin(pi*E(sign V0 sign Vj)/2)','q_guard':'R_null_1!=0; null_2/null_1 in (0,1) and null_3==null_1*q^2 in ideal law', 'determinant':'D_j=R_active,j+1-q R_active,j for j=1..3; H=D1 D3-D2^2', 'native':'H=0 exactly in limiting native partial and full Gaussian family','AUX_HB_OU':{'kappa_1':ver['working_precision']['kappa_1']['filtered_H_normalized'],'kappa_4':ver['working_precision']['kappa_4']['filtered_H_normalized']},'invalid_reference':'UNCLASSIFIED/REFERENCE_UNSUPPORTED','scope':'ensemble law of frozen stationary models; no one-record estimator or physical noise model'})
js('MM-003/IDEAL_MARGIN_CERTIFICATE.json',{'result':'STRICTLY_NEGATIVE_CONDITIONALLY_ON_GAMMA_INPUT_BOX','gamma_interval_assumption':['0.16257283698729','0.16257283698731'],'kappa_1':ver['interval_checks']['kappa_1'],'kappa_4':ver['interval_checks']['kappa_4'],'not_claimed':'This execution did not independently certify the A7 spectral input enclosure as original bytes; no finite-N margin or physical calibration is covered.'})

save('MM-004/FINITE_PARENT_LAW_TRANSFER.md',r'''
# MM-004 — From joint fast-path convergence to a frozen bounded sign law

**Dependent on MM-002 joint process theorem.** For fixed alpha=1/2, h=1 and T=4, let E_N be the full stationary labelled path law of either the original reflected native or the stationary reflected diagnostic. Observe only (sign V1(t_i),sign V3(t_i))_{i=0}^4 in {-1,0,+1}^{10}. Under the joint fixed-horizon weak path limit, every finite coordinate vector converges to its continuous Gaussian counterpart. Each coordinate has strictly positive limiting variance; its hyperplane {V=0} has probability zero. Therefore for each sign vector its cylinder is a continuity set; the complete finite-alphabet joint law converges, including finite-N ZERO outcomes whose total probability vanishes. No survivor postselection occurs.

Each of the seven pair means (null at lags 1,2,3 and active at lags 1,2,3,4) is a bounded expectation; convergence follows from that finite-alphabet law. The sine transform is continuous. The null denominator tends to exp(-1)>0 and the null single-rate consistency identities hold. Consequently the guarded q estimate tends to exp(-1), and the polynomial H_ref functional tends to 0 for native reflected (alpha=1/2 and alpha=1) and to a strictly negative limit for each fixed diagnostic kappa=1,4. Thus **for the two specified stationary sequences, their population bounded-observation laws have a nonzero asymptotic model contrast**. This is a model-law discrimination statement, not a finite-data classifier with guaranteed error.

Only native FIN may separately inherit bounded finite-window readout transfer from iid initialization after the admitted deterministic MT preparation. The diagnostic is initialized in Pi_D only: neither iid preparation nor unreflected diagnostic coupling is claimed. The old HM-004 critical slow theorem does not provide this fixed-fast readout result; it is a separate MM-002 gate.

**No finite-N error bars:** weak bounded-law convergence does not yield a useful N0, a rate, a sample-size guarantee, or total-variation path equivalence. The null reference is an ideal population-law check. A finite empirical estimate of q and H_ref requires new sampling/error/dependence bounds that are not executed here.
''')
js('MM-004/FROZEN_EXPERIMENT.json', {'models':['FIN_NATIVE_A7_LOO_HB','AUX_HB_OU_KAPPA_1','AUX_HB_OU_KAPPA_4'], 'observation':['mode_1_sign','mode_3_sign'],'sample_times':[0,1,2,3,4],'alpha':0.5,'negative_control_alpha':1,'zero_outcome_retained':True,'finite_alphabet_size':3**10,'statistics':'7 pair expectations, not all 3^10 probabilities required','preparations':{'native':'local stationary, plus only inherited late iid bounded-path transfer','diagnostic':'stationary reflected Pi_D only'},'theorem':'conditional law-level limits; no finite repeated-trial testing claim'})
js('MM-004/ERROR_AND_SCOPE_LEDGER.json',{'conditional_on':['HM accepted Pi_D static joint law and fixed-kappa diagnostics','KO fixed fraction hypergeometric CLT and MT local Gibbs high moments','MT exponential departure control','MM-002 new fixed-fast joint martingale path argument'],'paid_analytically':['zero-sign boundary has zero limiting Gaussian mass','bounded finite-alphabet weak transfer','continuous sine and guarded ratio','reference pole monotonicity and negative factorization'],'not_paid':['numerical finite-N path error','external detector noise model','iid-start diagnostic','growing observation time','physical time','finite-data statistical error or power']})

save('MM-004/NEGATIVE_CONTROLS.md',r'''
# MM-004 — Essential negative controls

1. **Native half-population:** two relaxation poles and ordinary positive Hankel curvature appear without any persistent mediator, but the reference-filtered determinant is exactly 0 in its limiting Gaussian law.
2. **Native full population:** one active pole; reference-filtered determinant remains 0.
3. **Exact-native hidden reset lifts at arbitrary eta:** every count path is exactly the native path, hence all block sign records and filtered laws are identical for all eta. Hidden clock nonidentifiability is exact.
4. **AUX-HB-OU kappa=1 and kappa=4:** the diagnostic law shares static Gibbs marginal but has two global feedback poles; after removing the private reference pole, both generate strictly negative H_ref in the ideal stationary limit.
5. **Reference failure:** absent a valid null-mode exponential, a zero or invalid denominator yields REFERENCE_UNSUPPORTED, not a patched contrast.
6. **Generic coloured readout memory:** may create the same extra poles and invalidate mechanism attribution; no noise or apparatus class was supplied to exclude that possibility. Do not claim unique physical mediation.
''')

save('MM-005/SPECIFICITY_VERDICT.md',r'''
# MM-005 — VISIBLE_FEEDBACK_CONTRAST_IN_DECLARED_CLASS (conditional)

The signed reference-filtered Gaussian-law functional is 0 for both native full/partial count observation controls and for every exact-native eta hidden reset lift, but negative for the two specified AUX-HB-OU-kappa stationary persistent-feedback diagnostics. The new fast-path argument derives the diagnostic limit from the existing mixed generator rather than reading it off from the static Gibbs measure or slow critical HM theorem. Its finite-model readout consequence follows conditionally by bounded path weak convergence with zero outcomes retained.

**What it actually identifies:** within the fixed model family, block, feature maps, preparation, lag, null-reference clock and perfect readout contract, it detects a failure of the native partial-observation Gaussian mechanism. It does not recover a unique kappa from this one scalar statistic with unknown alpha; it does not identify any hidden physical medium or rule out broader coloured observation laws. The HM-005 exact eta lifts show genuine hidden-clock non-identifiability from *all count paths*. The primary FIN law is untouched.

**Next specific physical prerequisite:** a separately licensed actual measurement/preparation/readout mapping including feature assignment, a valid null clock reference and bounded channel memory/errors. Without that, no empirical or natural mediator claim is activated. Return to architectural review; no further kappa, N, gain or time sweep.
''')

save('NONCONCLUSIONS.md',r'''
# MM — scientific nonconclusions

- AUX-HB-OU-kappa is an explicitly supplied diagnostic kinetic model, not a new FIN law or a sourced mediator.
- Shared Gibbs equilibrium does not establish shared count paths or the origin of the hidden OU clock.
- The two-pole native readout itself has memory; a generic nonexponential correlation is not mediator evidence.
- The negative filtered determinant distinguishes only the specified limiting statistical laws, not all hidden mechanisms or coloured detector-memory alternatives.
- No observed theta, direct field measurement, physical wave propagation, locality, fundamental particles, quantum instrumentation, time/energy units, gravitational dynamics or ToE.
- No iid diagnostic preparation theorem, unreflected diagnostic model result, alpha_N shrinking, T_N increasing, useful N0, path-TV equivalence, finite-data estimator accuracy, numerical model-selection power or detector calibration.
- The independent numerical certificate is conditional on a declared gamma3 input bracket and is not an independent exact A7 spectral certification.
- Historical reviews/AGENTS Git HEAD and missing archived wrappers are not claimed to have been inspected or verified; source predecessor statements are explicitly imported.
''')

save('NEXT_OPEN_ATOM.md',r'''
# MM next atom — a physically licensed observation-memory and null-clock reference

**Unpaid input:** concrete, independent source/preparation/readout evidence specifying the mapping from labelled FIN counts to simultaneous f1/f3 readouts, with a clock reference and an upper class for readout memory/error. The exact MM mathematical contrast is conditional on its ideal law; coloured instruments can mimic extra poles. Without a genuine typed input return OPEN_PHYSICAL_INTERFACE, not a new simulation or kappa sweep. Optional separate mathematical work: a prospectively frozen finite-data dependent-sign estimator with confidence/size bounds for one specified finite N and preparation, which would *not* by itself license physics.
''')

save('REPLAY.md',r'''
# MM-only reproducibility

At the directory of this MM campaign, run `PYTHONDONTWRITEBYTECODE=1 python3 verify_mm.py`, and check `verification.json`. The verifier uses `sympy`, `mpmath.iv`, and `numpy`; it performs exact polynomial/Lyapunov tests, fixed-state labelled jump sums, nominal spectral Gaussian calculations and outward arithmetic conditional on gamma3 lying in the explicitly declared interval. It does NOT simulate trajectories or prove the asymptotic martingale limit.

Run `python3 verify_package.py` to check manifest content coverage and hashes. No program rewrites AGENTS, the HM/MT/KO/DKO source archives, or repository state. `build_mm_documents.py` regenerates narrative scaffolding; its output may be reformatted but the scientific assumptions and checks are in the files. This file is not a script for reconstructing missing original uploaded wrappers or reviews.
''')

save('PROPOSED_AGENTS_PATCH.md',r'''
# Proposal only — not automatically applied

MM-001–005 (2026-10-16): conditional fixed-alpha stationary model-theory result. Actual AUX-HB-OU-kappa joint F/H/U generator and Gaussian fast limit submitted under inherited HM/MT/KO static/joint moments and local boundary hypotheses; a seven-pair bounded sign-law reference-filtered determinant is identically zero in native full/partial and exact-native eta-lift controls, strictly negative in diagnostic kappa1/4 limits. Conditional mpmath.iv intervals and 35 bounded tests, not full mechanized proofs or numerical finite-N rates. Verdict VISIBLE_FEEDBACK_CONTRAST_IN_DECLARED_CLASS; source of mediator, instrument memory and physical time remain OPEN. Before admission review the current local HM/MT proofs and MM-002 path UI/boundary payments; do not adopt AUX as FIN.
''')

save('DEPENDENCY_DAG.md',r'''
# New MM dependency DAG

```
existing A7 labelled FIN / HM AUX-HB-OU kappa diagnostic [two different supplied laws]
  -> HM static Pi_D & MT high moments, fixed critical stationary sector [inherited]
  -> MM-001 fixed block+features+five-time count-only contract [new freeze]
  -> MM-002 exact global/field/private characteristics; initial conditional Gaussian+HG CLT [new diagnostic analysis]
  -> MM-002 stationary fast joint martingale limit [new conditional proof gate]
  -> MM-003 two feedback poles + private pole and exact negative filtered determinant [new algebra]
  -> MM-004 bounded ternary joint law and 7 sign-pair functionals [new conditional transfer]
  -> MM-005 specificity within declared class / physical source and instrument premise remain open
```

Reused DKO native partial two-pole, KO sampling and HM-005 exact eta lifts are essential CONTROLS, not new MM theorems. M-level statistical result is not PM/PB/PP/EV physics evidence.
''')

save('HANDOFF.md',r'''
# FIN MM-001–MM-005 — MM-only research handoff

Date: 2026-10-16. Primary law: unchanged q12 A7 labelled unit-attempt LOO-HB. Comparator: separately declared stationary reflected AUX-HB-OU-kappa, kappa=1/4. No previous campaign documents copied or AGENTS modified.

## Verdict

**VISIBLE_FEEDBACK_CONTRAST_IN_DECLARED_CLASS — CONDITIONAL_ANALYTICAL_PROCESS_BRIDGE_SUBMITTED.**

MM-001 freezes ten count-only signs at times0..4 for null f1/active f3 on one predetermined half-population block. Native half-block observation already has two decay rates: generic memory is therefore a false positive. Native full, native half and exact-native eta-reset lifts are the controls.

MM-002 derives the diagnostic global/field/private fast dynamics from the labelled jump+OU generator rather than from matching one-time equilibrium: dF3=(-F3+H3)dt+sqrt(2)dB, dH3=kappa(cF3-H3)dt+sqrt(2 kappa c)dW, dU3=-U3dt+sqrt(2c_alpha)dV, with c=Lambda3/Lambda6 and c_alpha=(1-alpha)/alpha. Conditional static Gibbs/field normal and block hypergeometric laws give the Gaussian initial covariance; inherited MT moments, boundary estimates and Gaussian tails support integrated drift/bracket and martingale tightness arguments. This is a new **conditional** fast-process proof: active HM/MT review bytes and its strongest tightness/boundary assumptions were not independently available in this execution. No diagnostic iid-start or finite-N rate is claimed.

MM-003 filters the known null/private rate-one pole from normalized active correlations reconstructed using bounded sign-pair means. The resulting determinant is identically 0 for native partial/full and strictly negative for kappa1 and kappa4, by exact symbolic factorization with positive feedback weights. Computed working values at alpha=1/2,h=1 are approximately -0.0002940687670 and -0.00000416162704, respectively. Outward interval calculations establish their signs *conditional on a stated gamma3 bracket* inferred from inherited nominal spectral inputs; do not cite them as fresh exact A7 spectral certificates.

MM-004 states a finite-model asymptotic count-only contrast conditional on MM-002. Each 10-sign joint law, including finite zero outcomes, converges to a Gaussian sign law on the fixed horizon. The seven pair expectations therefore determine the ideal zero/negative functional. It is an **ensemble-law** theorem, not a single-shot decision or finite-sample calibration theorem.

MM-005 closes with VISIBLE_FEEDBACK_CONTRAST_IN_DECLARED_CLASS. Coloured instrumentation, an unsupplied preparation/readout/clock mapping, other hidden kinetics and natural mediator identity remain open. Exact-native eta reset lifts still render hidden clocks invisible to *all* count observations. No ToE/QM/GR promotion.

## Evidence and status

See MM-002/JOINT_FAST_LIMIT_OR_OPEN_LEMMA.md for the analytic proof and its exact imported hypotheses, MM-003/FILTERED_CONSTRAINT_PROOF.md for exact Gaussian algebra, MM-004/FINITE_PARENT_LAW_TRANSFER.md for the bounded-observation limit, verification.json for bounded calculations, INPUT_PROVENANCE.json for actual files accessible, and MANIFEST.sha256 for package integrity. `verify_mm.py` checks finite formulas; it does not independently certify the preceding stochastic theorems. No Git checkout at 3be481a7 or current local review audit is claimed.
''')

save('REPORT.md',r'''
# MM research report

Scientific question: can a count-only partial-observation record distinguish persistent AUX-HB-OU feedback from native FIN, after removing memory already created by partial sampling?

**Result:** Yes, in the fixed-alpha stationary *conditional asymptotic mathematical model* with a valid null reference: the derived Gaussian feedback process has a strictly negative filtered four-lag determinant while the native full/partial and exact-native hidden reset controls give zero. Bounded finite-record-law transfer is conditional on the separately derived full diagnostic fixed-fast-horizon martingale limit. Results are not a physical mediator discovery, not a finite-data statistical test and not accepted as peer-reviewed proof.

The scientific novelty is the actual diagnostic generator-to-fast-readout bridge and essential false-positive control, not the old generic two-exponential Hankel factor. See proof and numerical scope ledgers. Next gate is actual independent observation/source evidence, not extra numerical precision.
''')

save('RESEARCH_SUMMARY_PL.md',r'''
# FIN MM — podsumowanie badań

Porównano oryginalny FIN z wcześniej zdefiniowanym, osobnym modelem diagnostycznym, w którym pole Gaussa jest trwałe i sprzężone z licznikami. Obserwator dostaje tylko znaki dwóch składowych odczytu części populacji w pięciu chwilach. Natywny FIN może mieć pamięć obserwacyjną nawet bez trwałego pola; dlatego zwykłe dwie skale zaniku nie są dowodem mediatora.

Nowy test filtruje składnik o znanym tempie zaniku za pomocą referencyjnego trybu zerowego. W warunkowym modelu granicznym determinant jest dokładnie 0 dla natywnego pełnego/częściowego odczytu, a ujemny dla porównawczego trwałego pola z kappa=1 i 4. Wartości nominalne to -0,0002940687670 oraz -0,00000416162704. Liczby dotyczą modelu probabilistycznego, nie danych fizycznych.

Wyprowadzono także wspólną dynamikę szybką F/H/U z generatora diagnostycznego, kontrolując warunki początkowe, mieszane nawiasy martyngałów i brzeg odbijany. Dowód graniczny jest warunkowy wobec wcześniej zaakceptowanych oszacowań MT i wymaga niezależnego audytu pełnego lokalnego repo. Zbieżność rozkładu dziesięciu znaków i siedmiu odpowiednich korelacji wynika z tego twierdzenia oraz zaniku masy na granicach znaków.

Wniosek: matematycznie rozróżniamy dwie wymienione klasy dynamiki z obserwacji liczników, lecz nie dowodzimy istnienia ani unikatowości mediatora fizycznego. Inne źródła pamięci odczytu mogą imitować sygnał. Modele HM z ukrytym zegarem eta nadal zachowują wszystkie trajektorie liczników FIN bez możliwości rozróżnienia eta.
''')

save('CLAIM_REGISTER.json',json.dumps({'MM-001':{'status':'EXACT_READOUT_CONTRACT','new':'fixed comparison and explicit false-positive controls'},'MM-002':{'status':'CONDITIONAL_ANALYTICAL_FAST_PATH_THEOREM_SUBMITTED','hypotheses':'imported local MT moments, boundary, joint KO Gaussian and stationary diagnostic Pi_D; not independently audit-verified'},'MM-003':{'status':'EXACT_GAUSSIAN_FACTOR_WITH_INPUT_CONDITIONAL_INTERVAL','result':'native H=0; diagnostic H<0 for kappa=1,4'},'MM-004':{'status':'CONDITIONAL_FINITE_PARENT_BOUNDED_SIGN_LAW_LIMIT','requires':'MM-002 joint fast path and zero-boundary continuity'},'MM-005':{'status':'VISIBLE_FEEDBACK_CONTRAST_IN_DECLARED_CLASS','not_claimed':['unique mediator','physical validation','diagnostic iid start','finite N rate','QM GR ToE']}},indent=2))

for t in range(1,6):
 label=f'MM-{t:03d}'
 if t==1: status='EXACT_FROZEN_CONTRACT_AND_REUSED_CONTROLS';summary='Frozen native and diagnostic labelled count-only five-time design, including essential observation-memory negative controls.'
 elif t==2:status='CONDITIONAL_STATIONARY_FAST_PATH_PROOF_SUBMITTED';summary='Derived diagnostic F/H/U dynamics and initial joint Gaussian limit; tightness conditional on inherited MT high moment/boundary payment.'
 elif t==3:status='EXACT_GAUSSIAN_FILTER_IDENTITY_WITH_INPUT_CONDITIONAL_NUMERICS';summary='Native filtered determinant zero; auxiliary strictly negative for two fixed kappas; gain/clock/reference restrictions explicit.'
 elif t==4:status='CONDITIONAL_FIXED_WINDOW_FINITE_PARENT_LAW_LIMIT';summary='Bounded ten-sign joint law and seven pair functionals converge, conditional on MM-002; no finite-data rates.'
 else:status='VISIBLE_FEEDBACK_CONTRAST_IN_DECLARED_CLASS';summary='Mechanism discrimination only within frozen statistical models; no physical field attribution.'
 js(f'{label}/RESULTS.json',{'task':label,'execution':'COMPLETED_AS_SCOPED','scientific_status':status,'summary':summary,'current_historical_sources_modified':False,'primary_FIN_changed':False,'new_physics_postulate':False,'next_action':'Independent HM/MT stationarity/tightness audit and one licensed physical interface only if supplied.'})
 save(f'{label}/REPORT.md',f'# {label} — report\n\n**Status:** `{status}`.\n\n{summary}\n\nSee adjacent specific proof/contract files, root `verification.json`, `INPUT_PROVENANCE.json`, and `NONCONCLUSIONS.md`. No earlier research package is reproduced here.')
 save(f'{label}/PROOF_OR_OPEN_LEMMA.md',f'# {label} — proof status\n\n**Disposition:** `{status}`.\n\n{summary}\n\nMM-002 and MM-004 remain conditional on separately admitted MT/KO local stationary moment, conditional CLT, reflected boundary and martingale-tightness obligations. Exact finite-state and Gaussian formulas are distinguished from accepted physical claims. See root `CLAIM_REGISTER.json` and task-specific files.')
 save(f'{label}/NONCONCLUSIONS.md',f'# {label} — restrictions\n\n{summary}\n\nNo unique hidden mechanism, physical apparatus, general colored-noise robustness, iid diagnostic transfer, finite-N performance guarantee, physical time/space, quantum mechanics, gravity or ToE is established. Earlier source results are dependencies and not counted as new {label} theorems.')
 save(f'{label}/REPLAY.md',f'# {label} — replay\n\nRun `PYTHONDONTWRITEBYTECODE=1 python3 verify_mm.py` from the root MM directory and inspect `verification.json`. It checks exact finite/generic algebra and declared nominal or interval calculations only. The terminal analytical argument is documented in task-specific files; it requires the historically accepted MT/KO/HM assumptions and must not be presented as mechanized proof.')
 save(f'{label}/INPUTS.sha256','\n'.join(lines))

# hash complete research file set, excluding verifier's own generated manifest and ZIP packaging products
excluded={'MANIFEST.sha256','PACKAGE_VERIFICATION.json','verify_package.py','build_mm_documents.py'}
items=[]
for p in sorted(ROOT.rglob('*')):
 if p.is_file() and p.name not in excluded and '__pycache__' not in p.parts:
  items.append(f'{sha256(p.read_bytes()).hexdigest()}  {p.relative_to(ROOT)}')
save('MANIFEST.sha256','\n'.join(items))
print('Written',len(items),'manifest entries')
