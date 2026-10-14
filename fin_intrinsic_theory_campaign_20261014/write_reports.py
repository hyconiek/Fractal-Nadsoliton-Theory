from pathlib import Path
import json, hashlib, textwrap, os, math
ROOT=Path('/mnt/data/fin_intrinsic_theory_campaign_20261014')

def w(rel, txt):
    p=ROOT/rel; p.parent.mkdir(parents=True,exist_ok=True); p.write_text(textwrap.dedent(txt).lstrip(),encoding='utf-8')

def j(rel,obj):
    p=ROOT/rel; p.parent.mkdir(parents=True,exist_ok=True); p.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+"\n",encoding='utf-8')

remote={
"repository_baseline":"1b7cd447 (as supplied by the user plan; remote byte provenance only for URLs explicitly inspected)",
"project_sources":[
 {"id":"T02","url":"https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/1b7cd447/fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py","use":"q12 strict-kernel construction and leave-one-out finite-N formulas","status":"REMOTE_INSPECTED_FULL_TEXT"},
 {"id":"T03a","url":"https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_rank7_mathphysics_next/proofs/MP7-031_finite_N_variational_model.md","use":"finite-copy Gibbs model and variational limit","status":"REMOTE_INSPECTED_FULL_TEXT"},
 {"id":"T03b","url":"https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_rank7_mathphysics_next/proofs/MP7-032_gaussian_auxiliary_representation.md","use":"auxiliary Gaussian identity; distinction from actual kinetics","status":"REMOTE_INSPECTED_FULL_TEXT"},
 {"id":"T03c","url":"https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/1b7cd447/fin_rank7_mathphysics_next/proofs/MP7-033_probability_mediator_fluctuations.md","use":"probability/mediator fluctuation dictionary","status":"REMOTE_INSPECTED_FULL_TEXT"},
 {"id":"T04a","url":"https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/1b7cd447/fin_rank7_mathphysics_next/proofs/MP7-025_simple_fold_normal_form.md","use":"validated fold interval and a,b intervals","status":"REMOTE_INSPECTED_FULL_TEXT"},
 {"id":"T04b","url":"https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/1b7cd447/fin_rank7_mathphysics_next/proofs/MP7-026_controlled_fold_scaling.md","use":"local fold remainder and transverse bounds","status":"REMOTE_INSPECTED_FULL_TEXT"},
 {"id":"T04c","url":"https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_rank7_mathphysics_next/proofs/MP7-029_conditional_slowing_and_stochastic_laws.md","use":"comparison with previously supplied auxiliary dynamics","status":"REMOTE_INSPECTED_FULL_TEXT"},
 {"id":"T05a","url":"https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/FIN_RESEARCH_AFTER_POST_AFTER_CONTINUATION_HANDOFF_131_255_20260927/01_REPORTS/MICROSCOPIC_PROCESS_CONTRACT_138.md","use":"leave-one-out is primary microscopic process and convention tags","status":"REMOTE_INSPECTED_FULL_TEXT"},
 {"id":"T05b","url":"https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/FIN_RESEARCH_AFTER_POST_AFTER_CONTINUATION_HANDOFF_131_255_20260927/01_REPORTS/RECURSIVE_LOCAL_CHILD_TEST_131.md","use":"finite-N drift/covariance expansion","status":"REMOTE_INSPECTED_FULL_TEXT"}
],
"literature":[
 {"id":"L01","url":"https://arxiv.org/html/1502.05185v3","title":"Kraaij (2015), Large deviations for Markov jump processes with mean-field interaction","use":"Theorem 2.3 applied condition-by-condition"},
 {"id":"L02","url":"https://arxiv.org/abs/1312.7591","title":"Mielke, Peletier, Renger (2014)","use":"context for LDP/gradient structures; no direct theorem transfer"},
 {"id":"L03","url":"https://arxiv.org/abs/1204.4503","title":"Cuff et al. (2012), mean-field Potts Glauber dynamics","use":"spinodal scaling benchmark only"},
 {"id":"L10","url":"https://projecteuclid.org/journals/annals-of-probability/volume-19/issue-4/Solutions-of-a-Stochastic-Differential-Equation-Forced-Onto-a-Manifold/10.1214/aop/1176990225.full","title":"Katzenberger (1991)","use":"methodological pointer only; not invoked as a theorem for the jump chain"}
],
"hash_note":"Web references were inspected through the browsing layer. Full remote bytes were not materialized in the container, so this package does not invent SHA-256 values for them. Local user-supplied inputs are hash-bound separately."
}
j('REMOTE_SOURCES.json',remote)

model_contract={
"execution_status":"COMPLETED",
"scientific_verdict":"CONSISTENT_INTRINSIC_MICROSCOPIC_PARENT_WITH_SCOPED_SOURCE_GAPS",
"q":12,
"kernel":{"omega":0.18575,"phi":0.1625,"eta":1.8,"construction":"circulant distance kernel -> Laplacian -> Fourier k3,k4,k5,k6 features -> A=XX^T","rank":7},
"finite_copy_equilibrium":{"labels":"0..11","N":"integer >=1","weight":"exp[(g/(2N)) sum_ab A[x_a,x_b]]","all_pairs":True,"physical_role":"dimensionless model premise"},
"primary_microscopic_kinetics":{"id":"LEAVE_ONE_OUT_GIBBS_HEAT_BATH","attempt_rate":"1 per labelled copy per model-time unit","count_jump":"n -> n-e_i+e_j","rate":"n_i softmax_j[(g/N)A(n-e_i)]","status":"SUPPLIED_MODEL_PREMISE"},
"limiting_objects":{"p":"n/N","s":"softmax(g A p)","mu":"X^T p","theta_ret":"g mu","drift_p":"s-p","Phi":"||theta||^2/(2g)-log(mean_j exp(X_j.theta))"},
"fold":{"g_interval":[3.5156447068395917,3.5156447268395917],"a_interval":[-0.21257164008401387,-0.212571626020459],"b_interval":[0.11898252745049878,0.11898266710322366],"certificate_sha_from_MP7_026":"5537174281996415a20ea3dffce23c33258c51396b3850aedf4168d9280e352e","status":"SUPPLIED_INTERVAL_ASSISTED_CERTIFICATE"},
"initial_law_classes_for_TH002":["any initial empirical-law sequence satisfying a good LDP","deterministic type sequences p_N->p0 (proved special case)","finite-copy Gibbs equilibrium initial law (method-of-types LDP special case)"],
"distinctions":{"theta_aux":"equilibrium Gaussian latent variable from MP7-032; not microscopic trajectory","empirical_refresh":"separate finite-N comparison convention; no O(1/N) transfer without generator comparison","model_time":"attempt clock; not physical seconds"},
"input_discrepancy":"The materialized Library AGENTS snapshot predates the user-supplied 2026-10-14 theory reorientation. Binding current campaign guardrails therefore come from the user-supplied TH plan and derived-bridges files; the older AGENTS snapshot is retained only as historical context."
}
j('TH-001/MODEL_AND_AXIOM_CONTRACT.json',model_contract)

j('TH-001/DEPENDENCY_DAG.json',{
"nodes":[
 {"id":"A7_X","type":"SUPPLIED_MODEL_OBJECT"},{"id":"finite_Gibbs","type":"SUPPLIED_MODEL_EXTENSION","depends_on":["A7_X"]},
 {"id":"LOO_generator","type":"SUPPLIED_KINETIC_PREMISE","depends_on":["A7_X","finite_Gibbs"]},
 {"id":"path_LDP","type":"DERIVED_THEOREM","depends_on":["LOO_generator","initial_LDP","Kraaij_2_3"]},
 {"id":"coarse_contraction","type":"DERIVED_THEOREM","depends_on":["path_LDP","observable_map"]},
 {"id":"fold_certificate","type":"SUPPLIED_VALIDATED_INTERVAL_RESULT","depends_on":["A7_X","Phi"]},
 {"id":"microscopic_soft_coefficients","type":"DERIVED_IDENTITY_PLUS_CERTIFICATE_REUSE","depends_on":["LOO_generator","fold_certificate"]},
 {"id":"stochastic_fold_limit","type":"OPEN_CONVERGENCE_TARGET","depends_on":["microscopic_soft_coefficients","fast_mode_averaging_lemma"]}
]})

w('TH-001/DEDUPLICATION.md',r'''
# TH-001 deduplication map

The following are **reused**, not new claims in this campaign:

- finite-copy Gibbs variational model and its thermodynamic free-energy limit (MP7-031);
- exact Hubbard–Stratonovich / auxiliary Gaussian representation (MP7-032);
- leading deterministic mean-field drift `p_dot=s(p)-p`, retained drift `theta_dot=-g grad Phi`, and leading jump covariance;
- certified static simple fold, square-root branch splitting, local `epsilon^(3/2)` energy difference and soft Hessian scaling (MP7-025/026);
- finite-N leave-one-out versus empirical-refresh convention split (report 138 / report 131);
- earlier C/D/E/F static elimination, memory and N2 operational results;
- incidence/locality no-go and all physical-source/units/selector gates.

New targets/results in TH-001–005 are instead: (i) a theorem-level path LDP application to the actual leave-one-out chain; (ii) its non-quadratic flux action and basin-local quasipotential consequence; (iii) an exact sufficiency/non-sufficiency map for path fluctuations; (iv) actual generator-derived fold noise and finite-N correction; and (v) a stopped stochastic-fold convergence statement or a precise convergence blocker.
''')

w('TH-001/THEOREM_APPLICATION_TABLE.md',r'''
# TH-001 theorem application table

## L01 — Kraaij (2015), Theorem 2.3: directly applied in TH-002

For the q-state Glauber-type empirical measure chain, Theorem 2.3 assumes a continuous limiting edge intensity `v(a,b,p)`, uniform convergence of `p_a r_N(a,b,p)` to that intensity, positivity whenever `p_a>0`, a monotone boundary factorization when `p_a=0`, and an initial-law LDP.

For FIN define

`v_ij(p)=p_i s_j(p)`, `s(p)=softmax(g A p)`.

Condition-by-condition:

1. **Continuity:** PASS. `s_j` is analytic on the compact simplex.
2. **Interior positivity:** PASS. `s_j(p)>0` for all finite logits, hence `v_ij>0` whenever `p_i>0`.
3. **Boundary factorization:** PASS with `v_dagger(x)=x` and `v_ddagger(p)=s_j(p)`. The first factor is increasing; the second is continuous and strictly positive at the boundary.
4. **Uniform rate convergence:** PASS for each fixed finite `g`, and uniformly on compact gain intervals. The finite rate is `p_i softmax_j[gA(p-e_i/N)]`. Coordinate softmax has l-infinity Lipschitz constant at most 1/2, so the discrepancy is at most `g max|A_ki|/(2N)`.
5. **Initial LDP:** CONDITIONAL in the theorem, but two concrete classes are paid here: deterministic type sequences converging to `p0`; and the declared Gibbs equilibrium sequence, whose type probabilities obey the method-of-types LDP with rate `F_g-min F_g`.

Therefore TH-002 may invoke Theorem 2.3 for these initial classes. This is an application to an interacting mean-field chain, not the independent-walker result ST8548.

## L02 — Mielke–Peletier–Renger (2014)

Used as structural context only. Their paper connects path large-deviation costs with generalized gradient structures. TH-002 does **not** transfer an independent-particle theorem to interacting FIN. Instead, the FIN detailed-balance/reversal identity is derived directly from its Hamiltonian and free energy.

## L03 — Cuff et al. (2012)

Benchmark only. The paper proves an `n^(-2/3)` spinodal window and `n^(4/3)` mixing scale in single-update steps for mean-field Potts. Since FIN uses total attempt rate N per model-time unit, `n^(4/3)` single updates corresponds dimensionally to `n^(1/3)` FIN model time. No Potts coefficient or mixing theorem is transferred.

## L10 — Katzenberger (1991)

Not invoked as a completed theorem application. It is a possible method for strongly stable directions, but its SDE/manifold hypotheses are not automatically a theorem for the present density-dependent jump chain. TH-005 isolates the exact fast-mode/perturbed-test-function lemma that would have to be paid before using such a method.

## Scope boundary L04

Uniform-in-time mean-field convergence results requiring one exponentially stable nonlinear stationary law are not used at the fold: the fold has a zero mode and nearby multiple stationary branches, so that premise is unavailable exactly where TH-005 operates.
''')

w('TH-001/REPORT.md',r'''
# TH-001 report — model, axioms and deduplication

**Execution status:** COMPLETED.

**Scientific verdict:** PASS_WITH_PROVENANCE_QUALIFICATION.

One microscopic parent is frozen: q=12 rank-seven FIN, declared finite-copy Gibbs law and exact leave-one-out heat-bath with one attempt per copy. `theta_ret=g X^T p` is an observable of that chain; `theta_aux` is only the equilibrium Gaussian latent variable. The supplied fold certificate is numerically replayed against the same strict-kernel construction in TH-004.

The current user-supplied TH plan is newer than the materialized Library AGENTS snapshot, so the package records that discrepancy rather than silently claiming byte-identical repository state.
''')

j('TH-001/RESULTS.json',{"execution_status":"COMPLETED","scientific_verdict":"PASS_WITH_PROVENANCE_QUALIFICATION","supplied_objects":["strict kernel parameters","A=XX^T","finite-copy Gibbs extension","leave-one-out heat bath","attempt clock","fold interval certificate"],"derived_objects":["dependency DAG","theorem application map","deduplication map"],"fitted_objects":[],"arithmetic_status":"hash-bound local inputs; remote theorem/project references inspected but not locally byte-hashed","next_action":"TH-002"})
w('TH-001/PROOF_OR_OPEN_LEMMA.md','''# TH-001 proof/open-lemma status\n\nNo mathematical blocker was found in binding the finite parent. Remote-byte provenance for some repository proof files is not locally materialized; numerical/theorem claims therefore cite their inspected source URLs and supplied certificate identifiers rather than inventing hashes.\n''')
w('TH-001/NONCONCLUSIONS.md','''# TH-001 nonconclusions\n\nThis task does not source A7, the heat-bath law, N, gain, physical time, units, locality, a physical auxiliary field, QM/GR/SM, or empirical validation.\n''')
w('TH-001/REPLAY.md','''# TH-001 replay\n\n1. Inspect `source_context/FIN_INTRINSIC_THEORY_NEXT_CAMPAIGN_20261014_USER.md`.\n2. Inspect `REMOTE_SOURCES.json`.\n3. Run `python theory_core.py` and `python independent_verify.py`.\n4. Compare the fold replay and model spectrum with TH-004 certificates.\n''')

# TH002
w('TH-002/GENERATOR_TO_ACTION_PROOF.md',r'''
# TH-002 — interacting path action from the exact leave-one-out count generator

## 1. Exact finite-N nonlinear generator

Let `p=n/N` and `z_ij=e_j-e_i`. The exact count generator is

`L_N f(p)=N sum_ij p_i s_j^{(-i,N)}(p) [f(p+z_ij/N)-f(p)]`,

where

`s_j^{(-i,N)}(p)=softmax_j[g A(p-e_i/N)]`.

For a linear exponential test `f(p)=xi.p`,

`H_N(p,xi)=(1/N)e^{-Nf} L_N e^{Nf}`

is exactly

`sum_ij p_i s_j^{(-i,N)}(p) [exp(xi_j-xi_i)-1]`.

For bounded `g` and the finite matrix A, `s^{(-i,N)} -> s(p)=softmax(gAp)` uniformly. Thus

`H(p,xi)=sum_ij p_i s_j(p)[exp(xi_j-xi_i)-1]`.

It factorizes as

`H=(sum_i p_i e^{-xi_i})(sum_j s_j e^{xi_j})-1`,

and is invariant under `xi -> xi+c 1`.

## 2. Kraaij Theorem 2.3 applies

Set `v_ij(p)=p_i s_j(p)`. The coordinate softmax derivative satisfies

`sum_k |partial s_j/partial z_k|=2 s_j(1-s_j)<=1/2`.

Hence

`sup_p |p_i s_j[gA(p-e_i/N)]-p_i s_j(gAp)| <= g max_k |A_ki|/(2N)`.

The limiting intensity is continuous; it is strictly positive whenever `p_i>0`; and at a boundary face it factors as `p_i * s_j(p)` with the first factor increasing and the second positive continuous. These are exactly the boundary features required by Kraaij's finite-state Glauber theorem. Consequently, for every initial-law sequence satisfying a good LDP, the empirical path satisfies a path LDP with speed N and Lagrangian rate.

Two initial classes are explicitly paid here:

- deterministic type states `p_N -> p0`, giving `I0(p)=0` at p0 and infinity otherwise;
- the declared Gibbs equilibrium sequence. MP7-031 method-of-types bounds give its LDP with rate `F_g(p)-min F_g`, where `F_g=D(p||u)-(g/2)p^T A p` up to the fixed uniform convention.

## 3. Flux form of the exact limiting Lagrangian

Let `c_ij(p)=p_i s_j(p)` and `ell(q;c)=q log(q/c)-q+c`, with the standard conventions `ell(0;0)=0` and `ell(q;0)=infinity` for q>0. Convex duality for the edge-Poisson Hamiltonian gives

`L(p,v)= inf_{q_ij>=0: sum_ij q_ij z_ij=v} sum_ij ell(q_ij;c_ij(p))`.

This is the non-quadratic path action forced by the generator. At interior p all directed edge rates are positive, so every simplex-tangent velocity (`sum_i v_i=0`) has finite cost. At a boundary face, finite cost forbids outgoing flux from an empty state; equivalently the velocity must lie in the simplex tangent cone.

The zero-cost flux is `q_ij=c_ij`; its velocity is `s-p`, reproducing the mean-field ODE without assuming a quadratic action.

## 4. Exact reversal identity and a basin-local quasipotential

For interior p define, modulo an irrelevant constant,

`(grad F)_i=log p_i-g(Ap)_i`.

Then

`p_i s_j exp[(grad F)_j-(grad F)_i]=p_j s_i`.

Relabelling i,j yields the exact symmetry

`H(p,xi)=H(p,grad F-xi)`.

Legendre duality therefore gives

`L(p,v)=grad F(p).v + L(p,-v)`.

Since `L>=0`, every path from p_a to p_b satisfies

`Action >= F(p_b)-F(p_a)`.

If p_b lies in the deterministic basin of an interior stable minimum m, the time-reversed zero-cost relaxation trajectory from p_b toward m approaches this lower bound. Thus, for targets connected to m by the basin relaxation orbit, the basin-local quasipotential is

`V_m(p)=F(p)-F(m)`.

This is local/basin-scoped. It does not identify the global communication height, `Gamma=B4`, a capacity prefactor, or a physical least-action principle.

## 5. Gaussian truncation is not the path law

At the fold equilibrium, the actual soft theta jump exponent has second cumulant rate `2 g_f` but a nonzero fourth cumulant rate about `120.5038353`; the cubic rate vanishes there by the equilibrium i<->j symmetry. The quartic coefficient in the local Hamiltonian expansion is about `5.02099314` before any critical rescaling. Hence equal drift/covariance does not make the full jump Hamiltonian Gaussian. The diffusion in TH-005, if proved, is a critical scaling limit in which higher jump cumulants vanish—not an exact replacement of the path action.
''')

j('TH-002/HAMILTONIAN.json',{
"finite_N":"H_N(p,xi)=sum_ij p_i softmax_j[g A(p-e_i/N)](exp(xi_j-xi_i)-1)",
"limit":"H(p,xi)=sum_ij p_i s_j(p)(exp(xi_j-xi_i)-1)",
"factorization":"(sum_i p_i exp(-xi_i))(sum_j s_j exp(xi_j))-1",
"gauge":"H(p,xi+c1)=H(p,xi)",
"drift":"d_xi H|0=s-p",
"quadratic":"d_xi^2 H|0 = diag(p)+diag(s)-p s^T-s p^T",
"reversal":"H(p,xi)=H(p,grad F-xi)",
"lagrangian_flux":"inf over q_ij>=0 with sum q_ij(e_j-e_i)=v of sum [q log(q/c)-q+c], c_ij=p_i s_j",
"fold_soft_cumulants_nominal":{"second":7.031289433679185,"third":-6.108395039783332e-16,"fourth":120.50383529864993}
})

w('TH-002/LDP_ASSUMPTIONS.md',r'''
# TH-002 LDP assumptions and statuses

| Kraaij Theorem 2.3 premise | FIN status |
|---|---|
| finite q-state mean-field empirical chain | PASS, q=12 |
| continuous limiting `v_ij(p)` | PASS, analytic softmax |
| positive when source occupation positive | PASS |
| boundary monotone factorization | PASS: `v=p_i*s_j` |
| uniform finite-N rate convergence | PASS, O(1/N) for fixed/compact g |
| initial-law good LDP | CONDITIONAL generally; deterministic and equilibrium Gibbs classes paid |
| topology | Skorokhod path space on compact simplex, as in theorem |

The theorem proves a path LDP for the declared mathematical chain. It says nothing about physical time, a natural source of A7, or experimental data.
''')

j('TH-002/CUMULANT_DIAGNOSTICS.json',{"location":"fold equilibrium, soft retained theta direction","status":"NUMERICAL_DIAGNOSTIC_BOUND_TO_BINARY64_REPLAY","second_cumulant_rate":7.031289433679185,"expected_exact_identity":"2*g_f","third_cumulant_rate":-6.108395039783332e-16,"fourth_cumulant_rate":120.50383529864993,"quartic_H_coefficient":5.020993137443747,"interpretation":"non-quadratic jump action; not a non-Gaussian white-noise claim"})

w('TH-002/REPORT.md','''# TH-002 report\n\n**Execution status:** COMPLETED.\n\n**Scientific verdict:** PATH_LDP_PROVED_FOR_DECLARED_INITIAL_LDP_CLASSES.\n\nKraaij Theorem 2.3 applies to the exact leave-one-out chain because the limiting rate is `p_i s_j(p)`, the boundary factor is exactly the occupation `p_i`, and the self-exclusion changes the categorical softmax by O(1/N) uniformly. The resulting path cost is the non-quadratic edge-flux relative-entropy action. A direct detailed-balance symmetry yields a basin-local free-energy quasipotential statement without borrowing the unpaid global barrier theorem.\n''')
j('TH-002/RESULTS.json',{"execution_status":"COMPLETED","scientific_verdict":"PATH_LDP_PROVED_FOR_DECLARED_INITIAL_LDP_CLASSES","proof_scope":"Kraaij Theorem 2.3 + direct rate verification; deterministic and Gibbs initial LDP classes explicitly paid","supplied_objects":["LOO generator","A","g","initial law class"],"derived_objects":["limiting Hamiltonian","flux Lagrangian","reversal identity","basin-local quasipotential"],"fitted_objects":[],"arithmetic_status":"analytic theorem application plus nominal cumulant replay","next_action":"TH-003"})
w('TH-002/PROOF_OR_OPEN_LEMMA.md','''# TH-002 proof/open lemma\n\nThe path-LDP theorem is closed for any initial sequence satisfying the required good initial LDP; two concrete initial classes are proved. No global communication-height/capacity statement is claimed.\n''')
w('TH-002/NONCONCLUSIONS.md','''# TH-002 nonconclusions\n\n`H` is a stochastic large-deviation Hamiltonian, not conserved energy or a quantum Hamiltonian. The action is not a physical principle of least action. No global `Gamma=B4`, Eyring–Kramers prefactor, physical inertia, or empirical law is established.\n''')
w('TH-002/REPLAY.md','''# TH-002 replay\n\nRead `GENERATOR_TO_ACTION_PROOF.md` together with Kraaij Theorem 2.3 at https://arxiv.org/html/1502.05185v3. Run `python independent_verify.py` for the algebraic/numerical identities.\n''')

# TH003
w('TH-003/OBSERVABLE_SUFFICIENCY.md',r'''
# TH-003 — observable sufficiency and the minimal information lift

## 1. Mean closure is not path closure

The deterministic retained mean `mu=X^T p` closes because `A p=X mu`, so

`dot mu=X^T softmax(g X mu)-mu`.

The tilted retained Hamiltonian is instead

`H_mu(p,eta)=H(p,X eta)`

and contains the hidden factor `sum_i p_i exp[-X_i.eta]`. It is generally not a function of mu alone.

A finite exact witness uses the g=0 control and N=120. Let

`h_j=cos(2*pi*2*j/12)` (the hidden k=2 cosine),

`n_+=10+2h`, `n_-=10-2h`.

The vector `2h=(2,1,-1,-2,-1,1,2,1,-1,-2,-1,1)` is integer-valued, so these are genuine count states. Both have `mu=0` because k2 is orthogonal to the retained k3,k4,k5,k6 feature space.

For the retained c5 feature, `x_c5(j)^2=(Lambda5/12)[1+cos(4*pi*j/12)]`. Direct generator evaluation at g=0 therefore gives

`L(mu_c5^2)(p_+)-L(mu_c5^2)(p_-)=Lambda5/(60*120)=Lambda5/7200`.

With the supplied kernel `Lambda5≈2.298606272079097`, this is `0.0003192508711220968`. Thus equal retained means do not imply equal fluctuation generators.

## 2. An exact 11-dimensional linear lift of the occupation simplex

Write raw Fourier moments

`C_k=sum_j p_j cos(2*pi*k*j/12)`, `S_k=sum_j p_j sin(2*pi*k*j/12)`.

The seven retained feature means determine `C_3,S_3,C_4,S_4,C_5,S_5,C_6` after dividing by their known feature normalizations.

Four entries of `S(p)=X^T diag(p) X` determine four raw products:

- `P34c=E[c3 c4]=(C1+C5)/2`,
- `P34s=E[s3 c4]=(-S1-S5)/2`,
- `P35c=E[c3 c5]=(C2+C4)/2`,
- `P35s=E[s3 c5]=(-S2-S4)/2`.

Hence

`C1=2P34c-C5`,
`S1=-2P34s-S5`,
`C2=2P35c-C4`,
`S2=-2P35s-S4`.

All eleven nonconstant real Fourier coordinates are then known. Fourier inversion recovers the full occupation vector:

`p_j=(1/12)[1+2 sum_{k=1}^5 (C_k cos(k theta_j)+S_k sin(k theta_j))+C_6(-1)^j]`.

The numerical linear map consisting of normalization, seven retained means and these four products has rank 12 with smallest singular value about 0.8601.

### Minimality convention

On the continuous interior of the 11-dimensional simplex, no continuous injective coordinate chart can use fewer than 11 real coordinates. Therefore the **seven retained means plus four selected second moments** are dimension-minimal for exact continuous occupation recovery in this linear/Fourier convention. This is not a finite-type coding theorem: for a fixed finite N, an arbitrary real scalar can encode a finite set if one allows nonphysical discontinuous lookup encodings.

The lift recovers p, not the identities of labelled copies. It is therefore an exact predictive-state lift of the exchangeable occupation process, not eleven new material degrees of freedom.
''')

j('TH-003/MOMENT_LIFT_MAP.json',{
"retained_raw_modes":["C3","S3","C4","S4","C5","S5","C6"],
"additional_second_moments":["E[c3*c4]","E[s3*c4]","E[c3*c5]","E[s3*c5]"],
"recovery":{"C1":"2 E[c3*c4]-C5","S1":"-2 E[s3*c4]-S5","C2":"2 E[c3*c5]-C4","S2":"-2 E[s3*c5]-S4"},
"full_inversion":"12-point real Fourier inversion",
"matrix_rank_with_normalization":12,
"continuous_simplex_nonconstant_dimension":11,
"minimality_scope":"continuous exact occupation recovery / continuous linear observable convention; not arbitrary coding of a finite type set"
})

w('TH-003/PATH_CONTRACTION.md',r'''
# TH-003 path contraction

By TH-002, p_N obeys a path LDP with rate `I[p(.)]`. For the continuous observation `mu=X^T p`, the contraction principle gives

`I_mu[m(.)]=inf { I[p(.)] : X^T p(t)=m(t) for all t }`.

The infimum is over **entire hidden paths and the compatible initial-law cost**. An instantaneous minimization of the hidden state at each time is not generally equivalent and is not used.

Because `H(p,X eta)` varies inside a fibre of fixed mu (the N=120 witness already separates second derivatives), there is no universal exact Markov Hamiltonian on mu alone for arbitrary preparation. Three scopes must remain separate:

1. finite-N exact nonclosure;
2. local/stationary Gaussian reductions where hidden fluctuations can be slaved or averaged;
3. higher-order/nonstationary history dependence.

If exact path prediction is required without a reduction theorem, use full p or the exact 11-coordinate lift from `OBSERVABLE_SUFFICIENCY.md`.
''')

j('TH-003/SCOPE_AND_ERROR_REGISTER.json',{
"mean_only":{"deterministic_drift":"closed","path_fluctuations":"not closed in general","witness":"N=120 g=0 exact count-state pair"},
"lift":{"dimension":11,"status":"exact occupation reconstruction","compression":"none beyond simplex dimension"},
"contraction":{"status":"theorem via path-LDP contraction","hidden_object":"full hidden path and initial law"},
"memory":"may appear after projection; no fitted poles introduced"})

w('TH-003/REPORT.md','''# TH-003 report\n\n**Execution status:** COMPLETED.\n\n**Scientific verdict:** EXACT_MEAN_NONCLOSURE_PLUS_EXACT_FULL_OCCUPATION_LIFT.\n\nThe retained seven-dimensional mean is sufficient for the leading deterministic drift but not for the tilted path law. A genuine N=120 count-state witness separates fluctuation generators. Four existing second-moment observables complete the retained seven modes to an exact 11-dimensional occupation coordinate system.\n''')
j('TH-003/RESULTS.json',{"execution_status":"COMPLETED","scientific_verdict":"EXACT_MEAN_NONCLOSURE_PLUS_EXACT_FULL_OCCUPATION_LIFT","derived_objects":["N=120 hidden-fibre witness","path contraction formula","11D exact lift"],"supplied_objects":["X","LOO generator","existing harmonic identities"],"fitted_objects":[],"arithmetic_status":"analytic trigonometric identities; binary64 independent replay","next_action":"TH-004"})
w('TH-003/PROOF_OR_OPEN_LEMMA.md','''# TH-003 proof/open lemma\n\nNo open lemma is needed for the finite witness, Fourier lift, or contraction statement. A genuinely lower-dimensional stochastic sufficient state would require an additional reduction theorem; none is asserted.\n''')
w('TH-003/NONCONCLUSIONS.md','''# TH-003 nonconclusions\n\nThe second moments are not spacetime geometry. Exact occupation recovery does not recover labelled-copy identity, spatial incidence, particles, or quantum state. Mean nonclosure is not a proof that every asymptotic low-dimensional reduction fails.\n''')
w('TH-003/REPLAY.md','''# TH-003 replay\n\nRun `python theory_core.py` and `python independent_verify.py`. They rebuild the N=120 witness and check that the 12-row normalization/mean/product transform has rank 12.\n''')

# TH004
coeff={
"source_certified":{"g_fold":[3.5156447068395917,3.5156447268395917],"a":[-0.21257164008401387,-0.212571626020459],"b":[0.11898252745049878,0.11898266710322366],"transverse_lower_bound_MP7_026":0.00930411916943287,"H4_second_eigen_lower_MP7_026":0.0172259553883460,"tensor_T3_upper":0.222833872594371,"tensor_T4_upper":0.190171798000148},
"derived_interval_arithmetic":{"soft_noise_variance_2g":[7.0312894136791835,7.031289453679183],"soft_noise_sigma_sqrt_2g":[2.6516578613537574,2.651657868896209],"canonical_Z_scale":[3.2274443594317104,3.2274456221415457],"canonical_lambda_factor":[0.34303154768747873,0.34303170459029947],"canonical_time_factor":[1.4814348295357251,1.4814359971608897]},
"binary64_replay":{"theta_fold_even":[1.3643114281782407,1.4331027014285425,1.4080457012550374,1.0056820773201627],"g":3.5156447168395917,"v":[0.5073675396772498,0.5286856553138111,0.5537606940687492,0.39549810524420653],"a":-0.2125716330522363,"b":0.11898259727685968,"even_hessian_eigs":[1.700355114141964e-17,0.17800274933964866,0.21956286982117357,0.23310794393635673],"v_dot_theta":2.6273338943707807}
}
j('TH-004/COEFFICIENT_CERTIFICATE.json',coeff)

w('TH-004/SOFT_MODE_MICROSCOPIC_BRIDGE.md',r'''
# TH-004 — microscopic soft-mode bridge

## 1. Actual retained process

For the actual occupation chain define `theta_ret=g X^T p_N`. Its leading drift is

`dot theta_ret=-g grad Phi_g(theta_ret)`.

The leading covariance rate is

`(1/N) g^2 X^T D(p) X`,

where `D=diag(p)+diag(s)-p s^T-s p^T`. This is derived from the jump generator; no mobility/noise matrix is selected at the retained level.

At an equilibrium p=s, `D=2 Sigma`, `Sigma=diag(p)-pp^T`, so

`D_theta=2 g^2 M`, `M=X^T Sigma X`.

At the certified fold, the null-vector equation for the Hessian `H=I/g-M` gives `M_f v=v/g_f`. Therefore the **actual microscopic soft variance coefficient** is exactly

`v^T D_theta v = 2 g_f`.

Using the supplied fold interval gives

`2 g_f in [7.0312894136791835, 7.031289453679183]`,

and noise amplitude

`sqrt(2 g_f) in [2.6516578613537574, 2.651657868896209]`.

This proves the candidate TH-005 noise amplitude without importing MP7-029's auxiliary Langevin mobility.

## 2. Drift coefficients

The supplied validated fold coefficients are

`a=v.partial_g grad Phi in [-0.2125716400840139,-0.2125716260204590]`,

`b=D^3 Phi[v,v,v] in [0.1189825274504988,0.1189826671032237]`.

A direct binary64 solve of the same strict-kernel model gives

`theta_f≈(1.36431142818,1.43310270143,1.40804570126,1.00568207732)`,

`v≈(0.507367539677,0.528685655314,0.553760694069,0.395498105244)`,

`a≈-0.212571633052`, `b≈0.118982597277`, inside the supplied intervals.

## 3. Leave-one-out finite-N correction

Report 131 gives

`b_N(p)=s-p-(g/N) S_s A p+O(N^-2)`.

Multiplying by `g X^T` and projecting on v at the fold yields

`v.theta_dot_finite = -(v.theta_f)/N + O(N^-2)`

because `M_f v=v/g_f`. After the TH-005 amplitude/time rescaling this contributes

`-(v.theta_f) N^(-1/3)+o(N^(-1/3))`,

so it vanishes but is the leading convention correction. The nominal coefficient is `v.theta_f≈2.62733389437`. The finite-N diagnostic in TH-005 visibly converges to this coefficient.

## 4. Auxiliary versus actual

`theta_aux | p ~ Normal(g X^T p,(g/N)I)` is an exact equilibrium identity from MP7-032. It is not the path law of `theta_ret`. Only the soft fold channel happens to have the leading actual variance `2g_f/N`; other retained directions generally have `2g^2 M/N`, not isotropic `2g I/N`.
''')

w('TH-004/AUXILIARY_VS_ACTUAL_DICTIONARY.md',r'''
# TH-004 auxiliary-vs-actual dictionary

| Object | Meaning | Path status |
|---|---|---|
| `theta_ret=g X^T p_N` | observable of the leave-one-out jump chain | actual microscopic collective observable |
| `theta_aux` | Hubbard–Stratonovich latent Gaussian | equilibrium representation only |
| `-g grad Phi` | leading drift of theta_ret | derived from supplied heat bath |
| `g^2 X^T D(p)X/N` | leading covariance rate of theta_ret | derived from jumps |
| `2L/N` in MP7-029 | covariance of a separately supplied Langevin model | not transferable to actual labels |
| `sqrt(2g_f)` soft amplitude | fold projection of actual jump covariance | derived identity, source-certificate conditioned |
''')

w('TH-004/LOCAL_REMAINDER_SCOPE.md',r'''
# TH-004 local remainder scope

MP7-026 already pays the static fold Taylor/reduction remainder on its certified tube: Cartesian radius 0.002, `0<=epsilon<=1e-7`, transverse inverse bound, `T3/T4` tensor bounds, and `|r-r0|<=0.005965415221 epsilon` on the final branch bracket.

For TH-005, the stochastic path is not restricted to the stationary branch bracket. Smoothness of the finite log-sum-exp gives local Taylor remainders of schematic order

`O(epsilon*|xi| + |xi|^3 + |xi|*||y|| + ||y||^2 + epsilon||y||)`

in the soft drift. Under `epsilon=N^-2/3`, `xi=N^-1/3 Z`, these terms vanish after multiplication by `N^(2/3)` **provided the stable coordinate is controlled in the time-integrated/perturbed-test sense**. Paying that stochastic fast-mode estimate is the remaining TH-005 convergence lemma; the static MP7 remainder alone does not prove it.
''')

w('TH-004/REPORT.md','''# TH-004 report\n\n**Execution status:** COMPLETED.\n\n**Scientific verdict:** MICROSCOPIC_SOFT_DRIFT_AND_NOISE_COEFFICIENTS_DERIVED.\n\nThe actual leave-one-out generator, not the auxiliary Gaussian field, fixes the soft-channel noise. At the certified fold the exact algebraic identity is `variance coefficient = 2 g_f`; the supplied interval therefore yields a certified coefficient range. The a,b normal-form intervals are reused from MP7-025/026 and independently reproduced at binary64 precision.\n''')
j('TH-004/RESULTS.json',{"execution_status":"COMPLETED","scientific_verdict":"MICROSCOPIC_SOFT_DRIFT_AND_NOISE_COEFFICIENTS_DERIVED","supplied_objects":["fold interval root","a,b interval coefficients","LOO generator"],"derived_objects":["soft variance=2g_f","finite-N LOO correction in soft drift","auxiliary-vs-actual dictionary"],"fitted_objects":[],"arithmetic_status":"a,b,g reuse supplied interval-assisted certificate; derived scalar interval propagation + binary64 replay","next_action":"TH-005"})
w('TH-004/PROOF_OR_OPEN_LEMMA.md','''# TH-004 proof/open lemma\n\nThe coefficient identities are closed. The remaining stochastic theorem needs a fast stable-mode/hidden-mode averaging estimate on the N^(1/3) time scale; this is not supplied by the static fold certificate.\n''')
w('TH-004/NONCONCLUSIONS.md','''# TH-004 nonconclusions\n\nThe auxiliary Gaussian field is not fundamental, the relaxation rate is not mass, the attempt clock is not physical seconds, and isotropic retained noise is not established away from the soft fold projection.\n''')
w('TH-004/REPLAY.md','''# TH-004 replay\n\nRun `python theory_core.py`. Check `COEFFICIENT_CERTIFICATE.json` against MP7-025/026 source intervals listed in `REMOTE_SOURCES.json`. Run `python independent_verify.py` for an independent numerical construction.\n''')

# TH005
w('TH-005/STOCHASTIC_FOLD_LIMIT.md',r'''
# TH-005 — stopped stochastic-fold limit: derived candidate and exact blocker

## 1. Frozen double scaling

Let the supplied fold be `(theta_f,g_f,v)`. Fix a finite real lambda and set

`g_N=g_f+lambda N^(-2/3)`.

Define the soft observable

`Z_N(s)=N^(1/3) v.[theta_ret(N^(1/3)s)-theta_f]`,

stopped on exit from a fixed local fold neighbourhood. The microscopic process remains the same leave-one-out chain; only N and g are taken along an already-defined mathematical family.

## 2. Rescaled generator characteristics

A microscopic label jump i->j changes theta by

`Delta theta=(g_N/N)(X_j-X_i)`,

so `Delta Z=O(N^(-2/3))`. On s-time, the total event rate is O(N^(4/3)). Therefore:

- third and higher jump characteristics vanish (`third absolute characteristic = O(N^-2/3)` for bounded local states);
- the quadratic characteristic has an O(1) limit;
- the soft drift is multiplied by `N^(2/3)`.

Using TH-004 and the simple-fold Taylor expansion,

`N^(2/3) v.theta_dot -> -g_f[a lambda+(b/2)Z^2]`,

provided stable directions are negligible in the required averaged/perturbed-test sense. The quadratic characteristic tends to `2g_f` once the full occupation is locally driven to p_f.

The resulting candidate martingale problem is

`dZ=-g_f[a lambda+(b/2)Z^2] ds + sqrt(2g_f)dB_s`,

up to the declared stopping time.

## 3. What is already proved

The following ingredients are not conjectural:

1. simple fold and source intervals for g_f,a,b;
2. full path LDP of the microscopic chain (TH-002);
3. exact deterministic retained drift;
4. exact soft covariance identity `2g_f` at the fold;
5. leave-one-out O(1/N) drift correction, which becomes O(N^-1/3) under this scaling;
6. bounded jump sizes and vanishing third characteristic under the scaling;
7. a strictly stable six-dimensional transverse Hessian at the fold (source interval certificate plus numerical full-7D replay).

The Potts benchmark independently has an n^-2/3 spinodal window and n^(4/3) single-update critical scale. Converting its one-update clock by N updates per FIN model-time unit gives N^(1/3), matching the balance above. This is a universality diagnostic, not a transferred theorem.

## 4. Exact remaining lemma

A complete convergence theorem still requires a **FAST_MODE_AVERAGING_AND_CONTAINMENT lemma** for the actual jump chain:

For prepared local initial sequences and up to the stopped fold exit, prove that the six retained stable directions and the four occupation-fibre directions contribute `o_p(1)` to the rescaled soft martingale problem. Concretely, one needs a uniform time-integrated bound strong enough to make

`N^(2/3) integral [ |xi| ||y|| + ||y||^2 + hidden-dependent finite-N terms ] dt`

vanish after the s-time conversion, together with convergence of the soft quadratic characteristic to `2g_f`.

The mean-only path nonclosure from TH-003 prevents simply deleting the hidden occupation variables. Their deterministic fibre differences contract, and ordinary stable fluctuations are expected at O(N^-1/2), but that expectation is not substituted for the required stochastic estimate.

Katzenberger-style strong-drift methods may be useful, but the cited theorem is for SDEs and is not silently applied to this jump chain. A direct perturbed-test-function/martingale proof is the cleanest missing step.

## 5. Bounded finite-N diagnostic

For rounded fold occupations at N=120,240,480,960,1920,3840 the scaled soft variance approaches `2g_f`, the scaled absolute third jump characteristic decreases, and the drift residual behaves approximately as the predicted leave-one-out correction `(v.theta_f) N^-1/3`.

At lambda=0, `-drift_error*N^(1/3)` runs approximately

`2.663, 2.646, 2.637, 2.632, 2.630, 2.628`,

approaching the independently computed `v.theta_f≈2.627333894`. This is a diagnostic of the expansion, not a convergence proof.

## 6. Canonical dimensionless form

Let

`c=(4/b)^(1/3)`, `d=c^2/(2g_f)`, `Y=Z/c`, `tau=s/d`.

Then the candidate becomes

`dY=(Lambda-Y^2)d tau+dW`,

with

`Lambda=(-a c/2) lambda`.

Using the supplied a,b,g intervals,

- `c in [3.2274443594,3.2274456221]`,
- `d in [1.4814348295,1.4814359972]`,
- `(-a c/2) in [0.3430315477,0.3430317046]`.

Thus stopped first-exit and response laws, once the convergence lemma is paid, depend only on the dimensionless Lambda and normalized initial/boundary conditions. These are model-time combinations, not physical universal constants.
''')

j('TH-005/LIMIT_CONTRACT.json',{
"g_N":"g_f+lambda*N^(-2/3)","soft_coordinate":"Z_N=N^(1/3) v.(theta_ret-theta_f)","model_time":"t=N^(1/3)s","stop":"exit from fixed local interior fold neighbourhood","initial_scope":"prepared local sequences; soft coordinate converges and nonsoft occupation deviations are sufficiently controlled","candidate_SDE":"dZ=-g_f(a lambda+b Z^2/2)ds+sqrt(2g_f)dB","proved_components":["fold coefficients","microscopic drift","soft covariance identity","jump-size scaling","higher jump characteristic vanishing","LOO O(1/N) correction scaling"],"open_component":"FAST_MODE_AVERAGING_AND_CONTAINMENT","status":"FORMAL_GENERATOR_LIMIT_WITH_PRECISE_OPEN_LEMMA"})

w('TH-005/CONVERGENCE_OR_BLOCKER.md',r'''
# TH-005 convergence disposition

**Verdict: OPEN_ONE_LEMMA, NOT A FAILED SCALING.**

The rescaled drift, covariance, jump size and finite-N convention correction are all consistent with the proposed stopped diffusion. What is not yet rigorously paid is the singular-perturbation step that removes six stable retained coordinates plus the four occupation-fibre coordinates on an N^(1/3) time horizon.

A sufficient next theorem would establish, under a frozen prepared-initial-law class and stopped local neighbourhood, a perturbed-test-function estimate for every smooth compactly supported phi(Z):

`G_N(phi o Z + correction_N) -> G phi`

uniformly in probability/integrated occupation, where

`G phi = -g_f(a lambda+b z^2/2) phi' + g_f phi''`.

Together with tightness and uniqueness of the stopped one-dimensional martingale problem, this closes TH-005. Numerical plots, the static transverse Hessian alone, or assuming p=s(theta) at every time do not close it.
''')

j('TH-005/UNIVERSAL_COMBINATIONS.json',{
"canonical_SDE":"dY=(Lambda-Y^2)d_tau+dW",
"Z_scale_interval":[3.2274443594317104,3.2274456221415457],
"s_to_tau_interval":[1.4814348295357251,1.4814359971608897],
"Lambda_per_lambda_interval":[0.34303154768747873,0.34303170459029947],
"conditional_predictions":["normalized first-exit distribution","normalized survival curve","normalized response to lambda"],
"status":"CONDITIONAL_ON_FAST_MODE_CONVERGENCE_AND_FIXED_STOPPING_CONTRACT"})

w('TH-005/REPORT.md','''# TH-005 report\n\n**Execution status:** COMPLETED_WITH_OPEN_LEMMA.\n\n**Scientific verdict:** STOCHASTIC_FOLD_GENERATOR_BALANCE_DERIVED__FULL_CONVERGENCE_NOT_YET_PROVED.\n\nThe proposed N^-2/3 / N^-1/3 / N^1/3 scaling is internally forced by the actual microscopic drift and jump variance. The candidate SDE coefficients are derived, not fitted, and the finite-N leave-one-out correction vanishes at the predicted N^-1/3 rate. A complete stopped convergence theorem is blocked by one explicitly isolated fast-mode/hidden-mode averaging lemma.\n''')
j('TH-005/RESULTS.json',{"execution_status":"COMPLETED_WITH_OPEN_LEMMA","scientific_verdict":"STOCHASTIC_FOLD_GENERATOR_BALANCE_DERIVED__FULL_CONVERGENCE_NOT_YET_PROVED","supplied_objects":["LOO generator","fold certificate","attempt clock"],"derived_objects":["critical scaling balance","candidate stopped diffusion","finite-N correction order","canonical dimensionless form"],"fitted_objects":[],"proof_scope":"generator characteristics and diagnostics proved/checked; process convergence conditional on one fast-mode averaging lemma","arithmetic_status":"source interval coefficients + binary64 finite-N diagnostics","next_action":"architect review; if continued, prove FAST_MODE_AVERAGING_AND_CONTAINMENT only"})
w('TH-005/PROOF_OR_OPEN_LEMMA.md','''# TH-005 proof/open lemma\n\nOpen lemma: `FAST_MODE_AVERAGING_AND_CONTAINMENT`. It must control the stable retained coordinates and hidden occupation fibre in the stopped N^(1/3)-time martingale problem. No other new microscopic law is needed for the candidate limit.\n''')
w('TH-005/NONCONCLUSIONS.md','''# TH-005 nonconclusions\n\nNo global mixing theorem, global cubic equilibrium, physical time, propagating particle, spatial QFT, non-Gaussian white noise, physical universality, natural A7 selection, QM/GR/SM, or ToE closure is claimed. Brownian increments in the candidate limit are Gaussian; non-Gaussianity refers to the nonlinear state/exit law.\n''')
w('TH-005/REPLAY.md','''# TH-005 replay\n\nRun `python fold_diagnostics.py` to regenerate the finite-N characteristic table. Run `python independent_verify.py`. Compare coefficient intervals to TH-004. The missing fast-mode lemma is intentionally not replaced by simulation.\n''')

# root handoff
w('HANDOFF.md',r'''
# FIN intrinsic-theory campaign handoff — TH-001–TH-005

Date: 2026-10-14.

## Disposition

The first intrinsic-theory package has produced one theorem-level advance and one sharply isolated critical-limit obligation:

1. **TH-002 path LDP:** Kraaij Theorem 2.3 applies to the actual leave-one-out q12 empirical chain for initial sequences with a good initial LDP. Deterministic starts and the declared Gibbs equilibrium starts are explicitly paid examples. The resulting action is the non-quadratic edge-flux relative-entropy action.
2. **Local quasipotential:** the exact reversal symmetry gives a basin-local `F(p)-F(min)` quasipotential along reversible relaxation orbits. This does not pay the global Gamma/capacity gate.
3. **TH-003 information result:** seven retained means are not path-sufficient; an exact N=120 count-state witness separates their fluctuation generator. Seven means plus four existing second moments form an exact 11D occupation lift.
4. **TH-004 microscopic fold bridge:** the actual jump generator gives soft noise coefficient exactly `2g_f`, not an assumed auxiliary Langevin covariance. The supplied interval gives variance coefficient `[7.0312894136791835,7.031289453679183]`.
5. **TH-005 stochastic fold:** the critical scaling and candidate SDE are derived at generator-characteristic level. Full stopped process convergence is **not** claimed. One precise lemma remains: fast-mode/hidden-mode averaging and containment on the accelerated time scale.

## Strongest justified scientific statement

The supplied FIN heat-bath model now has a rigorous interacting path large-deviation theory in its declared mathematical scope, plus a source-derived critical diffusion candidate whose drift/noise coefficients are fixed by the microscopic generator and certified fold. The next proof-grade move, if the architect continues this lane, is **only** the fast-mode averaging lemma needed to turn the stochastic-fold generator balance into a stopped convergence theorem.

## Architecture boundary

TH-006–TH-009 were not executed. TH-010 is reserved for architect review; this handoff does not choose a new physical programme. Empirical promotion remains stopped without independent input.
''')

j('CLAIM_REGISTER.json',{
"claims":[
 {"id":"TH2-LDP","status":"THEOREM_APPLICATION","claim":"path LDP for leave-one-out empirical process under good initial LDP","scope":"dimensionless q12 model"},
 {"id":"TH2-FLUX","status":"ANALYTIC_IDENTITY","claim":"nonquadratic edge-flux Lagrangian"},
 {"id":"TH2-QP","status":"BASIN_LOCAL_THEOREM","claim":"quasipotential F-Fmin along reversible basin relaxation"},
 {"id":"TH3-NONCLOSURE","status":"EXACT_FINITE_WITNESS","claim":"same mu can have different fluctuation generator"},
 {"id":"TH3-LIFT","status":"EXACT_ALGEBRA","claim":"mu plus four selected second moments recovers full p"},
 {"id":"TH4-NOISE","status":"DERIVED_IDENTITY_PLUS_SUPPLIED_INTERVAL","claim":"fold soft variance coefficient is 2g_f"},
 {"id":"TH5-LIMIT","status":"OPEN_ONE_LEMMA","claim":"candidate stopped diffusion with derived coefficients; convergence not yet proved"}
]})

w('NONCONCLUSIONS.md','''# Campaign nonconclusions\n\nNo result here derives A7, the heat-bath kinetic law, physical time, SI units, spatial incidence/locality, a propagating field or particle, quantum measurement/unitary dynamics, gravity, Standard Model structure, selector/QW-2191 closure, legacy-role transfer, or a theory of everything. The empirical lane remains stopped.\n''')

w('TH-010_FOR_ARCHITECT.md',r'''
# TH-010 packet for architect — no automatic decision

The executor stops after TH-001–TH-005 as instructed.

Questions for architect review:

1. Is the new path-LDP theorem sufficient to promote the intrinsic finite model to a controlled mesoscopic statistical theory? (It is not a physical-source theorem.)
2. Is paying `FAST_MODE_AVERAGING_AND_CONTAINMENT` the single best next mathematical task, or should the lane stop at the path action/nonclosure result?
3. If the stochastic-fold theorem is pursued, preserve the prepared/stopped local scope and do not add mobility/noise/clock parameters.
4. TH-006–TH-009 remain deferred unless separately selected.
''')

w('RESEARCH_SUMMARY_PL.md',r'''
# Podsumowanie badań TH-001–TH-005

Najmocniejszym nowym wynikiem jest pełne zastosowanie twierdzenia Kraaija o dużych odchyleniach do zadanej dynamiki FIN. Otrzymany koszt trajektorii nie jest przyjętą zasadą działania, lecz wynika z generatora skokowego i ma postać kosztu przepływów Poissona po przejściach etykiet.

Drugi ścisły wynik pokazuje, że siedem średnich modów nie wystarcza do przewidywania fluktuacji trajektorii: dwa prawdziwe stany licznikowe N=120 mają tę samą średnią, lecz różny generator kwadratu obserwabli. Cztery dodatkowe momenty drugiego rzędu wystarczają natomiast do dokładnego odtworzenia całego wektora zajętości p.

W obszarze folda wyprowadzono z rzeczywistego generatora współczynnik szumu miękkiego modu: jego wariancja na jednostkę modelowego czasu wynosi dokładnie `2 g_f/N`, więc po krytycznym skalowaniu amplituda graniczna to `sqrt(2 g_f)`. Nie jest to szum dopisany na poziomie makro.

Skalowanie `g-g_f~N^-2/3`, amplituda `N^-1/3` i czas `N^1/3` prowadzą formalnie do zatrzymanej dyfuzji

`dZ=-g_f(a lambda+b Z^2/2)ds+sqrt(2g_f)dB`.

Nie ogłaszam jednak pełnego twierdzenia zbieżności. Pozostał jeden konkretny lemat: trzeba rygorystycznie wyeliminować szybkie mody stabilne i ukryte składowe zajętości w problemie martyngałowym na przyspieszonym czasie. Diagnostyka skończonego N wspiera dokładnie przewidywane poprawki, ale nie zastępuje tego dowodu.
''')

w('REPLAY.md','''# Campaign replay\n\n1. `python theory_core.py > NUMERICAL_REPLAY.json`\n2. `python fold_diagnostics.py`\n3. `python independent_verify.py`\n4. Inspect TH-001 theorem table and TH-002 proof against the Kraaij Theorem 2.3 source URL.\n5. Check `MANIFEST.sha256`.\n''')

w('README.md','''# FIN intrinsic theory TH-001–TH-005\n\nThis package contains proof notes, numeric cross-checkers and a scoped handoff. Start with `HANDOFF.md`, `RESEARCH_SUMMARY_PL.md`, and `TH-002/GENERATOR_TO_ACTION_PROOF.md`.\n''')
