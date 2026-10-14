**# Intake of FIN intrinsic TH-001–TH-005**



Date: 2026-10-14. Repository baseline: \`ec8baa15\`.

Source: [HANDOFF.md]\(../fin_intrinsic_theory_campaign_20261014/HANDOFF.md).

Predecessor: [theory synthesis]\(../fin_theory_synthesis_20261014/HANDOFF.md).

Plan: [intrinsic-theory campaign]\(../FIN_INTRINSIC_THEORY_NEXT_CAMPAIGN_20261014.md).



**## Verdict**



\*\*ACCEPT_INTERACTING_PATH_THEORY_AND_MICROSCOPIC_FOLD_COEFFICIENTS;

STOCHASTIC_PROCESS_CONVERGENCE_REMAINS_OPEN.\*\*



The specifically requested synthesis is already tracked at the baseline.

Its 13 manifest entries still agree. A new, initially untracked TH execution

package is present alongside it. This intake reviews that actual execution

as the continuation; it does not relabel the historical synthesis as having

performed TH or merge identical old findings a second time.



The chief new result is the application of a published interacting path-LDP

theorem to the already declared leave-one-out empirical chain. The finite

information witness and much of the coefficient algebra were present in the

synthesis: the new execution makes their conventions and implications more

explicit, but they are not all newly discovered results.



**## Provenance, dependencies and checks**



All **\*\*80/80\*\*** incoming manifest entries agree, including one producer \`.pyc\`

cache artifact. There are **\*\*79 non-cache entries\*\***. All eight supplied source

input entries and all 35 stage input entries agree. The cache is not scientific

evidence and need not be versioned or executed; its optional absence in a clean

checkout is reported separately by the intake verifier. No scientific source

file is changed. The two embedded predecessor ZIPs are hash-bound context,

not evidence of the executor's Git checkout or remote reading history.



The attached execution plan is not byte-identical to the repo plan. Removing

Markdown escaping/emphasis, HTML entities and whitespace expansion makes

their text agree. Both hashes are retained. The older attached AGENTS snapshot

does not establish that the executor had the live guardrail file.



The source scientific scripts and its package checker ran in a temporary copy.

\`independent_verify.py\` passed **\*\*41/41\*\***. In the current numerical environment

the three main regenerated outputs are **\*\*not byte-identical\*\*** to the incoming

outputs. Maximum coordinate drift in the 24-row generator diagnostic is

approximately \`7.44e-14\`. This admits numerical reproducibility, not a new

cross-environment byte-replay claim. \`write_reports.py\` was inspected where it

constructs contracts/coefficient/status narratives, but not run: it rewrites

scientific reports and hardcodes interval endpoints.



The separate intake [verifier]\(verify_intake.py) uses the versioned safe-seed

model, a different fold equation (smallest Hessian eigenvalue rather than a

determinant), the full occupation linearization and independent interval

scalar propagation. It passes **\*\*22 bounded checks\*\***. It also replays the

existing small R7P-031 interval root/inertia calculation, not a global census.



The actual local fold certificate hashes to the source's stated

\`5537174281996415a20ea3dffce23c33258c51396b3850aedf4168d9280e352e\`.

Independent interval evaluation places the exact-decimal kernel's four active

spectral values inside that certificate's accepted rational input boxes.

This closes a useful **\*\*local input binding\*\***, not a physical derivation of A7.

Interval checks retain the mpmath.iv correctness assumption. Floating spectrum,

cumulant and characteristic checks are not outward stochastic error bounds.



The full statement of [Kraaij's Theorem 2.3]\(https\://arxiv.org/html/1502.05185v3#S2.Thmtheorem3)

was freshly consulted. The captured text hash is in \`verification.json\`.

This access verifies the theorem used here; it does not verify the producer's

historical web access. The theorem's published title includes "via the

comparison principle for an associated Hamilton–Jacobi equation".



**## Admission by task**



**### TH-001 — model and dependencies**



Accept q=12, the supplied exact-decimal Fourier construction A=XX^T, the

Gibbs extension and the unit-attempt-rate leave-one-out heat bath as the

declared mathematical parent. A7, gain convention, kinetics and attempt clock

remain supplied model premises. The local input binding above improves the

source's missing-local-provenance status without retrospectively attesting

its remote environment. Existing static/ODE/CLT material is reused.



**### TH-002 — interacting path theory**



Accept the theorem application for each fixed finite gain (or its appropriately

bounded parameter setting), with an initial empirical-law LDP. The limiting

edge intensity is \`c_ij(p)=p_i softmax_j(gAp)\`. The uniform finite-N rate error

is at most \`|g| max_k|A_ki|/(2N)\`; the absolute value is needed if gains of

either sign are included. The declared physical-interest lane uses positive g.

Continuity, source-positive intensity and the monotone boundary factorization

\`p_i \* s_j(p)\` pay Theorem 2.3's actual hypotheses. Deterministic convergent

type starts have the trivial initial LDP; Gibbs starts have the already

available method-of-types rate \`F_g-min F_g\`.



The resulting Skorokhod path rate includes the **\*\*initial cost\*\*** and the time

integral of the Legendre Lagrangian, and is infinite outside absolutely

continuous paths. Its edge-flux representation is



&#x20;   L(p,v)=inf\_{q>=0, sum q_ij(e_j-e_i)=v}

&#x20;          sum_ij [q_ij log(q_ij/c_ij)-q_ij+c_ij].



Use zero-rate conventions and the simplex tangent cone at its boundary.

Self-refresh edges may be omitted: they have zero displacement and optimize

at their reference flux, contributing zero cost.



Accept the interior identity \`H(p,xi)=H(p,grad F-xi)\` and its consequence

\`L(p,v)=grad F.v+L(p,-v)\`. For an interior **\*\*nondegenerate stable minimum\*\***

and targets whose interior relaxation tends to it, time reversal gives the

basin-local quasipotential \`F(p)-F(m)\`. The small connection from m to the

truncated reversed orbit has vanishing cost: locally strictly positive rates

make the tangent Lagrangian smooth with a bounded quadratic comparison near

zero drift. This is an infimum over durations, not generally a finite-time

minimizer starting exactly at m. Boundary targets and degenerate minima require

separate treatment. Nothing here pays global \`Gamma=B4\`, capacity or prefactors.



Fixed-time path LDP is not the stochastic-fold theorem on a growing time

horizon and parameter sequence. It supplies no physical action or energy.



**### TH-003 — occupation information**



Accept the exact N=120/g=0 equal-mean witness with squared-feature generator

difference \`Lambda5/7200>0\`; it is an existing synthesis result replayed here.

Accept the constructive eleven-coordinate lift: seven means and four specified

second moments recover the occupation vector by Fourier inversion. Its

dimension-minimality concerns \*\*continuous injective recovery of the open

eleven-dimensional occupation simplex\*\*, not arbitrary finite-set coding or

the fewest variables for one particular observable or asymptotic law.



The continuous path contraction retains all compatible hidden paths and their

initial costs. No general mean-only instantaneous Markov action follows.

Recovering p does not recover labelled identities or introduce a new substrate.



**### TH-004 — microscopic fold coefficients**



Accept the generator-derived leading covariance and the conditional algebraic

identity \`v^T D_theta(p_f)v=2g_f\`. At finite N the variance rate is

\`2g_f/N+O(N^-2)\` at a corresponding near-fold type state, **\*\*not exactly\*\***

\`2g_f/N\` for the leave-one-out chain. This qualifies the wording of the source

Polish summary. Auxiliary Gaussian representation is not actual path dynamics.



The fold coefficients a,b are reused from R7P/MP7, not newly fitted. The

leading leave-one-out soft drift correction is \`-(v.theta_f)/N\`, becoming

\`O(N^-1/3)\` after critical rescaling. The nonzero fourth jump cumulant

approximately \`120.5038353\` is a numerical local diagnostic of the tilted

limiting jump law, not an additional conserved observable.



The source's decimal scalar endpoints are hardcoded floating representations.

Rather than promoting their last digits as outward, the intake interval

propagation exports conservative **\*\*outward decimal strings\*\***, including:



\| Object | Intake enclosure |

\|---|---|

\| \`2g_f\` | \`[7.03128941367918, 7.03128945367919]\` |

\| \`sqrt(2g_f)\` | \`[2.65165786135375, 2.65165786889621]\` |

\| \`c=(4/b)^(1/3)\` | \`[3.22744435943171, 3.22744562214155]\` |

\| \`d=c^2/(2g_f)\` | \`[1.48143482953572, 1.48143599716089]\` |

\| \`kappa=-ac/2\` | \`[0.34303154768747, 0.34303170459030]\` |



These depend on the accepted fold, spectral tuple, coordinate normalization

and interval backend; they are not measured constants.



**### TH-005 — open process convergence**



Accept the formal generator scaling and source-derived candidate



&#x20;   g_N=g_f+lambda N^(-2/3),

&#x20;   Z_N=N^(1/3) v.(theta_ret(N^(1/3)s)-theta_f),

&#x20;   dZ=-g_f(a lambda+b Z^2/2)ds+sqrt(2g_f)dB.



The 24 characteristic checks are diagnostics, not trajectories or convergence

proofs. The candidate's Gaussian increments do not make its nonlinear state

law Gaussian. Canonicalization is correct conditional algebra, not a proved

cross-model universality statement.



\*\*Do not accept the wording "only one precise lemma" as a completed theorem

contract.\*\* It packages several obligations: prepared initial law, stable-mode

occupation bounds, full occupation/covariance control, tightness, identification

and stopped martingale-problem uniqueness. Source "prepared sufficiently

controlled" is not a quantitative initial class. Stopping at a fixed original

theta-neighbourhood alone leaves \`|Z_N|\` growing as \`N^(1/3)\` and does not

license uniform compact-z Taylor bounds or global explosion/exit predictions.



The next plan supplies a **\*\*target contract\*\*** on fixed \`|Z|\<R\` and finite scaled

horizon S, with a shrinking fast tube and a provable negligible tube-exit event.

This is a proposed scope repair, not a retroactive source result. Neither

critical convergence nor global mixing/escape is promoted now.



**## Additional bounded architectural result**



From \`A=XX^T\`, let R project onto the four hidden occupation modes. Exactly

at the mean-field level,



&#x20;   r=Rp,   r_dot=R softmax(X theta)-r.



The retained drift is independent of r, so the full tangent linearization at

the fold is block triangular. Its eigenvalues are the seven retained rates

\`-g_f eig(H7)\` and four further rates \`-1\`. Together with the accepted inertia

this proves one soft and **\*\*ten stable\*\*** occupation directions for this supplied

model. It is not a uniform finite-N stochastic hidden-gap theorem. A full

numerical eleven-mode check independently confirms the algebraic structure.



The [proof route]\(FAST_MODE_PROOF_ROUTE.md) uses that derived stable structure

to propose a stopped Lyapunov estimate. Its stochastic estimates remain

targets; no new mobility, memory law or external fast environment is introduced.



**## Architectural outcome and limits**



[TH_010_DECISION.md]\(TH_010_DECISION.md) selects

**\*\*CLOSE_PREPARED_COMPACT_STOPPED_STOCHASTIC_FOLD\*\***. Finish this missing bridge

before optional phase-limit/universality tasks. One conditional first-passage

prediction is useful only after convergence, not another fit or large sweep.



The scientific potential has increased specifically toward a \*\*derived

mesoscopic statistical theory\*\* of the declared mean-field model. A7/kinetic

source, locality/incidence, physical units, conserved physical energy and the

operational/empirical interface remain missing. No spacetime, matter, QFT,

QM/GR, selector/role-transfer, Standard Model or ToE claim is admitted.
