**# FIN after UF: local equilibrium and a genuine finite-N trapping window**



Date: 2026-10-14. Repository snapshot: \`21a875c5\` plus the current UF intake.

Status: **\*\*RESEARCH_PLAN; MT tasks below have not been executed.\*\***

Execution directory: \`fin_critical_metastability_campaign_20261014/\`.

Language: English for the execution agent.

Architect decision: [MT_DECISION.md]\(fin_uniform_critical_review_20261014/MT_DECISION.md).



**## 1. One objective**



Determine whether the actual iid-uniform-start FIN process has a \*\*proved

finite-N time interval\*\* after local relaxation and before sector escape,

on which conditional local equilibrium is an accurate effective law.



UF established finite-scaled-horizon weak convergence and the quartic limit's

invariant density. The intake established exponentially small full Gibbs mass

near u. Neither establishes local relaxation time or sector residence time

for actual finite-N trajectories. This campaign pays that missing comparison.



No new microscopic law, graph, force, noise, drive, storage or clock is added.

Reflection, conditional Gibbs and stopping regions are explicit \*\*analysis

objects\*\*, not physical controllers or secretly modified parent dynamics.



**## 2. Reuse and exact input registry**



\| Alias | Input |

\|---|---|

\| T01 | AGENTS.md, current UF/KI/SF admission and unchanged physical gates |

\| T02 | fin_uniform_critical_review_20261014/INTAKE_20261014.md, CLAIM_REGISTER.json, verification.json |

\| T03 | fin_uniform_critical_review_20261014/TRAPPING_PROOF_ROUTE.md |

\| T04 | fin_uniform_critical_review_20261014/GIBBS_EXCLUSION_PROOF.md, GIBBS_EXCLUSION_CERTIFICATE.json and checker |

\| T05 | fin_uniform_critical_campaign_20261014/UF-001_UNIFORM_CONTRACT_AND_INITIAL_PROOF.md and UF-002_PARITY_COEFFICIENTS_AND_CONTROLS.md |

\| T06 | fin_uniform_critical_campaign_20261014/UF-003_CRITICAL_PROCESS_THEOREM_AND_ENERGY.md and UF-004_QUARTIC_VS_FULL_GIBBS.md |

\| T07 | fin_kinetic_invariants_campaign_20261014/KI-002/REVERSIBILITY_AND_MOBILITY.md and KINETIC_RATE_CONTRACTS.json |

\| T08 | fin_rank7_mathphysics_next/proofs/MP7-031_finite_N_variational_model.md; fin_rank7_followup/inputs/fin_handoff_audit/results.json |

\| T09 | FIN_RESEARCH_AFTER_POST_AFTER_CONTINUATION_HANDOFF_131_255_20260927/01_REPORTS/LOCAL_CAPACITY_GATE_134.md, STRENGTHENED_RECURSIVE_LOCAL_CHILD_TEST_131.md and MEMORY_AWARE_METASTABLE_UNIT_VERDICT_137.md |

\| T10 | fin_post_326_methodology_review/REVIEW_20261014.md, especially disconnected reflected-preparation and unpaid B4 gates |

\| T11 | fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py; fin_generative_operational_campaign_20261014/RG-E04/KINETIC_CONTRACTS.json |



Bind actual hashes. Generic conductance/cut/path ideas are already present

in reports63/134 and must be cited as reused. The new target is the \*\*degenerate,

N-dependent critical window, iid reference coverage and finite-N local law\*\*.

Do not reopen the old parity-bit isolation ratios or global Gamma=B4 certificate.



**## 3. Frozen parent and target laws**



Primary: actual q12 constant-diagonal A=XX^T, unit-clock LOO HB, parent

normalization N, all-pairs Gibbs convention and iid-uniform initial labels.

The already declared E04-g0 MET/BAR laws are diagnostic controls with the same

Gibbs weights. Keep their prefactors unchanged. For fixed bounded lambda,



&#x20;   g_N=g_c+lambda N^(-1/2),  g_c=12/Lambda6,

&#x20;   sigma_j=(-1)^j, m=sigma^T p, p=p_m+w,

&#x20;   p_m=u+m sigma/12, w orthogonal to 1 and sigma.



Exact spectral g_c is used; a fixed rounded gain is not an asymptotic input.

For large N the gain lies in [5,6], so the intake's equilibrium exclusion

applies. It gives no exit rate and must not be used as one.



Target analysis objects:



&#x20;   D_N = a specified connected local count-state region near u,

&#x20;   tau_D = first actual parent exit from D_N,

&#x20;   pi_D,N = Gibbs law conditioned on D_N,

&#x20;   t_loc(N) = a genuinely proved local relaxation scale.



The strongest desired conclusion is a correctly scoped interval



&#x20;   t_loc(N) <= t <= exp(cN),

&#x20;   TV(Law_iid(p_N(t)),pi_D,N) -> 0 uniformly on that interval.



This conclusion is **\*\*inactive\*\*** unless both exit and local-mixing gates are

paid. Trapping or static conditional laws may be accepted separately if

mixing remains open. No matching escape exponent or mean-exit asymptotic is

an acceptance requirement.



**## 4. Task cards**



**### MT-001 — A paid local critical domain and its count graph**



**\*\*Scientific question:\*\*** is there a fixed local region with positive boundary

free energy, high iid-start coverage and a well-defined communicating component?

**\*\*Why it matters:\*\*** finite-horizon convergence does not provide a boundary

barrier; a visually connected continuous basin does not ensure a connected

reflected count chain or coverage of its whole conditioned distribution.

**\*\*Exact inputs:\*\*** T01–T11, especially the exact parity free energy and fast gap.

**\*\*Dependencies:\*\*** UF intake; no new root atlas or physical premise.

**\*\*Method:\*\***



1\. Derive \`F_g(p_m)\` exactly and fast transverse convexity on one fixed small

&#x20;  domain, not a fitted N-dependent Hessian sample. Pay



&#x20;      F_gN(p)>=c_m m^4+c_w||w||^2-C Lambda^2/N.



2\. Select one epsilon and positive h from this bound, not to obtain a desired

&#x20;  lifetime. A candidate outer cylinder is \`|m|\<epsilon, ||w||\<epsilon^2\`.

3\. Define D_N as the communicating component of its sublevel \`F_gN<=h\`

&#x20;  containing a balanced-uniform type. Prove the artificial outer cylinder

&#x20;  cannot provide a cheaper boundary. Keep actual parent jumps unchanged.

4\. Prove a fixed central type ball connects to the seed through genuine

&#x20;  permitted count moves of energy below h. Deduce iid-start coverage with

&#x20;  an explicit tail argument. If needed first restrict to N divisible by12;

&#x20;  label that scope and do not claim the all-integer extension for free.

5\. Identify the inner exit layer B_N and prove \`inf_B F_gN>=b>0\` for all

&#x20;  sufficiently large N, including jump overshoot/lattice errors.

6\. Record connectivity/seed independence or component restrictions. Do not

&#x20;  enumerate large count lattices to replace the argument.



**\*\*Deliverables:\*\*** LOCAL_DOMAIN_CONTRACT.json, COERCIVITY_AND_BOUNDARY_PROOF.md,

COUNT_CONNECTIVITY_AND_IID_COVERAGE.md, INPUTS.sha256.

**\*\*Acceptance criterion:\*\*** one actual domain, positive paid boundary height,

proved start coverage and connected/component-specific graph, in stated N scope.

**\*\*Failure/refutation criterion:\*\*** boundary leakage lacks the positive height,

initial mass misses the component, or connectivity is only assumed.

**\*\*What becomes possible:\*\*** an actual iid-start trapping bound and a valid

conditional local equilibrium object.

**\*\*Must NOT conclude:\*\*** a physical wall, autonomous small-N bit, global

communication height, canonical sign selection or spatial locality.

**\*\*Estimated computational class:\*\*** M, inequalities and tiny certificates only.

**\*\*Priority:\*\*** P0.



**### MT-002 — Actual-process trapping from invariant-weight domination**



**\*\*Scientific question:\*\*** how long is exit unlikely for the unchanged process

started iid uniformly, rather than in a postselected Gibbs state?

**\*\*Why it matters:\*\*** global equilibrium exclusion does not tell whether an

initial critical sector is short-lived or metastable.

**\*\*Exact inputs:\*\*** MT-001, T03/T07/T08/T11, exact type weights and actual rate caps.

**\*\*Dependencies:\*\*** paid exit layer, N/g and iid-coverage contract.

**\*\*Method:\*\***



1\. Define unnormalized invariant W_N relative to iid-uniform type law mu0_N.

&#x20;  Use all-pairs A>=0 and g_N>=0 to prove \`mu0_N<=W_N\` and semigroup domination

&#x20;  \`mu0_N P_t<=W_N\`. Do not call mu0_N stationary or Gibbs.

2\. Derive the compensator of actual D_N-to-outside jump counts. Include

&#x20;  initial outside mass, all departure states and the real clock.

3\. Pay \`sum rates<=C_rate N\` from each rule; e.g. a common coarse constant2

&#x20;  is available for the frozen HB/MET/BAR definitions.

4\. Combine with the type bound to prove a sufficient form



&#x20;      P_iid(tau_D<=T)

&#x20;        <=mu0_N(D_N^c)+C_rate N T (N+1)^12 exp(-Nb).



&#x20;  Every B_N/domain and boundary error must match MT-001; do not skip an

&#x20;  unknown normalization or silently drop a labelled outcome.

5\. For each fixed c\<b deduce survival up to exp(cN) in the admitted N scope.

&#x20;  Distinguish a one-sided high-probability lower lifetime from mean-time

&#x20;  asymptotics, an exact activation barrier or a matching upper escape bound.

6\. Reuse prior cut tools honestly. This new application must not borrow the

&#x20;  unpaid global B4 cover or assume all other escape channels are irrelevant.



**\*\*Deliverables:\*\*** IID_TRAPPING_THEOREM.md, RATE_AND_FLUX_PAYMENT.json,

INITIAL_EXCEPTION_BUDGET.json, NONCONCLUSIONS.md.

**\*\*Acceptance criterion:\*\*** a full parent-process exit inequality from the

specified initial law, or an exact unpaid boundary/initial/rate term.

**\*\*Failure/refutation criterion:\*\*** only equilibrium occupancy or a reflected

simulation is shown, a boundary channel is missing, or an exponent is fitted.

**\*\*What becomes possible:\*\*** a controlled long residence window independent of

a full global saddle identification.

**\*\*Must NOT conclude:\*\*** \`tau\~exp(Nb)\`, a capacity prefactor, persistent sign

memory, global two-state autonomy, physical seconds or a matter excitation.

**\*\*Estimated computational class:\*\*** S/M.

**\*\*Priority:\*\*** P0.



**### MT-003 — Conditional local Gibbs to joint quartic/Gaussian fluctuations**



**\*\*Scientific question:\*\*** does the actual conditional local equilibrium have

the same collective law as the now-proved UF limiting diffusion?

**\*\*Why it matters:\*\*** a stationary law of a limiting SDE is not automatically

the limit of any finite-N parent measure.

**\*\*Exact inputs:\*\*** MT-001, T05/T06/T08 and the unchanged count Gibbs law.

**\*\*Dependencies:\*\*** a connected local domain containing the relevant critical

types and a paid coercivity estimate; not MT-004 mixing.

**\*\*Method:\*\***



1\. Define pi_D,N with its exact multinomial weights and actual component.

&#x20;  Do not replace the count law by labelled weights without multiplicity.

2\. Expand in \`Z=N^(1/4)m\` and \`Y=sqrt(N)E_fast^T w\`. Derive the quartic and

&#x20;  fast quadratic form, with vanishing mixed corrections from the exact

&#x20;  symmetry/normalization. Identify the actual fast covariance.

3\. Pay lattice Jacobian/mesh and uniform interior Stirling errors, scaled

&#x20;  tail domination and any removed component/boundary contribution.

4\. Prove weak convergence and the moment uniform integrability actually

&#x20;  needed later. Pointwise Taylor expansion or gamma numerics is insufficient.

5\. Match the marginal Z law to \`exp(chi lambda z^2/2-z^4/12)\` and explain

&#x20;  why pi_D,N remains different from full Gibbs. A local partition prefactor

&#x20;  is optional unless a downstream estimate uses it; do not refine digits.



**\*\*Deliverables:\*\*** CONDITIONAL_GIBBS_LIMIT.md, JOINT_FLUCTUATION_MAP.json,

TAIL_AND_LATTICE_ERROR_PROOF.md, REQUIRED_MOMENT_BOUNDS.md.

**\*\*Acceptance criterion:\*\*** conditional finite-N measure convergence, including

all used moments, under the actual local domain and model convention.

**\*\*Failure/refutation criterion:\*\*** unaccounted lattice/conditioning tail,

component ambiguity or a surviving coupling changes the predicted law.

**\*\*What becomes possible:\*\*** quantitative local-relaxation tests and an honest

finite-model meaning for the quartic effective equilibrium.

**\*\*Must NOT conclude:\*\*** full global Gibbs quartic law, exact independence at

finite N, physical particles or a preserved-information substrate.

**\*\*Estimated computational class:\*\*** M, analytic finite-alphabet asymptotics.

**\*\*Priority:\*\*** P0.



**### MT-004 — Pay or isolate the local relaxation/time-window gate**



**\*\*Scientific question:\*\*** does the actual process approach pi_D,N before its

controlled exit, on a provably smaller scale?

**\*\*Why it matters:\*\*** local conditional weights and long survival do not alone

show local thermalization; this is the genuine missing effective-law bridge.

**\*\*Exact inputs:\*\*** MT-001–003, T09/T10 and exact reversible Dirichlet forms.

**\*\*Dependencies:\*\*** connected D_N, iid coverage, paid trapping and conditional moments.

**\*\*Method:\*\***



1\. Define an auxiliary reflected chain by suppressing only outgoing jumps,

&#x20;  retaining internal rates/clocks. Prove its invariant pi_D,N and exact

&#x20;  coupling with the unchanged parent until first exit. No physical wall.

2\. Use a Rayleigh test with m and MT-003's variance to establish the critical

&#x20;  upper-gap obstruction, expected order<=C N^-1/2. A ten-fast-mode gap is

&#x20;  not a gap for the whole reflected chain.

3\. Attempt one lower-gap/relaxation proof, e.g. a discrete two-scale Poincare

&#x20;  comparison or count-path congestion estimate. Pay conditional-fibre

&#x20;  connectivity and boundary effects; exact parity lumpability is unavailable.

4\. Any proved polynomial bound is useful; optimal constants are not required.

&#x20;  Prove the initial density/entropy bound for iid starts rather than assume

&#x20;  a worst-case-start or globally strongly stable mean-field theorem applies.

5\. If paid, derive a specified t_loc(N), mix/coupling/tail budget and the

&#x20;  uniform window to exp(cN), with c below the actual trapping bound.

6\. If not paid, export LOCAL_RELAXATION_OPEN with the exact comparison lemma,

&#x20;  residual/graph obstruction and scope of already proved partial results.

&#x20;  Do not convert UF weak convergence into a spectral-gap theorem or fit poles.



**\*\*Deliverables:\*\*** REFLECTED_COMPARISON_CONTRACT.md, LOCAL_RELAXATION_OR_BLOCKER.md,

FINITE_N_WINDOW_THEOREM.md only if proved, ERROR_AND_SCOPE_LEDGER.json.

**\*\*Acceptance criterion:\*\*** a full local relaxation and parent-transfer theorem,

or a precise genuinely open lemma; the latter does not admit the window claim.

**\*\*Failure/refutation criterion:\*\*** disconnected reflection, unpaid Gibbs-to-start

coverage, a fake constant gap, or a time interchange inferred from UF alone.

**\*\*What becomes possible:\*\*** a controlled finite-N effective statistical regime,

or identification of the exact remaining obstacle instead of another PASS.

**\*\*Must NOT conclude:\*\*** global mixing, matching escape asymptotics, B4 closure,

finite-N universal particle dynamics, physical units or independent source.

**\*\*Estimated computational class:\*\*** M, proof difficulty; no huge matrices/sweeps.

**\*\*Priority:\*\*** P1, after the three P0 tasks.



**### MT-005 — Scope handoff and architect stop**



**\*\*Scientific question:\*\*** which of long trapping, conditional equilibrium and

actual finite-N local relaxation is now established?

**\*\*Why it matters:\*\*** these three statements are logically different, and a

long-lived sector need not retain its internal sign information.

**\*\*Exact inputs:\*\*** all MT results/failures and T01/T02/T09/T10.

**\*\*Dependencies:\*\*** terminal dispositions for preceding gates, not forced success.

**\*\*Method:\*\***



1\. Export a DAG from supplied model/initialization to coercivity/domain,

&#x20;  trapping, conditional law, relaxation and actual-process window.

2\. Distinguish exact identities, inherited input assumptions, analytic

&#x20;  asymptotics, numeric diagnostics, reused capacity tools and open estimates.

3\. Record sector residence and sign-memory/critical relaxation separately.

&#x20;  Do not identify either with a physical particle or primordial clock.

4\. Preserve failures and all source/input hashes, replay only tiny relevant

&#x20;  checks, and return one exact NEXT_OPEN_ATOM if the relaxation gate remains.

5\. Stop. Architect chooses any later global/source/empirical direction; no

&#x20;  automatic hierarchy, larger N/g sweep or new preparation controller.



**\*\*Deliverables:\*\*** HANDOFF.md, REPORT.md, CLAIM_REGISTER.json, NONCONCLUSIONS.md,

NEXT_OPEN_ATOM.md, REPLAY.md, MANIFEST.sha256 and proposed AGENTS patch.

**\*\*Acceptance criterion:\*\*** reproducible scoped results with unpaid windows

explicitly unpromoted; no physically expanded interpretation of local bounds.

**\*\*Failure/refutation criterion:\*\*** a lower trapping bound is reported as an

exact mean lifetime or a conditional law as attained parent equilibrium.

**\*\*What becomes possible:\*\*** one evidence-based next architecture decision.

**\*\*Must NOT conclude:\*\*** fundamental kernel/clock/initial-source, locality,

conservation/EV, selector/legacy-role, QM/GR/SM/L_total/ToE closure.

**\*\*Estimated computational class:\*\*** S.

**\*\*Priority:\*\*** P0.



**## 5. Execution limits**



Start MT-001, MT-002 and MT-003. MT-004 is a separate proof gate, not something

to replace by more simulations. No large state enumeration, full saddle

census, broad N/g sweep, new source inventory or new experiment.



Use \`fin_critical_metastability_campaign_20261014/\` only. Preserve AGENTS,

prior sources and reviews; return a proposed patch, not an automatic promotion.

Per task export REPORT.md, RESULTS.json, PROOF_OR_OPEN_LEMMA.md,

NONCONCLUSIONS.md, REPLAY.md and INPUTS.sha256. RESULTS separates execution,

science verdict, supplied/derived/fitted objects, initial/domain/time scopes,

arithmetic status, input dependencies, failures and next action.



If using a literature theorem, read it fully and pay its actual hypotheses.

Global unique-exponential-stability results do not apply here automatically.

No external spending, hardware, messages, FIN-M0/task335, selector/FAR replay

or quantum/gravity branch is authorized. No new primitive is allowed to repair

a failed estimate under the old model identifier.



**## 6. Kill tests and value for physics**



Kill the claimed window if boundary coercivity fails, the iid start is not

covered, reflected components are mixed without access, or local relaxation

cannot be separated from exit. Keep a correct partial trapping/static theorem

if the dynamic window remains open. A tiny conservative exponent alone does

not establish practical physical memory at accessible N.



Success would give an actual finite-N local statistical regime and an honest

time-scale hierarchy. It still would not derive physical space, propagation,

fundamental kinetics/gain, natural initial conditions, units or measured reality.
