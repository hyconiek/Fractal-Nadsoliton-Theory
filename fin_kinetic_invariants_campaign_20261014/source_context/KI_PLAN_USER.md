**# FIN after SF: kinetic invariants of a derived critical collective law**



Date: 2026-10-14. Baseline: \`cf1bf9ae\` plus the current SF intake.

Status: **\*\*RESEARCH_PLAN; KI tasks below have not been executed.\*\***

Execution directory: \`fin_kinetic_invariants_campaign_20261014/\`.

Language: English for the execution agent.

Architect decision: [KI_DECISION.md]\(fin_stochastic_fold_closure_review_20261014/KI_DECISION.md).



**## 1. Objective and boundaries**



Derive which parts of the admitted stochastic-fold law depend only on the

existing FIN equilibrium landscape and which depend on its specified kinetics.

Prove or refute a local kinetic-universality statement using already declared

reversible alternatives, with their actual mobility and fluctuations.



Primary FIN remains the supplied q12 A7 + exact LOO heat bath + unit attempt

rate per labelled copy. Diagnostic alternatives are not newly discovered FIN

laws or sources of nature's dynamics. No graph, field, bath, drive, clock,

physical scale, complex amplitude, extra conserved quantity or measurement

device is inserted. Model-family comparison is mathematical diagnosis.



The preceding SF result is prepared, compact-stopped and asymptotic. It does

not give quantitative finite-N TV error, arbitrary initial preparation,

unbounded mean-time convergence or global phase-switching/capacity control.

Do not reopen B4 or attempt to infer physical time from a normalization.



**## 2. Exact inputs and deduplication**



\| Alias | Repository input |

\|---|---|

\| K01 | AGENTS.md, current SF/TH admission and unchanged physical gates |

\| K02 | fin_stochastic_fold_closure_review_20261014/INTAKE_20261014.md, STOPPED_LIMIT_COMPLETION.md, CLAIM_REGISTER.json, verification.json |

\| K03 | fin_stochastic_fold_closure_campaign_20261014/SF-001/FAST_CHART_PROOF.md |

\| K04 | fin_stochastic_fold_closure_campaign_20261014/SF-002/GENERATOR_ENERGY_PROOF.md and SF-003/STOPPED_CONVERGENCE_THEOREM.md |

\| K05 | fin_stochastic_fold_closure_campaign_20261014/SF-004/FIRST_PASSAGE_PREDICTION.md and EXIT_SIDE_CERTIFICATE.json |

\| K06 | fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py |

\| K07 | fin_generative_operational_campaign_20261014/RG-E04/KINETIC_CONTRACTS.json and ERROR_BUDGET.md |

\| K08 | fin_physical_bridge_campaign_01/PHYS-004/REPORT.md and CALIBRATION_CONTRACT.json |

\| K09 | fin_rank7_followup/certificates/R7P-031_simple_fold.json and accepted spectral inputs; fin_rank7_mathphysics_next/proofs/MP7-025_simple_fold_normal_form.md |

\| K10 | FIN_RESEARCH_AFTER_POST_AFTER_CONTINUATION_HANDOFF_131_255_20260927/01_REPORTS/MICROSCOPIC_PROCESS_CONTRACT_138.md and RECURSIVE_LOCAL_CHILD_TEST_131.md |

\| K11 | fin_intrinsic_theory_campaign_20261014/TH-002/GENERATOR_TO_ACTION_PROOF.md and TH-003/PATH_CONTRACTION.md |

\| K12 | fin_theory_synthesis_20261014/LITERATURE_REGISTER.json, for primary theorem conditions only |



Bind actual hashes. Reuse the matching fold, heat-bath SF proof, prior mean-field

Lyapunov/Onsager identities and existing N2 kinetic contrasts; they are not

new discoveries. New targets concern a **\*\*general-N critical comparison\*\***, its

kinetic-adapted readout and invariants. No full source inventory is needed.



K07 and K08 use different stated clock conventions. Use \*\*K07's g0 equal

single-copy rule\*\* for this comparison; do not silently mix K08's N2-sector

normalization into primary unit-attempt FIN.



**## 3. The fixed comparison families**



For a count departure label i, define



&#x20;   h_j=(g/N)[A(n-e_i)]\_j,   delta_ij=h_j-h_i.



For i!=j the declared, g0-matched per-copy rates are



&#x20;   HB: r_ij=softmax_j(h),

&#x20;   MET: r_ij=(1/12) min(1,exp(delta_ij)),

&#x20;   BAR: r_ij=(2/12)/(1+exp(-delta_ij)).



The count rate is \`n_i r_ij\`; diagonal generator entries are minus actual

off-diagonal row sums. Self-attempt conventions do not authorize a new clock.

The MET/BAR definitions follow K07's distinct-label attempt rates 11/12 and

22/12. Freeze their extension to general N explicitly and prove reversibility

against the existing constant-diagonal Gibbs law; do not infer it from N2 data.

They are conditional diagnostic models, not replacements of primary HB.



**## 4. Useful proof targets, not already admitted cross-kinetic results**



Let \`F_g(p)=sum p_i log(12p_i)-(g/2)p^T A p\`. At the exact fold use

\`r=Sigma_f Xv\`, with the source-fixed normalization. Target:



&#x20;   H_F r=constant\*1 on the ambient space,

&#x20;   a_F=r^T partial_g grad F=a,

&#x20;   b_F=D^3 F[r,r,r]=b.



The last two identities should follow from

\`X^T r=v/g_f\` and \`r_i=p_f,i[(Xv)\_i-E_p_f(Xv)]\`.

A constrained Hessian, not the raw ambient matrix, is the relevant object.



For leading edge intensities \`c_ij=p_i r_ij(p)\`, reversibility should give



&#x20;   log(c_ij/c_ji)=partial_i F-partial_j F.



The logarithmic-mean conductance yields an exact interior gradient structure;

at equilibrium its mobility is



&#x20;   K_f=sum\_(i\<j) c_ij(p_f)(e_j-e_i)(e_j-e_i)^T,

&#x20;   D_f=2K_f.



On the tangent T, K_f is positive if all declared edge rates are positive.

The candidate kinetic soft clock is



&#x20;   m_K=[r^T(K_f|\_T)^(-1)r]^(-1),

&#x20;   ell_K^T=m_K r^T(K_f|\_T)^(-1),

&#x20;   dZ_K=-m_K(a lambda+b Z_K^2/2)ds+sqrt(2m_K)dB.



This inverse and left functional are restricted to T. A pseudoinverse must

use the same normalization gauge; an ordinary inverse of singular K_f is invalid.

HB has K_f=Sigma_f and should give m_K=g_f. Coefficients are derived from

the microscopic rule, not trajectory-fitted mobilities.



The coordinate \`ell_K^T(p-p_f)\` may differ from HB's thermal soft readout.

Do not assume a coordinate swap is harmless at path level. A promising extra

bound for rounded preparations is a stopped fourth-moment inequality for

\`V^2\`, giving



&#x20;   P(sup_before_exit ||w_N|| > epsilon N^(-1/3)) -> 0.



The preceding second-moment tube bound alone does **\*\*not\*\*** prove this assertion.



**## 5. Task cards**



**### KI-001 — Primal/dual thermodynamic fold matching**



**\*\*Scientific question:\*\*** are the soft coefficients a,b the same landscape

objects in the actual occupation free energy and the dual fold potential?

**\*\*Why it matters:\*\*** without matching chart/normalization, kinetic comparison

can introduce fictitious new couplings or a fitted canonical scale.

**\*\*Exact inputs:\*\*** K01–K06, K09–K11.

**\*\*Dependencies:\*\*** accepted SF intake and exact implicit fold, not midpoint data.

**\*\*Method:\*\***



1\. Freeze primary HB, comparison purpose, exact F/Phi/r/v and tangent gauge.

2\. Prove constrained Hessian nullity/inertia, using the accepted H7 certificate

&#x20;  and the entropy covariance congruence. The four extra occupation directions

&#x20;  must be included, not discarded by retained-mean closure.

3\. Derive \`a_F=a\` and \`b_F=b\` with the exact source-oriented r normalization.

4\. Explain why local transverse elimination does not alter the leading cubic

&#x20;  coefficient, while higher terms can depend on the reduction chart.

5\. Export which content is reused primal/dual theory and what is newly matched.

&#x20;  Re-run only a tiny arithmetic/input check whose hashes changed.



**\*\*Deliverables:\*\*** MODEL_AND_COMPARISON_CONTRACT.json, THERMODYNAMIC_MATCHING.md,

HESSIAN_AND_COEFFICIENT_MAP.json, INPUTS.sha256.

**\*\*Acceptance criterion:\*\*** exact constrained identities and leading matching

with all normalizations fixed; no numerical a/b refit.

**\*\*Failure/refutation criterion:\*\*** a mismatch survives correct coordinates or

a required thermodynamic object is supplied rather than derived.

**\*\*What becomes possible:\*\*** common landscape coefficients across kinetic controls.

**\*\*Must NOT conclude:\*\*** conserved physical energy, units or a source of A7.

**\*\*Estimated computational class:\*\*** S/M.

**\*\*Priority:\*\*** P0.



**### KI-002 — Derive the reversible mobility and fluctuation matrix**



**\*\*Scientific question:\*\*** what does each already declared microscopic kinetic

rule change at the same equilibrium fold?

**\*\*Why it matters:\*\*** equilibrium probabilities do not fix rates; arbitrary

macro mobility would defeat the purpose of a source-derived collective theory.

**\*\*Exact inputs:\*\*** KI-001, K06–K08 and the finite/limiting generator formulas.

**\*\*Dependencies:\*\*** exact frozen comparison contract.

**\*\*Method:\*\***



1\. Prove finite-N Gibbs reversibility for all three explicit count rules;

&#x20;  separate actual finite-N and limiting edge intensities.

2\. Derive the interior force-ratio and logarithmic-mean mobility directly.

&#x20;  Reuse existing Onsager theory where applicable, but pay every new rate law.

3\. At p_f prove \`D_f=2K_f\`, tangent positivity and \`J=-K_f H_F\`.

4\. Establish one null/ten stable directions via the appropriate constrained

&#x20;  congruence, and derive r, ell_K and m_K in one gauge.

5\. Recover HB's \`K_f=Sigma_f\`, \`m_K=g_f\`. Compute only nominal MET/BAR mobility

&#x20;  comparisons if useful; mark numerics separately from analytic identities.

6\. Retain the fixed g0 convention, invalid denominators and controller distinctions.

&#x20;  A clock calibrated by fitting the final fold curve is forbidden.



**\*\*Deliverables:\*\*** KINETIC_RATE_CONTRACTS.json, REVERSIBILITY_AND_MOBILITY.md,

DERIVED_SOFT_CLOCKS.json, NEGATIVE_CONTROLS.md.

**\*\*Acceptance criterion:\*\*** drift/noise clocks follow from the actual declared

rules, with tangent positivity and matching HB recovery.

**\*\*Failure/refutation criterion:\*\*** reversibility, positivity or matching fails,

or rates require an extra bath/clock parameter after the fact.

**\*\*What becomes possible:\*\*** a proved or falsified critical comparison candidate.

**\*\*Must NOT conclude:\*\*** all reversible kinetics are physically equivalent, or

that Gibbs equilibrium derives the actual fundamental kinetic rule.

**\*\*Estimated computational class:\*\*** S/M, at most 12/11/10-dimensional matrices.

**\*\*Priority:\*\*** P0.



**### KI-003 — Process theorem and kinetic-versus-thermal observable gate**



**\*\*Scientific question:\*\*** do the controls have the same canonical critical law,

and does it govern the same observable rather than only a convenient coordinate?

**\*\*Why it matters:\*\*** normalization may hide kinetic dependence; mean closure

and a left-null coordinate do not automatically coincide outside HB.

**\*\*Exact inputs:\*\*** KI-001/002, K03/K04, and the intake's stopped-limit completion K02.

**\*\*Dependencies:\*\*** derived positive mobility and full stable restriction.

**\*\*Method:\*\***



1\. Define \`x_K=ell_K^T(p-p_f)\` and its exact theta-readout relation; use the

&#x20;  same fold right vector for actual rounded preparations.

2\. Derive the projected drift/noise candidate with m_K; compare with SF's HB

&#x20;  theorem. Do not copy HB's closed seven-mean drift to MET/BAR.

3\. Check regularity of the actual rates: MET has equal-logit kinks. Prove

&#x20;  sufficient locally Lipschitz flux/mobility estimates multiplied by the

&#x20;  vanishing fold force, or leave the MET theorem open. No fake C3 premise.

4\. Adapt actual-generator fast energy and integrated remainder estimates;

&#x20;  retain finite-N corrections, auxiliary stop and exact root uncertainty.

5\. Pay the stopped martingale/clock absorption argument for each admitted lane.

&#x20;  Transfer a general theorem only after every hypothesis is explicitly checked.

6\. Prove a strong enough maximal fast bound to compare x_K with the actual

&#x20;  thermal readout at N^(1/3) amplitude scale. One route is \`L_N V^2\` with

&#x20;  every quadratic/higher jump term. Rounded starts have the required higher

&#x20;  moment; do not silently extend that to arbitrary second-moment-only starts.

7\. Return either a full same-observable canonical theorem, a coordinate-only

&#x20;  theorem with explicit mismatch, or the exact surviving error/regularity gate.



**\*\*Deliverables:\*\*** KINETIC_CRITICAL_LIMIT.md, METROPOLIS_REGULARITY_GATE.md,

OBSERVABLE_TRANSFER_PROOF.md, ESTIMATE_LEDGER.json and stopped-scope contracts.

**\*\*Acceptance criterion:\*\*** process/readout/absorption transfer is paid analytically,

or its failure is precisely isolated. Generator balance is not acceptance.

**\*\*Failure/refutation criterion:\*\*** an observable correction survives, a kink

destroys the bound, another soft direction appears or mobility/noise is fitted.

**\*\*What becomes possible:\*\*** a scoped kinetic-universality or nonuniversality result.

**\*\*Must NOT conclude:\*\*** arbitrary kinetic universality, global mixing/escape,

new Markov closure at all N, physical time or empirical validation.

**\*\*Estimated computational class:\*\*** M, proof work without trajectory sweeps.

**\*\*Priority:\*\*** P0.



**### KI-004 — Separate universal response from kernel and kinetic fingerprints**



**\*\*Scientific question:\*\*** what experimentally meaningful content could remain

after the lawful normalizations, and what no longer discriminates FIN?

**\*\*Why it matters:\*\*** a universal fold number is useful class physics but cannot

independently select a particular A7 or physical ontology.

**\*\*Exact inputs:\*\*** accepted KI-003 scopes, K05, fixed comparison clocks K07.

**\*\*Dependencies:\*\*** only controls whose process/readout theorem passed.

**\*\*Method:\*\***



1\. Canonicalize with \`c=(4/b)^(1/3)\`, \`d_K=c^2/(2m_K)\` and \`kappa=-ac/2\`.

&#x20;  Transform initial law, boundaries and observation clock together.

2\. Classify predicted exponents, exit-side laws, time factors and gain map:

&#x20;  universal after declared normalizations / landscape-dependent /

&#x20;  kinetic-dependent / unresolved. Different categories must remain distinct.

3\. Reuse the existing rational ±1/lambda0 certificate. If the theorem makes

&#x20;  it identical in the controls, report a **\*\*non-specificity result\*\***, not a

&#x20;  FIN-selected physical constant or the need for more decimal precision.

4\. Identify at most one dimensionless contrast that does not disappear under

&#x20;  the declared lawful clock/coordinate normalization. Specify its assumptions.

&#x20;  If none is justified, return NO_NEW_FIN_SPECIFIC_CONTRAST.

5\. State what real preparation/readout/source information would eventually be

&#x20;  needed to use this as physics, without inventing records or activating a lab.



**\*\*Deliverables:\*\*** INVARIANT_AND_FINGERPRINT_TABLE.md, SPECIFICITY_VERDICT.json,

ONE_CONTRAST_OR_NO_GO.md, NONCONCLUSIONS.md.

**\*\*Acceptance criterion:\*\*** theorem-backed invariance/non-specificity or one

honest surviving contrast; no numbers promoted by analogy or normalization.

**\*\*Failure/refutation criterion:\*\*** invariance depends on retuning controller,

kernel, boundary or observable after looking at the result.

**\*\*What becomes possible:\*\*** architect choice of a source-independent robustness

question or identification of the exact indispensable physical premise.

**\*\*Must NOT conclude:\*\*** natural kernel selection, universal constants, units,

locality, quantum/gravity or a generated physical world.

**\*\*Estimated computational class:\*\*** S/M.

**\*\*Priority:\*\*** P1, dependent on KI-003.



**### KI-005 — Scope-aware handoff and architect checkpoint**



**\*\*Scientific question:\*\*** what is now kinetic-invariant, and which physical

premises cannot be removed by mathematical universality?

**\*\*Why it matters:\*\*** useful negative controls must not be turned into another

unending sequence of toy-model successes or nonexistent source closures.

**\*\*Exact inputs:\*\*** terminal KI-001–004 artifacts, K01/K02 and preserved source packages.

**\*\*Dependencies:\*\*** preceding tasks have scoped terminal dispositions; failed

KI-003 controls cannot be admitted through a KI-004 comparison table.

**\*\*Method:\*\***



1\. Build the dependency graph from supplied rates/equilibrium to mobility,

&#x20;  full chart, stochastic bounds, observable transfer and canonical law.

2\. Separate analytic theorem, conditional premise, interval result, numerical

&#x20;  diagnostic and open/no-go. Mark every reused SF/TH/MP7 result.

3\. Preserve all adverse or inconclusive controls, source hashes and replay.

4\. Summarize whether the supplied kinetic law is only a clock choice in this

&#x20;  critical scope or leaves a further observable dependence.

5\. Stop and return to the architect. Do not automatically alter A7, add a

&#x20;  spatial composition rule, launch an empirical campaign or extend all horizons.



**\*\*Deliverables:\*\*** HANDOFF.md, REPORT.md, CLAIM_REGISTER.json, NONCONCLUSIONS.md,

NEXT_OPEN_ATOM.md, REPLAY.md, MANIFEST.sha256 and proposed AGENTS patch.

**\*\*Acceptance criterion:\*\*** one small reproducible package with complete proof

and negative-result scope; no automatic next campaign.

**\*\*Failure/refutation criterion:\*\*** reporting counts replace proofs, controls

become new FIN axioms or numerical comparisons masquerade as process theorems.

**\*\*What becomes possible:\*\*** one evidence-driven next architecture decision.

**\*\*Must NOT conclude:\*\*** full finite-model uniqueness, empirical or fundamental closure.

**\*\*Estimated computational class:\*\*** S.

**\*\*Priority:\*\*** P0.



**## 6. Execution budget and rules**



\- Start KI-001 -> KI-002 -> KI-003. KI-004 is conditional, then KI-005 handoff.

\- Work only in the new execution directory. Do not edit AGENTS or history;

&#x20; export a proposed accepted-claim patch for intake.

\- Primarily analytic finite-dimensional work; no large-N trajectories, full

&#x20; count-state enumeration, broad gain/kernel sweep or static-root refinement.

\- Tiny direct-rate/matrix checks are allowed. Certify arithmetic only when

&#x20; claiming numerical enclosures; exact existential proof constants need no

&#x20; invented numerical N0. Avoid labeling a theorem "machine checked".

\- For each executed task: REPORT.md, RESULTS.json, PROOF_OR_OPEN_LEMMA.md,

&#x20; NONCONCLUSIONS.md, REPLAY.md and INPUTS.sha256. RESULTS separates supplied,

&#x20; derived and fitted objects, process/readout scope, errors and dependencies.

\- Read full primary theorem statements before applying them. Whitt's UI

&#x20; martingale-stability result is not a deterministic-clock FCLT for random stops.

\- No new physical inputs, experiments, purchases, messaging, FAR/selector,

&#x20; FIN-M0/B4/task335 or QM/GR/SM/ToE activation.



**## 7. Kill tests and scientific potential**



Stop a candidate if its actual reversible mobility loses tangent positivity,

its fast chart develops another soft direction, its readout correction survives,

or MET regularity remains unpaid. A coordinate-only law cannot silently be

presented as a law of the chosen physical readout. If all canonical effects

are common to the controls, accept class-level universality and reject

FIN-specificity of those effects.



Success would identify a collective statistical structure robust to some

microscopic kinetic freedom. It would not derive nature's A7, the actual

microscopic kinetics, local spacetime, energy units or measured reality.

Those remain separate bridge questions; neither a positive nor negative KI

result licenses adding their answers by hand.
