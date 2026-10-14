**# FIN: next chain after the RG-C06 repository checkpoint**



Date: 2026-10-14

Status: RESEARCH_PLAN; this campaign has not been executed.

Checkpoint: fin_generative_review_20261014/RG_C06_DECISION.md.

Selected branch: RETAIN_EXPLICIT_MEMORY.

Execution directory: fin_generative_memory_campaign_20261014/.

Scope: RG-D01–RG-D05, followed by the architect-only RG-D06.



**## 1. What the previous campaign settled**



The inspected source package and repository intake support:



\- a non-pairwise static marginal for four-copy FIN;

\- a validated canonical pair-truncation bound TV<=0.051301995505;

\- an exact-decimal local-rate witness excluding strong lumpability of the three-copy first-two-label observation;

\- numerical finite-time deviations and a material dependence on hidden preparation;

\- the same qualitative mechanisms in trace-matched full-Potts and flat-P7 controls.



The earlier static event comparison is an internal consistency check.

The dynamic finite-window prediction task remains incomplete.



Do not redo the full static campaign or search for another unprotected two-label autonomous Markov model. Derive the additional information from the original process.



**## 2. The question**



Can a small memory representation, with coefficients and preparation maps derived from the three-copy microscopic generator, predict observed finite-window laws more accurately than its autonomous two-label projection?



“Small” is an operational compression objective, not a claim about the fundamental number of physical degrees of freedom. The complete microscopic chain is always the reference.



**## 3. Frozen parent and model scope**



Use:

\- q=12;

\- three labelled copies, 1728 states;

\- exact-decimal source kernel/gain convention accepted in the intake;

\- zero production control field;

\- leave-one-out heat-bath at attempt rate one per copy;

\- observation Y=(sigma1,sigma2), 144 values;

\- the exact two-label marginal of this same three-copy parent.



The earlier four-to-three static distribution is a separate example.

It must not be used as the stationary law of this three-to-two process.



Gain is the already-used g=5.145228719489142; g=0 is the control.

No task335 gain is opened, and no protocol333 parameter is changed.



The claim target is a mathematically controlled finite model.

No new physical source, field, particle, spatial interpretation or clock is asserted.



**## 4. Memory notation and construction**



Let L be the microscopic backward generator and pi its equilibrium law.

In the reversible representation,



S = diag(sqrt(pi)) L diag(1/sqrt(pi))



is symmetric negative semidefinite.



For observed labels y, define the orthonormal resolved embedding



U\_(x,y) = sqrt(pi_x / pi_Y(y)) \* 1{Y(x)=y}.



Then U^T U=I. Set Q=I-U U^T and define



A = U^T S U,

C = Q S U,

H = -Q S Q restricted to Ran(Q).



On resolved and hidden coordinates the microscopic system is



a' = A a + C^T b,

b' = C a - H b.



Eliminating b gives



a'(t) = A a(t)

&#x20;      \+ C^T exp(-tH) b(0)

&#x20;      \+ integral_0^t C^T exp[-(t-s)H] C a(s) ds.



The inhomogeneous preparation term is part of the model.

A stationary memory kernel alone does not describe a delta-prepared hidden label.



For K(t)=C^T exp(-tH) C, keep distinct:

\- K(0), instantaneous kernel strength;

\- K'(0), initial derivative;

\- M0=int_0^infinity K(t)dt;

\- M1=int_0^infinity t K(t)dt.



The old package's M0_coupled and M1_derivative fields are K(0) and K'(0).

They do not supply integrated moments or a tail certificate.



An important local consistency identity is



U^T S^2 U - A^2 = C^T C.



It quantifies the second-order semigroup defect. Checking it is a useful control; equal stationary laws or equal first derivatives cannot replace it.



**## 5. Reference arithmetic and propagation**



For the three-copy process, L=3(P-I), where P selects one copy uniformly and refreshes its label. This gives a stochastic uniformization kernel without fitting an event clock.



For probability rows:



p(t) = exp(-3t) sum\_(m>=0) (3t)^m/m! \* p(0) P^m.



Build P from interval-enclosed categorical conditional probabilities or another audited input/propagation scheme. A numerical SciPy Poisson tail is a useful diagnostic, not an all-operation certificate.



For z=3t and cutoff m with m+2>z, one usable tail bound is



Tail_m <= exp(-z) \* z^(m+1)/(m+1)! / [1-z/(m+2)].



Evaluate the bound outward. Selected observation events can be propagated without forming a dense 1728x1728 exponential.



Required reference allowance for the frozen test: <=1e-6 TV, or an explicitly revised allowance frozen before opening validation. This is a target to demonstrate, not a granted error bar.



**## 6. Dependencies and exact inputs**



\| Alias | File |

\|---|---|

\| I01 | AGENTS.md |

\| I02 | fin_generative_review_20261014/INTAKE_20261014.md |

\| I03 | fin_generative_review_20261014/verification.json |

\| I04 | fin_generative_review_20261014/verify_intake.py |

\| I05 | fin_generative_review_20261014/RG_C06_DECISION.md |

\| I06 | fin_generative_campaign_20261014/RG-C01/MODEL_CONTRACT.json |

\| I07 | fin_generative_campaign_20261014/RG-C03/PROCESS_MAP.json |

\| I08 | fin_generative_campaign_20261014/RG-C03/MEMORY_REPORT.json |

\| I09 | fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py |

\| I10 | fin_physical_bridge_campaign_01/PHYS-004/COUNTERMODELS.json |

\| I11 | fin_physical_bridge_campaign_01/PHYS-004/CALIBRATION_CONTRACT.json |

\| I12 | fin_generative_campaign_20261014/interval_certify.py |

\| I13 | FIN_GENERATIVE_NEXT_CAMPAIGN_20261014.md |



Source artifacts remain immutable. Use their actual file hashes.

A changed input invalidates any assumption of byte-identical replay.



**## 7. Task cards**



**### RG-D01 — One-parent static/dynamic contract and corrected quantities**



Class [A], priority P0, computational class S.

Question: are the effective stationary law, kinetics and preparation map bound to the same microscopic parent?

Inputs I01–I09, I12–I13.

Dependencies: the accepted RG-C06 checkpoint.



Steps:

1\. Freeze the parent and the Y map from section3.

2\. Derive its exact stationary two-label marginal.

3\. Bind the averaged resolved generator and memory coordinates to that marginal.

4\. Export harmonized K(0)/K'(0) notation without rewriting the older files.

5\. Freeze initial distributions: conditional equilibrium at Y=(0,0), hidden delta0 and delta1.

6\. State the complete readout and the role of initial hidden coordinates.

7\. Record the g=0 exact independent-refresh control.

8\. Preserve numerical and interval input roles separately.



Deliverables: MODEL_CONTRACT.json, PREPARATIONS.json, NOTATION_MAP.md.

PASS: all static and dynamic objects use the three-copy parent; no four-copy stationary law or integrated moment is silently substituted.

Kill-test: the description requires fitting an initial hidden state rather than deriving it from the declared preparation.

Forbidden shortcuts: renaming K(0) as M0; choosing a clock independently for the coarse generator; changing gain or task333.

After PASS: RG-D02 and RG-D03.

After FAIL: repair the contract before propagation.

Must not conclude: physical source or autonomous closed-system ontology.



**### RG-D02 — Independent finite-window reference with paid errors**



Class [A]/[B], priority P0, computational class M capped at 1728 states.

Question: can the reference observed probabilities be enclosed sufficiently tightly?

Inputs I04, I09, I12 and RG-D01.

Dependencies: RG-D01 PASS.



Steps:

1\. Enclose the exact mathematical kernel and required conditional rates.

2\. Construct the labelled P and account for self-refresh events correctly.

3\. Verify normalization and constant preservation analytically and numerically.

4\. Implement Poisson uniformization or another documented rigorous propagation route.

5\. Check already-opened training times t=0.15,0.5,1,1.5 and source probabilities.

6\. Propagate all 144 observed outcomes and include arithmetic and series tails.

7\. Check both stationary conditional and the two delta-hidden preparations.

8\. Reserve t=0.35 and t=0.8 and hidden delta2 for RG-D05; do not generate their validation outputs yet.

9\. Where source/input errors remain, export an OPEN bound rather than a numerical tolerance labelled certified.



Deliverables: REFERENCE_CHECKER.py, INPUT_ENCLOSURES.json, TRAINING_REFERENCE.json, ERROR_CERTIFICATE.json.

PASS: demonstrated reference uncertainty within the frozen allowance and agreement with numerical controls.

Kill-test: the error budget contains only a Poisson tail while input/matrix arithmetic is unbounded.

Forbidden shortcuts: clipping negative entries and hiding the mass correction; treating a residual as a full error theorem.

After PASS: RG-D04 and validation support for RG-D05.

After FAIL: bound a smaller set of observations or retain numerical evidence explicitly; architectural review before widening cost.

Must not conclude: the numerical parent is independently confirmed physics.



**### RG-D03 — Derived memory reduction with explicit preparation**



Class [B], priority P0, computational class M.

Question: how many derived hidden response directions are needed for the frozen observation class?

Inputs RG-D01, I07–I09 and section4.

Dependencies: RG-D01 PASS; RG-D02 may proceed independently.



Steps:

1\. Build A,C,H and verify U^T S^2 U-A^2=C^T C.

2\. For each training preparation express forward density as z0=p0/sqrt(pi), then a0=U^T z0 and b0=Qz0.

3\. Construct a hidden subspace from C a0 and the specified b0 vectors using a declared block Krylov/Galerkin rule.

4\. Start with orders1,2,4; cap the candidate at 32 hidden basis vectors beyond the 144 resolved coordinates.

5\. Derive all coefficients from S and the subspace. Do not fit decay rates or amplitudes to the training trajectory.

6\. Include the reduced preparation forcing and its omitted-coordinate error.

7\. Verify dissipativity and the equilibrium constant mode in the actual representation.

8\. Test output normalization/positivity for the claimed probabilities. A symmetric dissipative linear system is not automatically a CTMC.

9\. Compare with the averaged Markov baseline using the same clock and preparations.

10\. If a scalar/channel response reduction succeeds but the full probability law does not, narrow the claim rather than clipping probabilities.



Deliverables: MEMORY_REDUCTION.json, BASIS_CERTIFICATE.md, TRAINING_PREDICTIONS.json, PROBABILITY_SCOPE.md.

PASS: a derived reduction whose scope, coefficients and omitted residual are explicit.

Kill-test: meeting the target requires a separately fitted memory pole or preparation for each test.

Forbidden shortcuts: memory auxiliaries=particles/fields; dissipativity=Markov positivity; stationary K(t) substitutes for b0.

After PASS or a scoped NO_COMPACT_REDUCTION: RG-D04.

After FAIL: retain explicit microscopic hidden variables or the exact memory equation.

Must not conclude: memory is unique to FIN or is physical inertia.



**### RG-D04 — Freeze the total prediction and its countermodels**



Class [A]/[B], priority P0, computational class S/M.

Question: does the candidate improve the observed-process law with a defensible bound?

Inputs RG-D02/03, I10–I11.

Dependencies: complete training scope and reference allowances.



Steps:

1\. Choose one candidate using only the stated training outputs.

2\. Freeze the number of hidden directions and all preparation/readout maps.

3\. Sum kernel/input, numerical reference, subspace, memory and initial-condition errors.

4\. Use FIN, full-Potts and flat-P7 with matched trace and the same observation contract.

5\. Distinguish better dynamic prediction from an A7-specific mechanism.

6\. Freeze a two-time observation in the reserved window. Include the effect of filtering after the first readout.

7\. Check that a projected observation update has its own controlled error; one-time correlations do not establish multi-time prediction.

8\. Freeze primary target: <=0.01 TV on the selected full observed law, with the reference allowance from RG-D02.

9\. If the primary preparation class cannot meet that target, report it before validation. A narrower alternative needs a new declared scope and architectural approval.

10\. Hash the complete predictor and validation protocol.



Deliverables: FROZEN_PREDICTOR.json, PREREGISTRATION.json, TOTAL_ERROR_BUDGET.md.

PASS: a fixed model, preparation class, times and bound, with all unsupported claims marked OPEN.

Kill-test: the bound is fitted from the validation errors, or an invalid approximate probability law is silently repaired.

Forbidden shortcuts: a good scalar correlation claimed as a complete path law; the same free parameters independently reselected for controls.

After PASS: RG-D05.

After FAIL: stop for RG-D06 and preserve the bounded negative result.

Must not conclude: generic improvement sources A7 or proves universality.



**### RG-D05 — Complete the deferred dynamic test**



Class [A], priority P0, computational class S/M.

Question: does the frozen memory-aware model predict the reserved finite-window law?

Inputs frozen RG-D04 artifacts and RG-D02 reference generator.

Dependencies: RG-D04 PASS with hashes written before opening the test.



Steps:

1\. Generate the reserved reference outputs at t=0.35 and t=0.8 and the hidden-delta2 preparation.

2\. Compare the same finite-window law and all retained outcomes.

3\. Compare the averaged Markov baseline and the derived memory candidate.

4\. Evaluate the predeclared uncertainty and 0.01-TV criterion.

5\. Distinguish genuine improvement, scope-limited improvement, invalid probability predictions and failure.

6\. Do not revise the candidate after opening these outputs under the same ID.

7\. Export the complete data and an independently runnable verifier.



Deliverables: VALIDATION.json, REPORT.md, RAW_MODEL_OUTPUTS and REPLAY.md.

PASS: the predeclared law and preparation class meet the target with paid errors.

Kill-test: the new preparation requires a new fitted b0, clock or memory term.

Forbidden shortcuts: an internal computational test presented as laboratory validation.

After either outcome: stop and return the handoff for RG-D06.



**### RG-D06 — Architectural audit and next physical bridge**



Class [S], priority P0, computational class S.

Question: did the model produce a controlled, compact observed process?

Inputs all D01–D05 outputs and current AGENTS.md.

Dependencies: scoped terminal verdicts or an explicit early kill-test.



Decide:

1\. Which result is mathematically certified, numerical or still conditional?

2\. Is the approximation a probability model, a restricted response model or explicit full memory?

3\. Does it improve over the memoryless baseline without fitted kinetic parameters?

4\. Which observables distinguish supplied A7 from controls?

5\. Is the next step a physical multi-time readout design, a source atom or a necessary model expansion?



Deliverable: one CHECKPOINT_VERDICT.md with a bounded next action.

PASS: a clear decision retaining failures and the physical-source distinction.

Kill-test: repeated calibration or new auxiliary variables are introduced at every observation without a controlled rule.

After verdict: scope the next campaign separately.

No automatic source, geometry, particle, quantum or gravity activation.



**## 8. Parallel work and what remains deferred**



The static three-body witness and static TV theorem need no repetition unless their input model changes.



Global B4 still requires its own arithmetic/coverage repair. It is not a dependency of this three-copy dynamic test.



FIN-M0 remains OPEN_INPUT. A different relational model may be declared only as a new frozen modelling hypothesis, not as replay of unavailable historical runs.



PHYS-014 real precalibration is an independent engineered-physics route. Its hardware kinetics must be identified before applying the leave-one-out dynamic predictions. Stationary agreement alone is insufficient.



Task337 and task335 retain their original preparation and freeze gates.



No static or dynamic comparison here independently explains nature's choice of A7. That source question remains separate and demands a new independent typed object.



**## 9. Copyable executor instruction**



\~\~\~text

Read AGENTS.md, fin_generative_review_20261014/INTAKE_20261014.md,

RG_C06_DECISION.md and FIN_GENERATIVE_MEMORY_NEXT_CHAIN_20261014.md.



Execute RG-D01–RG-D05 under fin_generative_memory_campaign_20261014/.

Use one three-labelled-copy parent, q=12, the existing gain and the

first-two-label observation. No larger parent or new-g test.



Correct the old memory notation in your derived files:

M0_coupled is K(0), M1_derivative is K'(0).

Do not rewrite the source campaign.



Derive any auxiliary memory directions and their coefficients from

the microscopic generator. Include preparation forcing explicitly.

Pay arithmetic, propagation and observation-filter errors.

A dissipative linear response model is not automatically a CTMC.



Reserve the declared validation times/preparation until RG-D04 is

fully frozen with hashes. Retain every outcome and all failures.

No post-test retuning under the same identifier.



Report execution_status and scientific_verdict separately.

Export HANDOFF.md, CLAIM_REGISTER.json, NONCONCLUSIONS.md,

source/certificate hashes and the minimal replay, then stop.

RG-D06 is reserved for the architect. Do not alter AGENTS.md,

open task335, buy hardware or claim source/QM/GR/ToE closure.

\~\~\~
