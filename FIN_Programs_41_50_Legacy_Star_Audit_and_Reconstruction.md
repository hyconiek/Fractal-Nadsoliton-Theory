# FIN Programs 41–50 Legacy Star Audit and Reconstruction

## Metadata

- Date UTC: `2026-07-19T18:46:10.397627+00:00`
- Git HEAD: `fd31e59e5a9b703f837b54e87176117719ffcfbb`
- Python: `3.14.4`
- NumPy: `2.5.1`
- Seed convention: `20260719`

## Executive summary

This audit treats the Program 42a and Programs 41–50 artifacts as inputs to be checked, not as authorities. The independent reconstruction confirms the main algebraic correction while narrowing several claims.

- **[Refuted]** The double-damped product `K_rej(d)=exp(-2.9d)*(1+0.2d)/(1+d)*cos(pi*d/4+pi/6)` is not a faithful historical reconstruction.
- **[Proven]** Under the historical role assignment and the path-sum correction, the admissible reconstructed class is `K*_legacy(d)=A*cos(pi*d/4+pi/6)/(1+beta*d)`.
- **[No-go on current artifacts]** The data do not derive a unique tuple `(A,beta)`; `A=2.9`, `A=4 ln 2`, `beta=0.01`, `beta=0.05`, and `beta=1` are freezes, gauges, or comparisons unless an extra source theorem is supplied.
- **[Strong evidence]** Absolute-value/positivity-repaired `K*` supports useful finite dual-dynamics operators on `C12`; the raw signed object must not be silently called a Markov generator.
- **[No-go on current artifacts]** `K*` does not close physical units, `QW-2191`, orientation/polarity, legacy-to-strict bridge completion, role transfer, `L_total`, or ToE.

## Repository guardrails and scope

The report keeps the legacy/strict split explicit: `K_legacy_ont(d)=alpha_geo*cos(omega*d+phi)/(1+beta_tors*d)` is an intermediate bridge kernel, while `K_strict_gate(d)=cos(omega*d+phi)/(1+beta*d^eta)` is a later operational strict working kernel. No silent identity or role transfer is used. The ontology remains `nadsoliton -> light -> matter -> emergent observer`, and scalar Shannon identities do not discharge `QW-2191`.

## Program 42a independent reconstruction

Starting axiom:

```text
K_total = K_geo*K_res*(1+0.2*K_tors)*K_topo
```

| Factor | Audited historical role |
| --- | --- |
| K_geo | exp(-2.9 d) as historical viscosity before transformation |
| K_res | approximately constant phase-sync factor 0.8-1.2 |
| K_tors | cos(pi d/4 + pi/6) turbulent oscillatory current |
| K_topo | path-sum/topological factor transformed to hyperbolic 1/(1+beta d) |

### Algebraic errors A/B/C

- **Error A:** `2.9` is retained as the historical exponential rate when the pre-transform `K_geo` is written. The fatal error is multiplying that exponential by a second hyperbolic damping factor after the path-sum transform.
- **Error B:** If `N(d) ~ d^1.6` and the target total tail is `d^-1`, then a single-path amplitude must scale as `d^-2.6`; the earlier `d^-0.6` gives growth `d^1.0`, not decay.
- **Error C:** `cos(pi*d/4+pi/6)=0` exactly at `d=4/3+4n`; the integer sequence `2,5,8,11,...` is not a zero sequence.

### Accepted class and parameter status

The reconstructed class is:

```text
K*_legacy(d)=A*cos(pi*d/4+pi/6)/(1+beta*d)
```

| Quantity | Status | Comment |
| --- | --- | --- |
| omega=pi/4 | derived/frozen from historical phase | Do not replace with strict omega. |
| phi=pi/6 | derived/frozen from historical phase | Do not replace with strict phi. |
| hyperbolic envelope | derived as transformed path-sum class | Not independent extra damping. |
| A=2.9 | historical diagram-scale freeze | Not unique theorem. |
| A=4 ln 2 | informational normalization candidate | Exact Shannon identity, not strict amplitude source. |
| beta=0.01 | legacy canonical freeze/post-hoc classical comparison | Not derived here. |
| beta=0.05 | historical working calibration | Not strict source. |
| beta=1 | Z12/unitless Laplacian gauge/comparison | Not legacy historical derivation. |

## Product candidate rejection

The independent profile check on `d=1..12` finds `|K_rej|<1e-3` for all `d>=2`: `True`. Its correlation with the `A=4 ln 2, beta=0.01` classical profile is `0.209715`. This is profile-level strong evidence in addition to the algebraic refutation.

| Reason | Verdict |
| --- | --- |
| role swap: cosine belongs to torsion/current, not K_res | Refuted product assumption |
| K_tors=d is not the historical oscillatory current | Refuted product assumption |
| exp(-2.9d) multiplied by 1/(1+d) double-counts attenuation after transform | Refuted product assumption |
| (1+0.2d)/(1+d) tends to 0.2 rather than generating a 1/d path-sum tail by itself | Refuted product assumption |
| profile is killed exponentially by d=2 on d=1..12 | Refuted product assumption |

## Relation to Programs 31–40

Programs 31–40 remain compatible with this audit only in the conditional sense: negative information coupling can be a bridge component, but it is not hidden in the static legacy formula and cannot complete the legacy-to-strict bridge alone. The corrected Program 42a object improves the intermediate kernel but still lacks phase/frequency transformation, nonlinear damping/compression source, state/time/channel interpretation, information functional, unit bridge, and sign reference.

## `A = 4 ln 2` information-amplitude audit

`4 ln 2 = ln 16 = 2.772588722239781` is **[Proven]** as the Shannon entropy of a uniform 16-state/four-bit source: `exp(4 ln 2)=15.999999999999998` and `bits=4.0`. It is **[Speculative]** as a unique physical amplitude law without an additional source theorem. It does not export SI units, hbar, non-premise selector, orientation polarity, role transfer, or ToE.

## Operator and dual-dynamics audit

The finite operator check constructs `W_ij=K(|i-j|)` on `C12` and `L=diag(W 1)-W`. For signed kernels this is a signed Laplacian candidate; for absolute-value repaired kernels it is a nonnegative-weight graph Laplacian. Only the repaired kernels should be used for Markov-cone claims.

| Kernel | min eig | gap | neg weights | PSD tol |
| --- | --- | --- | --- | --- |
| K_rej_signed | 0.000e+00 | 0.00144564 | 96 | True |
| K_rej_abs | 0.000e+00 | 0.00313343 | 0 | True |
| Kstar_hist_abs | 0.000e+00 | 14.1791 | 0 | True |
| Kstar_info_abs | 0.000e+00 | 15.1854 | 0 | True |
| Kstar_Z12_abs | 0.000e+00 | 1.50776 | 0 | True |
| Kstrict_abs | 0.000e+00 | 0.754121 | 0 | True |

## Programs 41–50 detailed audit

| Program | Independent audit verdict | Key numeric/status from JSON |
| --- | --- | --- |
| Program 41 | Support/normalization dependence confirmed; C12-cyclic PSD row cannot be promoted to a support-independent bridge theorem. | json_theorem_scope=A_s - A_ℓ ≽ 0 is support- and normalization-dependent; not a universal passive-dissipation bridge. |
| Program 42 | Affine phase transport can reduce carrier mismatch but leaves envelope residual and exports no source law. | affine_residual_l2=8.1773032027145e-16, exports_strict_phase_source=False |
| Program 43 | Held-out hazard behavior is not strong enough to infer a microscopic law or eta source. | train_relative_l2_d1_4=0.13956336349983897, holdout_relative_l2_d5_6=0.9988656436697538 |
| Program 44 | Markov relative information is monotone in the declared test; unitary Born information is not a Markov loss law. | markov_ck_residual=1.7826341206873999e-16, markov_backflow_steps=0, unitary_born_backflow_steps=91 |
| Program 45 | Apparent system loss is compatible with environment storage; not fundamental deletion. | environment_holds_at_least_lost_info=True, classical_record_recovery_fidelity=1.0 |
| Program 46 | Polarity separation remains apparatus/reference conditioned and does not discharge QW-2191. | exports_nonpremise_selector=False |
| Program 47 | Influence-functional proxy is receiver-like; no bath/source theorem from kernel alone. | log_retention_relative_residual=0.17708641220032095, exports_bath_from_kernel_alone=False |
| Program 48 | Stable feedback is not automatically variational; L_total is not exported. | stable=True, normalized_integrability_defect=0.9048187022009939 |
| Program 49 | Process-tensor challenge is useful synthetic discrimination, not experimental validation. | mean_intervention_margin=0.004619757690376036 |
| Program 50 | Multi-size selection favors strict if included; this is benchmark evidence, not a legacy-to-strict derivation. | winner_consistency=['strict', 'strict', 'strict'], strict_always_wins_if_present=True |

## Cross-document consistency audit

- The existing Program 42a JSON class/verdict agrees with the independent reconstruction at the class level.
- The Programs 41–50 JSON reports `used_strict_in_reconstruction=False` and `uniqueness=False`, consistent with this audit's non-unique-class verdict.
- The prior Markdown/TeX/PDF monograph is broadly aligned with the JSON verdicts, but this report makes the parameter-status and source-theorem limitations more explicit.
- The current environment used NumPy `2.5.1`; pre-existing JSON metadata may record an older NumPy. Numeric conclusions checked here are stable at the level used by the report.

## Leaf-cut analysis

| Claim | Status | Missing premise |
| --- | --- | --- |
| physical units | No-go on current artifacts: all kernels remain dimensionless; no unit-bearing bridge is exported. | unit-bearing bridge / conversion constants |
| selector / QW-2191 | No-go on current artifacts: radial K* is inversion-blind and Program 46 remains reference-conditioned. | non-premise selector source |
| symmetry breaking | No-go on current artifacts: no internal non-premise orientation/polarity source is exported. | internal orientation/polarity source and coupling theorem |
| legacy->strict bridge | Partial progress only: better intermediate object and benchmarks; missing phase/frequency source, nonlinear compression source, units, selector, and role-transfer theorem. | completion map plus role-transfer audit |

## Final verdict

| Object | Verdict |
| --- | --- |
| program42a | Accepted with corrections: the class K*_legacy is algebraically justified, but not a unique tuple. |
| K_rej | Refuted as historical reconstruction. |
| Kstar_class | Proven as corrected path-transformed class under historical assumptions; parameters remain free/frozen. |
| A_4ln2 | Proven Shannon identity and plausible informational normalization; not a strict source theorem. |
| programs41_50 | Reproducible finite/synthetic suite with generally guarded conclusions; no strict closure. |
| ToE_or_strict_core | No-go on current artifacts. |

## Recommended next research moves

No generic replay is unlocked. A legitimate next move must introduce exactly one genuinely new strict typed object/source/theorem/blocker-cut/provider class. In this lane, the most relevant admissible candidates would be: (i) a source theorem fixing `A` or `beta` without post-hoc target insertion; (ii) an explicit phase/frequency source localizer; (iii) a nonlinear compression-source theorem connecting linear `beta*d` to strict `beta*d^eta`; or (iv) an independent orientation/polarity source with coupling theorem. Without one of these, preserve the no-new-live-frontier/conditional-bridge status.

## Reproducibility appendix

Commands run for the audit session included:

```text
python fin_programs_41_50_legacy_star.py
python fin_kstar_info_4ln2_research.py
python fin_programs_41_50_legacy_star_audit.py
python -m py_compile fin_programs_41_50_legacy_star_audit.py
git status --short
```

Input file status and SHA256 hashes are stored in `FIN_Programs_41_50_Legacy_Star_Audit_and_Reconstruction_Results.json`.
