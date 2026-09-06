# FIN: Thirty-step research ledger

## 01. ST8549 — Can both raw kernels be Markov conductances?

Strict is positive; the declared raw Legacy* profile is signed and fails the rate test.

Status: Proven (conditional; see proof)

Next-step rationale: Instead of taking absolute values silently, retain signs in a doubled state space.

Numerical/check evidence:

```json
{
  "strict_weights": [
    0.4699856726450201,
    0.1920435516901028,
    0.09142861427792495,
    0.0470291687456504,
    0.02413122336363006,
    0.011070817321442113
  ],
  "legacy_weights": [
    0.7104938272793256,
    -1.3591121187449902,
    -2.6001117014458375,
    -2.3087810266402764,
    -0.6834273957639219,
    1.307824868981029
  ],
  "negative_legacy": [
    [
      2,
      -1.3591121187449902
    ],
    [
      3,
      -2.6001117014458375
    ],
    [
      4,
      -2.3087810266402764
    ],
    [
      5,
      -0.6834273957639219
    ]
  ],
  "strict_row_sum": 1.6603072787660986,
  "legacy_raw_min_eigenvalue": -21.377059810423745
}
```

## 02. ST8550 — Can legacy signs be retained in a positive stochastic cover?

A 24-state positive double cover has unsigned-even and signed-odd Laplacian sectors.

Status: Proven (conditional; see proof)

Next-step rationale: A shifted signed amplitude is not a probability: test whether a positive gauge removes the distinction.

Numerical/check evidence:

```json
{
  "intertwining_errors": [
    6.153480596427404e-15,
    6.153480596427404e-15
  ],
  "odd_shift": [
    27.805728970380105,
    27.805728970380105,
    27.805728970380102,
    27.805728970380102,
    27.805728970380102,
    27.805728970380102,
    27.805728970380102,
    27.805728970380102,
    27.80572897038011,
    27.805728970380105,
    27.805728970380105,
    27.805728970380105
  ],
  "odd_min_eigenvalue": 6.428669159956358
}
```

## 03. ST8551 — Can a diagonal sign or positive Doob gauge make raw legacy positive?

A negative triangle product obstructs sign switching to all-positive weights; positive diagonal similarity preserves each sign.

Status: Proven (conditional; see proof)

Next-step rationale: The sign degree of freedom must remain explicit; determine what a coarse observer loses.

Numerical/check evidence:

```json
{
  "frustrated_triangles": 112,
  "witness": [
    0,
    1,
    2
  ],
  "witness_product": -0.686081807128401
}
```

## 04. ST8552 — Does the sign-cover automatically complete legacy into strict?

Coarse mass is blind to every sign assignment, but its weights are |legacy|, not strict. The construction is not a completion map.

Status: Proven (conditional; see proof)

Next-step rationale: A mathematical dilation is insufficient: test the claimed dual dynamics at the level of vertex records.

Numerical/check evidence:

```json
{
  "coarse_error": 8.702335715267317e-15,
  "different_cover_norm": 37.187493956634206,
  "best_absolute_legacy_scale": 0.061208206344650525,
  "relative_strict_mismatch": 0.8834120392397349
}
```

## 05. ST8553 — Are heat and unitary vertex records the same process?

Off-diagonal heat is t W_ij+O(t^2); a finite Hamiltonian gives t^2 |H_ij|^2+O(t^3). No positive linear clock identifies them.

Status: Proven (conditional; see proof)

Next-step rationale: Try an explicit repeated-measurement limit instead of Wick continuation.

Numerical/check evidence:

```json
{
  "short_time": [
    {
      "t": 0.01,
      "heat_over_t": 0.46337267843756896,
      "born_over_t2": 0.2208789560157308
    },
    {
      "t": 0.001,
      "heat_over_t": 0.4693194231361113,
      "born_over_t2": 0.22088645672568305
    },
    {
      "t": 0.0001,
      "heat_over_t": 0.469918997900532,
      "born_over_t2": 0.22088653173393274
    }
  ],
  "W01": 0.4699856726450201,
  "W01_squared": 0.220886532491592
}
```

## 06. ST8554 — Does rapid projective measurement recover the strict heat operator?

With H scaled as delta^(-1/2), repeated vertex measurements converge to rates |H_ij|^2. For H=A these are W_ij^2, not W_ij.

Status: Proven (conditional; see proof)

Next-step rationale: Try to compile the target rates explicitly, then test uniqueness of the compiled Hamiltonian.

Numerical/check evidence:

```json
{
  "convergence": [
    {
      "steps": 20,
      "error": 0.003514604900640912
    },
    {
      "steps": 80,
      "error": 0.000873991756636173
    },
    {
      "steps": 320,
      "error": 0.00021820932243590355
    },
    {
      "steps": 1280,
      "error": 5.453433178545506e-05
    }
  ],
  "squared_to_original_edge_ratios": [
    0.4699856726450201,
    0.1920435516901028,
    0.09142861427792495,
    0.0470291687456504,
    0.02413122336363006,
    0.011070817321442113
  ]
}
```

## 07. ST8555 — Can one compile strict heat from a Hamiltonian, uniquely?

Off-diagonal H_ij=sqrt(W_ij) exp(i theta_ij) compiles W in that limit. Cycle phases and diagonals remain free and change coherent dynamics.

Status: Proven (conditional; see proof)

Next-step rationale: A compiled measurement limit introduces a Hamiltonian and instrument; test an exact quantum channel extension instead.

Numerical/check evidence:

```json
{
  "rate_error": 1.4015878649856313e-16,
  "phase_changed_spectral_difference": 1.1123686151331336,
  "triangle_fluxes": [
    0,
    3.141592653589793
  ]
}
```

## 08. ST8556 — Is there an exact completely positive heat embedding?

Rank-one jumps sqrt(W_ij)|j><i| give a trace-preserving Lindblad semigroup with exactly autonomous strict heat populations for all density matrices.

Status: Proven (conditional; see proof)

Next-step rationale: Exact embedding exists; test whether vertex records determine its coherence law.

Numerical/check evidence:

```json
{
  "population_error": 6.546139416680338e-17,
  "trace_derivative": 0.0
}
```

## 09. ST8557 — Do complete heat population records select a quantum completion?

Arbitrary nonnegative vertex dephasing changes coherences while preserving every autonomous population trajectory.

Status: Proven (conditional; see proof)

Next-step rationale: Test the desired original coherent generator H=A in the same extension.

Numerical/check evidence:

```json
{
  "population_difference": 0.0,
  "full_derivative_difference": 2.8352180759972483
}
```

## 10. ST8558 — Can H=A be added without altering exact autonomous heat populations?

Not in the declared rank-one-jump plus vertex-dephasing class: autonomous heat for all density matrices forces H to be vertex-diagonal.

Status: Proven (conditional; see proof)

Next-step rationale: Do not identify two semigroups by a density-state Wick rotation; test normalization and linearity explicitly.

Numerical/check evidence:

```json
{
  "coherent_population_witness": [
    -0.017042566983990495,
    0.11866752540287523,
    -0.029518591907667417,
    -0.013927722839304713,
    -0.007413401652082918,
    -0.004177277420779322,
    -0.002315147952574829,
    -0.0011354684432248324,
    -0.0026251310275000102,
    -0.005392841961055529,
    -0.010995851623281969,
    -0.024123523591413232
  ],
  "witness_norm": 0.1274976234146897
}
```

## 11. ST8559 — Does normalized imaginary-time filtering define a universal physical channel?

It is nonlinear in the density matrix except for scalar filters. Unnormalized P rho P is a trace-decreasing CP operation, not deterministic heat on states.

Status: Proven (conditional; see proof)

Next-step rationale: Both stochastic and quantum lifts need extra premises. Test whether population consistency removes the stochastic freedom.

Numerical/check evidence:

```json
{
  "two_level_affinity_defect": 0.5385283921883663
}
```

## 12. ST8560 — Does exchangeability plus exact projective consistency force independent walkers?

No. Common vertex transpositions give a consistent exchangeable family with exactly the same single-particle generator Q and synchronized multi-particle jumps.

Status: Proven (conditional; see proof)

Next-step rationale: The fully synchronous law is reducible. Try irreducibility, detailed balance and the full equilibrium law as possible selectors.

Numerical/check evidence:

```json
{
  "single_marginal_error": 4.945168153041043e-15,
  "different_two_particle_generator": 8.846942060507196
}
```

## 13. ST8561 — Do irreducibility, reversibility and a product-uniform equilibrium select the kinetics?

No. Every mixture 0<=theta<1 is irreducible and reversible with the identical product-uniform equilibrium and exact singleton heat.

Status: Proven (conditional; see proof)

Next-step rationale: Static Shannon entropy is therefore insufficient. Locate a dynamical observable separating these realizations.

Numerical/check evidence:

```json
{
  "theta": 0.35,
  "symmetry_error": 0.0,
  "gap": 0.754121154207076,
  "uniform_stationary_error": 4.145082815114946e-15
}
```

## 14. ST8562 — Which infinitesimal observation distinguishes the hidden common dynamics?

Pair coincidence / cross quadratic variation distinguishes the generators even when every single-particle trajectory and equilibrium law agrees.

Status: Proven (conditional; see proof)

Next-step rationale: Ask exactly which extra dynamical premise makes singleton rates determine the full generator.

Numerical/check evidence:

```json
{
  "independent_pair_rate": 0.0,
  "common_pair_rate": 0.16449498542575702,
  "predicted": 0.16449498542575702
}
```

## 15. ST8563 — What minimal local premise yields a unique population generator?

Autonomous singleton generator Q at every configuration plus no simultaneous coordinate changes forces G_N=sum_k Q^(k). Exchangeability is unnecessary.

Status: Proven (conditional; see proof)

Next-step rationale: Replace the qualitative no-simultaneous-jump premise by a measurable, nonnegative witness.

Numerical/check evidence:

```json
{
  "three_coordinate_replay_errors": [
    0.0,
    0.0,
    0.0
  ],
  "configurations": 27
}
```

## 16. ST8564 — Can absence of common jumps be expressed without cancellation?

The pair budget B(x)=sum_y g_xy binom(K(x,y),2) is nonnegative and is zero exactly when every jump changes at most one coordinate.

Status: Proven (conditional; see proof)

Next-step rationale: A zero witness gives an exact theorem; derive a stable version for small nonzero coincidence rates.

Numerical/check evidence:

```json
{
  "pair_budget": [
    0.41999999999999993,
    0.41999999999999993,
    0.41999999999999993,
    0.41999999999999993,
    0.41999999999999993,
    0.41999999999999993,
    0.41999999999999993,
    0.41999999999999993
  ],
  "independent_budget": [
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
  ]
}
```

## 17. ST8565 — Does approximate pair locality bound the full finite-time model error?

With identical autonomous singleton rates, sup_x B(x)<=epsilon implies sup_x TV(P_t(x,.),P_ind,t(x,.))<=min(1,2 epsilon t).

Status: Proven (conditional; see proof)

Next-step rationale: Challenge the bound and its assumptions using many exchangeable particles with small jumps.

Numerical/check evidence:

```json
{
  "generator_tv_norm": 0.41999999999999993,
  "epsilon": 0.41999999999999993,
  "time": 0.8,
  "measured_tv": 0.04512605898471805,
  "upper_bound": 0.6719999999999999
}
```

## 18. ST8566 — Do exchangeability and a deterministic large-population limit select shot noise?

No. Uniformly chosen pairs with rate theta W/(N-1), mixed with single jumps, have the same mean heat and O(1/N) empirical variance but an altered noise prefactor.

Status: Proven (conditional; see proof)

Next-step rationale: This family is not projectively consistent across N. Test whether moment data can certify the missing event structure directly.

Numerical/check evidence:

```json
{
  "pair_model": [
    {
      "N": 2,
      "mean_per_particle": 0.7,
      "variance_per_particle": 0.9449999999999998,
      "variance_of_fraction": 0.4724999999999999
    },
    {
      "N": 4,
      "mean_per_particle": 0.7,
      "variance_per_particle": 0.9449999999999998,
      "variance_of_fraction": 0.23624999999999996
    },
    {
      "N": 16,
      "mean_per_particle": 0.7,
      "variance_per_particle": 0.9449999999999998,
      "variance_of_fraction": 0.05906249999999999
    },
    {
      "N": 64,
      "mean_per_particle": 0.7,
      "variance_per_particle": 0.9449999999999998,
      "variance_of_fraction": 0.014765624999999998
    }
  ],
  "variance_factor": 1.35
}
```

## 19. ST8567 — Which raw-event inequalities survive arbitrary nonnegative jump rates?

Local integer-count cumulants obey kappa2>=|kappa1| and kappa4>=kappa2; equality in the latter excludes all marks with |k|>=2.

Status: Proven (conditional; see proof)

Next-step rationale: These are instantaneous raw-count statements. Attempt to falsify their use with ordinary detector binning.

Numerical/check evidence:

```json
{
  "compound_mean": 1.0,
  "compound_variance": 2.0,
  "compound_fourth": 11.6,
  "strict_quadratic_variance_over_drift": 0.8048592087636892
}
```

## 20. ST8568 — Can an ordinary detector make valid Poisson events appear to violate the local count bound?

Yes. A one-click-per-bin threshold maps Poisson arrivals to Bernoulli clicks, with variance/mean=exp(-lambda Delta)<1. Finite-bin data cannot be inserted into the instantaneous theorem.

Status: Proven (conditional; see proof)

Next-step rationale: Separate removable loss from irreversible thresholding and test whether a robust coincidence protocol exists.

Numerical/check evidence:

```json
{
  "intensity": 2.0,
  "bin_width": 0.4,
  "click_probability": 0.5506710358827784,
  "click_variance_over_mean": 0.44932896411722156
}
```

## 21. ST8569 — Which detector imperfection preserves the zero-coincidence theorem?

Independent resolution-preserving Bernoulli thinning gives B_observed=eta^2 B. For eta>0 it preserves the exact zero set; a quantitative bound needs a lower efficiency bound.

Status: Proven (conditional; see proof)

Next-step rationale: Finite bins mix multiple events. Bound that contamination instead of assuming it vanishes.

Numerical/check evidence:

```json
{
  "thinning": [
    {
      "mark": 1,
      "detected_pair_mean": 0.0,
      "prediction": 0.0
    },
    {
      "mark": 2,
      "detected_pair_mean": 0.36,
      "prediction": 0.36
    },
    {
      "mark": 3,
      "detected_pair_mean": 1.08,
      "prediction": 1.08
    },
    {
      "mark": 4,
      "detected_pair_mean": 2.16,
      "prediction": 2.16
    },
    {
      "mark": 5,
      "detected_pair_mean": 3.6,
      "prediction": 3.5999999999999996
    },
    {
      "mark": 6,
      "detected_pair_mean": 5.3999999999999995,
      "prediction": 5.3999999999999995
    },
    {
      "mark": 7,
      "detected_pair_mean": 7.56,
      "prediction": 7.56
    },
    {
      "mark": 8,
      "detected_pair_mean": 10.08,
      "prediction": 10.08
    }
  ]
}
```

## 22. ST8570 — Can finite time resolution still test absence of simultaneous changes?

With a finite escape bound Lambda, a single-change null has Pr(at least two detected changes)<= (Lambda Delta)^2/2. A resolved joint event gives a positive order-Delta lower bound.

Status: Proven (conditional; see proof)

Next-step rationale: The rate orders separate analytically; compute a finite-sample test without synthetic-data fitting.

Numerical/check evidence:

```json
{
  "delta": 0.01,
  "escape_bound": 1.4,
  "null_upper_bound": 9.799999999999998e-05,
  "alternative_probability": 0.0005310458821762946,
  "alternative_lower_bound": 0.0004969931623084825
}
```

## 23. ST8571 — Does the proposed rate-order test have a finite resource specification?

For independently reset trials a predeclared binomial upper-tail test controls the null and has high calculated power in the supplied two-particle model. No physical trial was performed.

Status: Proven conditional test rule; numerical tails supplemented by exact rational bounds

Next-step rationale: Attempt to destroy identifiability with an uncalibrated but nondegenerate detector.

Numerical/check evidence:

```json
{
  "trials": 200000,
  "threshold": 40,
  "numerical_minimum_threshold": 36,
  "nominal_alpha": 0.001,
  "calculated_size": 3.472385031549474e-05,
  "calculated_power": 0.9999999999999405,
  "null_probability_bound": 9.799999999999998e-05,
  "alternative_probability": 0.0005310458821762946,
  "exact_size_upper_bound": "672749994932560009201/1208925819614629174706176",
  "exact_miss_upper_bound": "239299329230617529590083/324518553658426726783156020576256",
  "size_upper_bound_float": 0.0005564857529033613,
  "miss_upper_bound_float": 7.373979901392417e-10
}
```

## 24. ST8572 — Can unknown event merging exactly hide common dynamics?

Yes. One-click-per-event reporting of the common-flip mixture and independently thinned single flips can produce the same entire nonzero Poisson record law despite different microscopic generators.

Status: Proven (conditional; see proof)

Next-step rationale: Even a mark-preserving detector may leave a clock/efficiency ambiguity. Test the full record likelihood, not just its mean.

Numerical/check evidence:

```json
{
  "theta": 0.3,
  "null_efficiency": 0.85,
  "common_observed_rate": 1.19,
  "independent_observed_rate": 1.19
}
```

## 25. ST8573 — Does full resolved-count likelihood determine both clock and efficiency?

No. The two-particle count law has only c1,c2; a simultaneous change of w,theta,eta preserves both and hence every counting distribution.

Status: Proven (conditional; see proof)

Next-step rationale: Operational nonuniqueness is explicit. Audit the earlier nonlinear/prism and inverse-action alternatives before the final synthesis.

Numerical/check evidence:

```json
{
  "original": {
    "w": 0.7,
    "theta": 0.2,
    "eta": 0.6,
    "c1": 0.7392,
    "c2": 0.05039999999999999
  },
  "indistinguishable": {
    "w": 1.4,
    "theta": 0.4,
    "eta": 0.3,
    "c1": 0.7392,
    "c2": 0.05039999999999999
  }
}
```

## 26. ST8574 — Can the failed pure-parity kinetic counterexample be repaired honestly?

Adding an explicitly declared symmetric field gives two flip-energy magnitudes, so Metropolis and heat-bath cease to differ only by a clock. The original pure-parity proportionality is retained.

Status: Proven (conditional; see proof)

Next-step rationale: Equilibrium geometry fails to select kinetics. Test the analogous inverse-propagator-to-action claim.

Numerical/check evidence:

```json
{
  "field": 0.17,
  "rate_ratios": [
    1.390627835359,
    1.771051585804
  ],
  "stationary_error": 6.359601310784502e-17,
  "old_proportionality_error": 5.874748045952207e-16
}
```

## 27. ST8575 — Does a supplied propagator determine the nonlinear variational theory?

No. Every S_lambda=phi^T M phi/2-J^T phi+lambda sum phi_i^4, lambda>=0, has the identical zero-source Hessian and linear response M^-1 but different cubic response.

Status: Proven (conditional; see proof)

Next-step rationale: Test whether hidden-state compression can select the missing dynamical response from the same static Green data.

Numerical/check evidence:

```json
{
  "coupling": 0.2,
  "cubic_checks": [
    {
      "source_scale": 0.1,
      "cubic_expansion_error": 3.7689956957822254e-05
    },
    {
      "source_scale": 0.05,
      "cubic_expansion_error": 1.204853567991613e-06
    },
    {
      "source_scale": 0.025,
      "cubic_expansion_error": 3.786994696757599e-08
    }
  ]
}
```

## 28. ST8576 — Does exact static Green/Schur compression select a hidden memory law?

No. Positive two-variable actions have identical coarse static Green function but different rational response and exponential memory; iteration does not remove the unsourced hidden parameters.

Status: Proven (conditional; see proof)

Next-step rationale: The surviving ambiguities occur in signs, quantum instruments, collective jumps, interactions and memory. Formulate only the no-go actually supported by explicit models.

Numerical/check evidence:

```json
{
  "positive_completions": [
    {
      "c": 1.0,
      "d": 1.0,
      "static_Green": 0.5,
      "effective_slope": 2.0,
      "memory_amplitude": 1.0,
      "memory_decay": 1.0
    },
    {
      "c": 2.0,
      "d": 3.0,
      "static_Green": 0.5000000000000001,
      "effective_slope": 1.4444444444444444,
      "memory_amplitude": 4.0,
      "memory_decay": 3.0
    }
  ]
}
```

## 29. ST8577 — Can the audited finite FIN data logically determine a unique collective physical realization?

Not in the stated finite Markov completion class: same W, all singleton paths, product equilibrium, detailed balance, projective consistency and FIN graph symmetries coexist with different pair records.

Status: Proven (conditional; see proof)

Next-step rationale: This is not a universal impossibility of future physics. Close with a premise-removal and robustness audit of the smallest demonstrated positive bridge.

Numerical/check evidence:

```json
{
  "time": 0.02,
  "independent_pair_probability": 8.349151200901338e-05,
  "common_pair_probability": 0.0031556560516842556
}
```

## 30. ST8578 — Which rigorously narrowed bridge survives premise removal and perturbations?

In a declared finite CTMC, autonomous singleton Q plus zero pair-coincidence budget uniquely determines the generator; a small budget gives a stable TV bound. Initial law, clock, detector and the source of zero budget remain separate obligations.

Status: Proven (conditional; see proof)

Next-step rationale: Next: seek a strict-derived composition law forcing the nonnegative pair budget to vanish, or prove its independence in a precisely declared enlarged source class.

Numerical/check evidence:

```json
{
  "random_graph_stress": [
    {
      "marginal_error": 2.0014830212433605e-16,
      "tv": 0.03438957099236574,
      "bound": 0.04978534187136315
    },
    {
      "marginal_error": 5.438959822042073e-16,
      "tv": 0.07825897581081953,
      "bound": 0.23104234502539595
    },
    {
      "marginal_error": 4.440892098500626e-16,
      "tv": 0.2527544268721992,
      "bound": 0.5285344692696886
    },
    {
      "marginal_error": 3.6821932062951477e-16,
      "tv": 0.05373237389252959,
      "bound": 0.10199851454531177
    },
    {
      "marginal_error": 3.8459253727671276e-16,
      "tv": 0.18789112116294615,
      "bound": 0.34339612680354376
    },
    {
      "marginal_error": 7.021666937153402e-16,
      "tv": 0.21544032247319858,
      "bound": 0.5318662276654861
    },
    {
      "marginal_error": 4.710277376051325e-16,
      "tv": 0.002273559608273544,
      "bound": 0.006677513130809377
    },
    {
      "marginal_error": 2.482534153247273e-16,
      "tv": 0.17638821998614002,
      "bound": 0.36683171271401843
    },
    {
      "marginal_error": 4.965068306494546e-16,
      "tv": 0.32778065542129764,
      "bound": 0.6917855541322852
    },
    {
      "marginal_error": 8.455206652451151e-16,
      "tv": 0.44124951252296685,
      "bound": 0.9347138971710798
    },
    {
      "marginal_error": 5.20740757162067e-16,
      "tv": 0.3201008975437382,
      "bound": 0.5358318712882749
    },
    {
      "marginal_error": 3.1401849173675503e-16,
      "tv": 0.3647545195087196,
      "bound": 0.8122774272289731
    }
  ],
  "removed_no_joint_witness": "rounds 12--14",
  "removed_singleton_rate_witness": "2 times the independent generator",
  "remaining_source": "No FIN theorem presently forces the pair budget to vanish."
}
```
