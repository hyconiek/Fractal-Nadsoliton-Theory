# Exact-kernel numerical certificate

The input values are **exact decimals**; `pi` is the mathematical circle constant enclosed by mpmath.iv. The k=3,4,5,6 Fourier coefficients of the finite q12 kernel are evaluated with interval functions, then A7 is built by the exact circulant Fourier sum. All 1728 logweights are outward evaluated. The maximal exact discrepancy between each interval endpoint and the binary64 producer value (regarded as an exact dyadic rational) is ε. Computing the third Hoeffding component using **exact rational arithmetic** on those dyadics yields ω0. Since ||(I-E1)(I-E2)(I-E3)||∞ ≤ 8, true residual oscillation belongs to [ω0−16ε,ω0+16ε]. The ensuing tanh(ω/4) upper bound uses outward interval tanh. This is a mathematically validated bound conditional on the mpmath.iv implementation (not a formal Lean/Coq proof).

{
  "execution_status": "EXECUTED",
  "scientific_verdict": "VALIDATED_INTERVAL_ENCLOSURE_CONDITIONAL_ON_MPMATH_IV_BACKEND",
  "mathematical_input": "exact decimal omega=0.18575 phi=0.1625 eta=1.8 g=5.145228719489142, exact pi",
  "source_certifier": "independent circulant Fourier A7 calculation with outward interval arithmetic",
  "source_A7_max_discrepancy_from_binary64": "1.79427020874074039359010566465266854506180183e-15",
  "F_exact_max_discrepancy_from_binary64_upper": "7.63461396215704738950128575409210815822378565e-15",
  "F_exact_max_interval_width": "6.5983878780401080014094494032802208883250173e-51",
  "nominal_exact_dyadic_ANOVA_osc": "0.205388295447682482179156219596431280175844828",
  "validated_residual_osc_lower": "0.205388295447560328355761706838199259603779355",
  "validated_residual_osc_upper": "0.205388295447804636002550732354663300747910302",
  "validated_TV_bound_upper": "0.05130199550486165514201468836168008105861",
  "third_contrast_finite_numerical": null,
  "method": "1728 outward interval logweights + rational (binary64-as-exact) ANOVA + projection norm bound ||Pi3||\u221e\u22648; no uncertified floating argmin used",
  "runtime_seconds": 0.5488784313201904,
  "limitations": [
    "mpmath.iv interval arithmetic implementation is assumed sound, not formally verified",
    "certifies the exact decimal mathematical tuple, not a separate implementation whose gain is interpreted as an IEEE binary64 constant",
    "physical source for kernel is outside this proof"
  ]
}
