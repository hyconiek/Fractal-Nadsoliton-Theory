# Arithmetic and uncertainty boundary

- FIN_A7 static 4→3: all **1728** configurations, exact-decimal kernel `(0.18575,0.1625,1.8)` and gain `5.145228719489142`; outward mpmath.iv input/formula enclosures; exact rational ANOVA on dyadic producer logweights; residual oscillation bounded by interval-projection operator norm. **Validated TV upper:** `0.05130199550486165514201468836168008105861`. This is conditional on correct `mpmath.iv` implementation and on **this** exact decimal model definition. The source physical A7 remains unproved.
- The true residual ω is inside [`0.205388295447560328355761706838199259603779355`, `0.205388295447804636002550732354663300747910302`]. Full extrema need not be assumed numerically exact; their error is enclosed by 16ε.
- The independent verifier passes 14/14 algebraic/model checks. Independent verifier does **not** constitute formal proof of scipy sparse propagation error.
- Controls' higher-body and TV results are exhaustive **binary64** finite enumerations, not outward-certified against separately defined transcendental control matrices.
- Dynamic CTMC: finite 1728-state sparse generator; reversibility, stationarity and row sums checked numerically. Krylov projected distributions and memory K(t) are numerically propagated without a formal adaptive error envelope.
- RG-C05: static computational holdout uses a producer-based prediction plus a differently enumerated 12^4 parent reference; both depend on the **same** model specification. The rigorous static input uncertainty is covered for FIN only in RG-C02, not an independent data-source error.
- No metrological calibration, device readout, physical time, universal statistical guarantee or external scientific prediction is inferred.
